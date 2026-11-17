"""CZ drive sweep with Förster coupling and explicit Rydberg decay.

The data loss channels leave the computational space. Ancilla loss is retained
by integrating its jump trajectories, so the result is the reduced-data
Lindblad channel fidelity rather than a subtracted decay estimate.
"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import solve_ivp
from scipy.linalg import expm

from models import forster_operators


BASE = 2 * np.pi * 10
V = 2 * np.pi * 200
GAMMA_DATA = 1 / 100  # inverse microseconds
GAMMA_ANCILLA = 1 / 200
TARGET = np.diag([-1, -1, -1, 1])
AREA = np.trapezoid(
    np.exp(-0.5 * ((np.linspace(0, 1, 2001) - .5) / .25) ** 2) - np.exp(-2),
    np.linspace(0, 1, 2001),
)


def envelope(t, start, duration):
    x = (t - start) / duration
    if x < 0 or x > 1:
        return 0.0
    return (np.exp(-.5 * ((x - .5) / .25) ** 2) - np.exp(-2)) / AREA


def gate_point(mhz, shape="square", noisy=True, quadrature=32,
               rtol=2e-9, atol=2e-11,
               gamma_data=GAMMA_DATA, gamma_ancilla=GAMMA_ANCILLA):
    """Return (average fidelity, maximum computational leakage).

    ``shape='ancilla'`` smooths the only pulse during which Förster-coupled
    pair states are populated. Both shapes have the same duration and area.
    """
    omega = 2 * np.pi * mhz
    td, ta = np.pi / omega, 2 * np.pi / omega
    T = 2 * td + ta
    (_, _, _, _, _, _, _, _, data_q, ancilla_q, pair_q, inputs) = (
        forster_operators(BASE, BASE, V))
    hd = data_q.full() * (omega / BASE)
    ha = ancilla_q.full() * (omega / BASE)
    pair = pair_q.full()
    states0 = np.column_stack([state.full().ravel() for state in inputs])
    comp = np.array([0, 1, 4, 5])
    # Index order is ancilla (A,R,R'), data 1, data 2.
    data_excited = np.array([int(i in (2, 3)) for i in range(4)])
    data_decay = (data_excited[:, None] + data_excited[None, :]).ravel()
    decay = np.tile(data_decay, 3).astype(float) * gamma_data
    decay += gamma_ancilla
    heff_idle = pair - .5j * np.diag(decay if noisy else np.zeros(48))
    # After ancilla loss, pair interaction and ancilla drive vanish.
    data_idle = -.5j * gamma_data * np.diag(data_decay) if noisy else np.zeros((16, 16), complex)
    data_drive = hd[:16, :16] + data_idle
    stages = [(0, td, hd), (td, td + ta, ha), (td + ta, T, hd)]
    x = states0
    trajectories = []
    nodes, weights = np.polynomial.legendre.leggauss(quadrature)
    for stage, (start, end, drive) in enumerate(stages):
        sample_t = start + (nodes + 1) * (end - start) / 2
        if shape == "ancilla" and stage == 1:
            def rhs(t, y):
                h = heff_idle + envelope(t, start, end - start) * drive
                return (-1j * h @ y.reshape(48, 4)).ravel()
        else:
            h = heff_idle + drive
            def rhs(t, y):
                return (-1j * h @ y.reshape(48, 4)).ravel()
        sol = solve_ivp(rhs, (start, end), x.ravel(), t_eval=np.r_[sample_t, end],
                        rtol=rtol, atol=atol, method="DOP853")
        if not sol.success:
            raise RuntimeError(sol.message)
        trajectories.append((sample_t, weights * (end - start) / 2,
                             sol.y[:, :-1].T.reshape(quadrature, 3, 16, 4)))
        x = sol.y[:, -1].reshape(48, 4)

    def score(k):
        return (np.vdot(k, k).real + abs(np.trace(TARGET.conj().T @ k)) ** 2) / 20

    fidelity = sum(score(x[a * 16 + comp, :]) for a in range(3))
    if noisy:
        for stage, (sample_t, sample_w, states) in enumerate(trajectories):
            for t, w, psi in zip(sample_t, sample_w, states):
                if stage == 0:
                    post = expm(-1j * data_drive * td) @ expm(-1j * data_idle * ta) @ expm(-1j * data_drive * (td - t))
                elif stage == 1:
                    post = expm(-1j * data_drive * td) @ expm(-1j * data_idle * (td + ta - t))
                else:
                    post = expm(-1j * data_drive * (T - t))
                for a in range(3):
                    k = (post @ psi[a])[comp, :]
                    fidelity += gamma_ancilla * w * score(k)
    computational_return = np.sum(abs(x.reshape(3, 16, 4)[:, comp, :]) ** 2,
                                  axis=(0, 1))
    leakage = float(np.max(1 - computational_return))
    return float(fidelity), leakage


def main():
    drives = np.array([5, 7.5, 10, 12.5, 15, 20, 25, 27.5, 30,
                       32.5, 35, 37.5, 40])
    rows = []
    for mhz in drives:
        square_u, square_l = gate_point(mhz, noisy=False)
        smooth_u, smooth_l = gate_point(mhz, shape="ancilla", noisy=False)
        square_n, _ = gate_point(mhz)
        smooth_n, _ = gate_point(mhz, shape="ancilla")
        rows.append([mhz, square_u, smooth_u, square_n, smooth_n,
                     square_l, smooth_l])
        print(f"{mhz:g} MHz: square={square_n:.7f}, smooth={smooth_n:.7f}, "
              f"leakage={square_l:.5f}/{smooth_l:.5f}", flush=True)
    rows = np.array(rows)
    Path("results").mkdir(exist_ok=True)
    Path("figures").mkdir(exist_ok=True)
    np.savetxt("results/drive_sweep.csv", rows, delimiter=",", comments="",
               header="drive_mhz,square_unitary,ancilla_smooth_unitary,square_noisy,ancilla_smooth_noisy,square_max_leakage,ancilla_smooth_max_leakage")
    baseline = rows[np.where(drives == 10)[0][0]]
    data_only, _ = gate_point(10, gamma_ancilla=0)
    ancilla_only, _ = gate_point(10, gamma_data=0)
    np.savetxt("results/reproduction_summary.csv",
               [[baseline[1], data_only, ancilla_only, baseline[3]]],
               delimiter=",", comments="",
               header="forster_no_decay,forster_data_decay_only,forster_ancilla_decay_only,forster_both_decay")
    fig, axes = plt.subplots(2, 1, sharex=True, figsize=(7, 7))
    axes[0].plot(drives, 100 * rows[:, 3], "o-", label="Square")
    axes[0].plot(drives, 100 * rows[:, 4], "s-", label="Smooth ancilla")
    axes[0].set_ylabel("Average data-gate fidelity (%)")
    axes[0].legend()
    axes[1].semilogy(drives, rows[:, 5], "o-", label="Square")
    axes[1].semilogy(drives, rows[:, 6], "s-", label="Smooth ancilla")
    axes[1].set_ylabel("Maximum computational leakage")
    axes[1].set_xlabel("Data and ancilla Rabi frequency (MHz)")
    axes[1].legend()
    fig.tight_layout()
    fig.savefig("figures/drive_sweep_noisy.png", dpi=180)


if __name__ == "__main__":
    main()
