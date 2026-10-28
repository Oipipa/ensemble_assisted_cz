import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path
from qutip import basis, expect, qeye, sesolve, tensor


zero = basis(3, 0)
one = basis(3, 1)
r = basis(3, 2)
Omega = 2 * np.pi * 10  # radians per microsecond
H = Omega / 2 * (one * r.dag() + r * one.dag())
times = np.linspace(0, 2 * np.pi / Omega, 401)

result_zero = sesolve(H, zero, times)
result_one = sesolve(H, one, times)
p_zero = [abs(zero.overlap(state)) ** 2 for state in result_zero.states]
p_one = [abs(one.overlap(state)) ** 2 for state in result_one.states]
p_r = [abs(r.overlap(state)) ** 2 for state in result_one.states]

print(f"|0> unchanged: {p_zero[-1]:.6f}")
print(f"Pi-pulse transfer: {p_r[200]:.6f}")
print(f"2pi-pulse return amplitude: {one.overlap(result_one.states[-1]):.6f}")

plt.plot(times, p_zero, label=r"$P_0$ from $|0\rangle$")
plt.plot(times, p_one, label=r"$P_1$ from $|1\rangle$")
plt.plot(times, p_r, label=r"$P_r$ from $|1\rangle$")
plt.xlabel(r"Time ($\mu$s)")
plt.ylabel("Population")
plt.legend()
plt.tight_layout()

H_pair = tensor(H, qeye(3)) + tensor(qeye(3), H)
n_r = r * r.dag()
rydberg_count = tensor(n_r, qeye(3)) + tensor(qeye(3), n_r)
inputs = [("00", tensor(zero, zero)), ("01", tensor(zero, one)),
          ("10", tensor(one, zero)), ("11", tensor(one, one))]

plt.figure()
for label, initial in inputs:
    result = sesolve(H_pair, initial, times, e_ops=[rydberg_count])
    print(f"|{label}>: Rydberg excitations after pi pulse = {result.expect[0][200]:.6f}")
    plt.plot(times, result.expect[0], label=rf"$|{label}\rangle$")

plt.xlabel(r"Time ($\mu$s)")
plt.ylabel("Mean number of Rydberg excitations")
plt.legend()
plt.tight_layout()

A = basis(2, 0)
R = basis(2, 1)
Omega_R = 2 * np.pi * 10
H_ancilla = Omega_R / 2 * (A * R.dag() + R * A.dag())
ancilla_times = np.linspace(0, 2 * np.pi / Omega_R, 401)
ancilla_result = sesolve(H_ancilla, A, ancilla_times)
p_A = [abs(A.overlap(state)) ** 2 for state in ancilla_result.states]
p_R = [abs(R.overlap(state)) ** 2 for state in ancilla_result.states]
print(f"Ancilla pi transfer: {p_R[200]:.6f}")
print(f"Ancilla 2pi return amplitude: {A.overlap(ancilla_result.states[-1]):.6f}")

plt.figure()
plt.plot(ancilla_times, p_A, label=r"$P_A$")
plt.plot(ancilla_times, p_R, label=r"$P_R$")
plt.xlabel(r"Time ($\mu$s)")
plt.ylabel("Population")
plt.legend()
plt.tight_layout()

H_data_joint = tensor(qeye(2), H, qeye(3)) + tensor(qeye(2), qeye(3), H)
H_ancilla_joint = tensor(H_ancilla, qeye(3), qeye(3))
A10 = tensor(A, one, zero)
R10 = tensor(R, one, zero)
Ar0 = tensor(A, r, zero)
print(f"Data drive |A10> -> |Ar0>: {Ar0.overlap(H_data_joint * A10):.3f}")
print(f"Ancilla drive |A10> -> |R10>: {R10.overlap(H_ancilla_joint * A10):.3f}")
print(f"Cross action: {R10.overlap(H_data_joint * A10):.3f}, "
      f"{Ar0.overlap(H_ancilla_joint * A10):.3f}")

joint_result = sesolve(H_data_joint + H_ancilla_joint, A10, times)
joint_p_R = [abs(R10.overlap(state)) ** 2 for state in joint_result.states]
joint_p_Ar0 = [abs(Ar0.overlap(state)) ** 2 for state in joint_result.states]
plt.figure()
plt.plot(times, joint_p_R, label=r"$P_{R10}$")
plt.plot(times, joint_p_Ar0, label=r"$P_{Ar0}$")
plt.xlabel(r"Time ($\mu$s)")
plt.ylabel("Population")
plt.legend()
plt.tight_layout()

data_pi_time = np.pi / Omega
ancilla_2pi_time = 2 * np.pi / Omega_R
print(f"Paper-parameter gate duration: {2 * data_pi_time + ancilla_2pi_time:.3f} microseconds")
first_times = np.linspace(0, data_pi_time, 201)
middle_times = np.linspace(data_pi_time, data_pi_time + ancilla_2pi_time, 401)
last_times = np.linspace(data_pi_time + ancilla_2pi_time,
                         2 * data_pi_time + ancilla_2pi_time, 201)
first = sesolve(H_data_joint, A10, first_times)
middle = sesolve(H_ancilla_joint, first.states[-1], middle_times)
last = sesolve(H_data_joint, middle.states[-1], last_times)
gate_times = np.r_[first_times, middle_times[1:], last_times[1:]]
gate_states = first.states + middle.states[1:] + last.states[1:]
P_R_joint = tensor(R * R.dag(), qeye(3), qeye(3))
data_rydberg_count = tensor(qeye(2), n_r, qeye(3)) + tensor(qeye(2), qeye(3), n_r)
gate_p_R = expect(P_R_joint, gate_states)
gate_n_r = expect(data_rydberg_count, gate_states)
print(f"No interaction: ancilla P_R midway = {gate_p_R[400]:.6f}")
print(f"No interaction: final |A10> amplitude = {A10.overlap(last.states[-1]):.6f}")

plt.figure()
plt.plot(gate_times, gate_p_R, label=r"$P_R$")
plt.plot(gate_times, gate_n_r, label="Data Rydberg count")
plt.axvline(data_pi_time, color="gray", linestyle=":")
plt.axvline(data_pi_time + ancilla_2pi_time, color="gray", linestyle=":")
plt.xlabel(r"Time ($\mu$s)")
plt.ylabel("Population or excitation count")
plt.legend()
plt.tight_layout()

V = 2 * np.pi * 200
H_blockade = V * (tensor(R * R.dag(), n_r, qeye(3))
                  + tensor(R * R.dag(), qeye(3), n_r))
first_blocked = sesolve(H_data_joint + H_blockade, A10, first_times)
middle_blocked = sesolve(H_ancilla_joint + H_blockade,
                         first_blocked.states[-1], middle_times)
last_blocked = sesolve(H_data_joint + H_blockade,
                       middle_blocked.states[-1], last_times)
blocked_states = (first_blocked.states + middle_blocked.states[1:]
                  + last_blocked.states[1:])
blocked_p_R = expect(P_R_joint, blocked_states)
print(f"Effective blockade: peak ancilla P_R = {max(blocked_p_R):.6f}")

plt.figure()
plt.plot(gate_times, gate_p_R, label="No interaction")
plt.plot(gate_times, blocked_p_R, label=rf"Effective blockade $V/\Omega={V/Omega:g}$")
plt.axvline(data_pi_time, color="gray", linestyle=":")
plt.axvline(data_pi_time + ancilla_2pi_time, color="gray", linestyle=":")
plt.xlabel(r"Time ($\mu$s)")
plt.ylabel(r"Ancilla $P_R$")
plt.legend()
plt.tight_layout()

input_curves = {}
full_state_curves = {}
for label, data_state in inputs:
    initial = tensor(A, data_state)
    step1 = sesolve(H_data_joint + H_blockade, initial, first_times)
    step2 = sesolve(H_ancilla_joint + H_blockade, step1.states[-1], middle_times)
    step3 = sesolve(H_data_joint + H_blockade, step2.states[-1], last_times)
    full_states = step1.states + step2.states[1:] + step3.states[1:]
    full_state_curves[label] = full_states
    p_R_curve = expect(P_R_joint, full_states)
    input_curves[label] = p_R_curve
    print(f"|A{label}>: peak P_R={max(p_R_curve):.6f}, "
          f"final return={abs(initial.overlap(step3.states[-1]))**2:.6f}")

plt.figure()
for label, curve in input_curves.items():
    plt.plot(gate_times, curve, label=rf"$|A{label}\rangle$")
plt.axvline(data_pi_time, color="gray", linestyle=":")
plt.axvline(data_pi_time + ancilla_2pi_time, color="gray", linestyle=":")
plt.xlabel(r"Time ($\mu$s)")
plt.ylabel(r"Ancilla $P_R$")
plt.legend()
plt.tight_layout()

fig, axes = plt.subplots(3, 1, figsize=(8, 8), sharex=True)
for ax, label in zip(axes, ["00", "01", "11"]):
    data = [zero if bit == "0" else one for bit in label]
    excited = [zero if bit == "0" else r for bit in label]
    starting_state = tensor(A, data[0], data[1])
    raised_state = tensor(A, excited[0], excited[1])
    blocked_state = tensor(R, excited[0], excited[1])
    states = full_state_curves[label]
    ax.plot(gate_times, [abs(starting_state.overlap(s)) ** 2 for s in states],
            label=rf"$|A{label}\rangle$")
    if label != "00":
        ax.plot(gate_times, [abs(raised_state.overlap(s)) ** 2 for s in states],
                label="Data excited, ancilla A")
    ax.plot(gate_times, [abs(blocked_state.overlap(s)) ** 2 for s in states],
            label="Ancilla R")
    ax.axvline(data_pi_time, color="gray", linestyle=":")
    ax.axvline(data_pi_time + ancilla_2pi_time, color="gray", linestyle=":")
    ax.set_ylabel("Population")
    ax.legend(loc="center right")
axes[-1].set_xlabel(r"Time ($\mu$s)")
fig.suptitle("Figure 3(c) comparison: effective blockade model")
fig.tight_layout()
Path("figures").mkdir(exist_ok=True)
fig.savefig("figures/fig3c_effective.png", dpi=150)
plt.show()
