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
final_gate_states = {}
for label, data_state in inputs:
    initial = tensor(A, data_state)
    step1 = sesolve(H_data_joint + H_blockade, initial, first_times)
    step2 = sesolve(H_ancilla_joint + H_blockade, step1.states[-1], middle_times)
    step3 = sesolve(H_data_joint + H_blockade, step2.states[-1], last_times)
    full_states = step1.states + step2.states[1:] + step3.states[1:]
    full_state_curves[label] = full_states
    final_gate_states[label] = step3.states[-1]
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

final_amplitudes = {}
for label, data_state in inputs:
    initial = tensor(A, data_state)
    final_amplitudes[label] = initial.overlap(final_gate_states[label])
reference = final_amplitudes["00"]
relative_phases = [np.angle(final_amplitudes[label] / reference)
                   for label, _ in inputs]
returns = [abs(final_amplitudes[label]) ** 2 for label, _ in inputs]
for (label, _), phase, returned in zip(inputs, relative_phases, returns):
    print(f"|A{label}>: return={returned:.6f}, phase relative to 00={phase:.6f} rad")

phase_fig, phase_axes = plt.subplots(1, 2, figsize=(8, 3))
phase_axes[0].bar([label for label, _ in inputs], returns)
phase_axes[0].set_ylabel("Return probability")
phase_axes[0].set_ylim(0, 1.05)
phase_axes[1].bar([label for label, _ in inputs], relative_phases)
phase_axes[1].set_ylabel("Phase relative to 00 (rad)")
phase_fig.tight_layout()

computational_states = [tensor(A, data_state) for _, data_state in inputs]
effective_gate = np.array([
    [bra.overlap(final_gate_states[label]) for label, _ in inputs]
    for bra in computational_states
])
Z1Z2 = np.diag([1, -1, -1, 1])
paper_operation = Z1Z2 @ np.diag([-1, 1, 1, 1])
CZ = np.diag([1, 1, 1, -1])
leakage = 1 - np.sum(abs(effective_gate) ** 2, axis=0)
print("Effective computational operation:")
print(np.round(effective_gate, 4))
print(f"Maximum difference from paper Eq. (4): "
      f"{np.max(abs(effective_gate - paper_operation)):.4f}")
print(f"Maximum computational leakage: {max(leakage):.6f}")
assert np.array_equal(paper_operation, -CZ)

plus_plus = sum(computational_states[1:], computational_states[0]) / 2
super_first = sesolve(H_data_joint + H_blockade, plus_plus, first_times)
super_middle = sesolve(H_ancilla_joint + H_blockade,
                       super_first.states[-1], middle_times)
super_last = sesolve(H_data_joint + H_blockade,
                     super_middle.states[-1], last_times)
ideal_super_output = -plus_plus + computational_states[3]
super_fidelity = abs(ideal_super_output.overlap(super_last.states[-1])) ** 2
print(f"|++> output fidelity to ideal Eq. (4): {super_fidelity:.6f}")

for Hamiltonian in [H_data_joint, H_ancilla_joint, H_blockade]:
    assert Hamiltonian.isherm
assert np.isclose(Omega * data_pi_time, np.pi)
assert np.isclose(Omega_R * ancilla_2pi_time, 2 * np.pi)
max_norm_error = max(abs(state.norm() - 1)
                     for states in full_state_curves.values() for state in states)
print(f"Maximum state-norm error: {max_norm_error:.2e}")
assert max_norm_error < 1e-8
assert min(returns) > 0.999
assert super_fidelity > 0.998

fine_first_times = np.linspace(0, data_pi_time, 401)
fine_middle_times = np.linspace(data_pi_time, data_pi_time + ancilla_2pi_time, 801)
fine_last_times = np.linspace(data_pi_time + ancilla_2pi_time,
                              2 * data_pi_time + ancilla_2pi_time, 401)
fine_initial = computational_states[1]
fine_first = sesolve(H_data_joint + H_blockade, fine_initial, fine_first_times)
fine_middle = sesolve(H_ancilla_joint + H_blockade,
                      fine_first.states[-1], fine_middle_times)
fine_last = sesolve(H_data_joint + H_blockade,
                    fine_middle.states[-1], fine_last_times)
convergence_error = (fine_last.states[-1] - final_gate_states["01"]).norm()
print(f"Fine-grid final-state difference: {convergence_error:.2e}")
assert convergence_error < 1e-4


def run_effective_gate(initial, V_ratio):
    blockade = V_ratio * Omega * (H_blockade / V)
    start = sesolve(H_data_joint + blockade, initial, first_times)
    middle = sesolve(H_ancilla_joint + blockade, start.states[-1], middle_times)
    end = sesolve(H_data_joint + blockade, middle.states[-1], last_times)
    return end.states[-1]


interaction_ratios = [2, 5, 10, 20, 40, 80]
leakage_by_V = []
phase_error_by_V = []
gate_error_by_V = []
for V_ratio in interaction_ratios:
    outputs = [run_effective_gate(state, V_ratio) for state in computational_states]
    matrix = np.array([[bra.overlap(ket) for ket in outputs]
                       for bra in computational_states])
    leakage_by_V.append(max(1 - np.sum(abs(matrix) ** 2, axis=0)))
    phase_error_by_V.append(max(abs(np.angle(np.diag(matrix) / np.diag(paper_operation)))))
    gate_error_by_V.append(np.max(abs(matrix - paper_operation)))
    print(f"V/Omega={V_ratio:>2}: leakage={leakage_by_V[-1]:.2e}, "
          f"phase error={phase_error_by_V[-1]:.4f} rad, "
          f"matrix error={gate_error_by_V[-1]:.4f}")

error_fig, error_axes = plt.subplots(1, 2, figsize=(9, 3))
error_axes[0].plot(interaction_ratios, leakage_by_V, "o-")
error_axes[0].set_ylabel("Max computational leakage")
error_axes[1].plot(interaction_ratios, phase_error_by_V, "o-", label="Phase (rad)")
error_axes[1].plot(interaction_ratios, gate_error_by_V, "s--", label="Matrix deviation")
error_axes[1].set_ylabel("Error")
error_axes[1].legend()
for ax in error_axes:
    ax.set_xlabel(r"$V/\Omega$")
error_fig.tight_layout()

from qutip import mesolve

zero4, one4, r4, lost4 = [basis(4, i) for i in range(4)]
I4 = qeye(4)
n4 = r4 * r4.dag()
loss4 = lost4 * lost4.dag()
H_data4 = Omega / 2 * (one4 * r4.dag() + r4 * one4.dag())
H_joint_data4 = tensor(qeye(2), H_data4, I4) + tensor(qeye(2), I4, H_data4)
H_joint_ancilla4 = tensor(H_ancilla, I4, I4)
H_joint_blockade4 = V * (tensor(R * R.dag(), n4, I4)
                          + tensor(R * R.dag(), I4, n4))
Gamma_r = 1 / 100  # per microsecond; all r decay goes to an explicit loss state
jump4 = lost4 * r4.dag()
data_jumps = [np.sqrt(Gamma_r) * tensor(qeye(2), jump4, I4),
              np.sqrt(Gamma_r) * tensor(qeye(2), I4, jump4)]
initial4 = tensor(A, one4, one4)
rho_initial4 = initial4 * initial4.dag()
noisy_first = mesolve(H_joint_data4 + H_joint_blockade4, rho_initial4,
                      first_times, c_ops=data_jumps)
noisy_middle = mesolve(H_joint_ancilla4 + H_joint_blockade4,
                       noisy_first.states[-1], middle_times, c_ops=data_jumps)
noisy_last = mesolve(H_joint_data4 + H_joint_blockade4,
                     noisy_middle.states[-1], last_times, c_ops=data_jumps)
noisy_states = noisy_first.states + noisy_middle.states[1:] + noisy_last.states[1:]
loss_count4 = tensor(qeye(2), loss4, I4) + tensor(qeye(2), I4, loss4)
noisy_loss = expect(loss_count4, noisy_states)
final_survival = expect(initial4 * initial4.dag(), noisy_last.states[-1])
print(f"Data decay: final |A11> return={final_survival:.6f}, "
      f"mean lost atoms={noisy_loss[-1]:.6f}")
assert abs(noisy_last.states[-1].tr() - 1) < 1e-8
assert final_survival < returns[3]

plt.figure()
plt.plot(gate_times, noisy_loss, label="Mean lost data atoms")
plt.plot(gate_times, expect(initial4 * initial4.dag(), noisy_states),
         label=r"$P_{A11}$ with decay")
plt.xlabel(r"Time ($\mu$s)")
plt.ylabel("Population or loss count")
plt.legend()
plt.tight_layout()

A3, R3, lostA3 = [basis(3, i) for i in range(3)]
I3 = qeye(3)
H_ancilla3 = Omega_R / 2 * (A3 * R3.dag() + R3 * A3.dag())
H_data43 = tensor(I3, H_data4, I4) + tensor(I3, I4, H_data4)
H_ancilla43 = tensor(H_ancilla3, I4, I4)
H_block43 = V * (tensor(R3 * R3.dag(), n4, I4)
                 + tensor(R3 * R3.dag(), I4, n4))
data_jumps43 = [np.sqrt(Gamma_r) * tensor(I3, jump4, I4),
                np.sqrt(Gamma_r) * tensor(I3, I4, jump4)]
Gamma_ancilla = 1 / 200
ancilla_jumps43 = [np.sqrt(Gamma_ancilla) * tensor(lostA3 * A3.dag(), I4, I4),
                   np.sqrt(Gamma_ancilla) * tensor(lostA3 * R3.dag(), I4, I4)]


def evolve_open_gate(initial, jumps):
    rho = initial * initial.dag()
    first = mesolve(H_data43 + H_block43, rho, first_times, c_ops=jumps)
    middle = mesolve(H_ancilla43 + H_block43,
                     first.states[-1], middle_times, c_ops=jumps)
    last = mesolve(H_data43 + H_block43,
                   middle.states[-1], last_times, c_ops=jumps)
    return first.states + middle.states[1:] + last.states[1:]


initial_A11 = tensor(A3, one4, one4)
initial_A00 = tensor(A3, zero4, zero4)
ancilla_only = evolve_open_gate(initial_A11, ancilla_jumps43)
combined_decay = evolve_open_gate(initial_A11, data_jumps43 + ancilla_jumps43)
ancilla_00 = evolve_open_gate(initial_A00, ancilla_jumps43)
P_A43 = tensor(A3 * A3.dag(), I4, I4)
P_R43 = tensor(R3 * R3.dag(), I4, I4)
P_lostA43 = tensor(lostA3 * lostA3.dag(), I4, I4)
P_lost_data43 = tensor(I3, loss4, I4) + tensor(I3, I4, loss4)
for label, states in [("A11", ancilla_only), ("A00", ancilla_00)]:
    loss_from_A = Gamma_ancilla * np.trapezoid(expect(P_A43, states), gate_times)
    loss_from_R = Gamma_ancilla * np.trapezoid(expect(P_R43, states), gate_times)
    print(f"Ancilla-only |{label}>: loss from A={loss_from_A:.6f}, "
          f"loss from R={loss_from_R:.6f}")

return_ancilla_only = expect(initial_A11 * initial_A11.dag(), ancilla_only[-1])
return_combined = expect(initial_A11 * initial_A11.dag(), combined_decay[-1])
print(f"|A11> return: ancilla decay={return_ancilla_only:.6f}, "
      f"data + ancilla decay={return_combined:.6f}")
assert return_combined < return_ancilla_only < returns[3]

plt.figure()
plt.plot(gate_times, expect(P_lostA43, ancilla_only), label="Ancilla loss")
plt.plot(gate_times, expect(P_lost_data43, combined_decay), label="Data loss (both)")
plt.xlabel(r"Time ($\mu$s)")
plt.ylabel("Loss probability or count")
plt.legend()
plt.tight_layout()


def evolve_open_operator(rho, jumps):
    short_first = [0, data_pi_time]
    short_middle = [data_pi_time, data_pi_time + ancilla_2pi_time]
    short_last = [data_pi_time + ancilla_2pi_time,
                  2 * data_pi_time + ancilla_2pi_time]
    options = {"normalize_output": False, "rtol": 1e-9, "atol": 1e-11}
    first = mesolve(H_data43 + H_block43, rho, short_first,
                    c_ops=jumps, options=options)
    middle = mesolve(H_ancilla43 + H_block43, first.states[-1], short_middle,
                     c_ops=jumps, options=options)
    last = mesolve(H_data43 + H_block43, middle.states[-1], short_last,
                   c_ops=jumps, options=options)
    return last.states[-1]


def data_channel(rho_data):
    rho_full = tensor(A3 * A3.dag(), rho_data)
    final_full = evolve_open_operator(rho_full, data_jumps43 + ancilla_jumps43)
    return final_full.ptrace([1, 2])


data_comp_states = [tensor(zero4, zero4), tensor(zero4, one4),
                    tensor(one4, zero4), tensor(one4, one4)]
channel_basis = {}
for i, ket in enumerate(data_comp_states):
    for j, bra in enumerate(data_comp_states):
        channel_basis[i, j] = data_channel(ket * bra.dag())

data_plus_plus = sum(data_comp_states[1:], data_comp_states[0]) / 2
rho_plus_plus = data_plus_plus * data_plus_plus.dag()
direct_output = data_channel(rho_plus_plus)
linear_output = sum(channel_basis.values()) / 4
linearity_error = (direct_output - linear_output).norm()
P_comp_data = tensor(zero4 * zero4.dag() + one4 * one4.dag(),
                     zero4 * zero4.dag() + one4 * one4.dag())
channel_leakage = 1 - expect(P_comp_data, direct_output)
print(f"Noisy channel: trace={direct_output.tr():.6f}, "
      f"linearity error={linearity_error:.2e}, leakage={channel_leakage:.6f}")
assert abs(direct_output.tr() - 1) < 1e-6
assert direct_output.isherm
assert linearity_error < 1e-6
assert channel_leakage > 0
plt.show()
