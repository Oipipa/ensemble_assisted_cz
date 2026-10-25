import matplotlib.pyplot as plt
import numpy as np
from qutip import basis, qeye, sesolve, tensor


zero = basis(3, 0)
one = basis(3, 1)
r = basis(3, 2)
Omega = 1.0
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
plt.xlabel("Time")
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

plt.xlabel("Time")
plt.ylabel("Mean number of Rydberg excitations")
plt.legend()
plt.tight_layout()

A = basis(2, 0)
R = basis(2, 1)
Omega_R = 1.0
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
plt.xlabel("Time")
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
plt.xlabel("Time")
plt.ylabel("Population")
plt.legend()
plt.tight_layout()
plt.show()
