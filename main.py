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
plt.show()
