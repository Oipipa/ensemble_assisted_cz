from qutip import basis, qeye, tensor


def forster_operators(omega_data, omega_ancilla, interaction):
    ancilla_a, ancilla_r, ancilla_rp = [basis(3, i) for i in range(3)]
    zero, one, rydberg, rydberg_prime = [basis(4, i) for i in range(4)]
    identity = qeye(4)
    data_drive = omega_data / 2 * (one * rydberg.dag() + rydberg * one.dag())
    ancilla_drive = omega_ancilla / 2 * (ancilla_a * ancilla_r.dag()
                                       + ancilla_r * ancilla_a.dag())
    data_full = (tensor(qeye(3), data_drive, identity)
                 + tensor(qeye(3), identity, data_drive))
    ancilla_full = tensor(ancilla_drive, identity, identity)
    pair1 = tensor(ancilla_rp * ancilla_r.dag(),
                   rydberg_prime * rydberg.dag(), identity)
    pair2 = tensor(ancilla_rp * ancilla_r.dag(),
                   identity, rydberg_prime * rydberg.dag())
    forster_full = interaction * (pair1 + pair1.dag() + pair2 + pair2.dag())
    computational = [tensor(ancilla_a, zero, zero),
                     tensor(ancilla_a, zero, one),
                     tensor(ancilla_a, one, zero),
                     tensor(ancilla_a, one, one)]
    return (ancilla_a, ancilla_r, ancilla_rp, zero, one, rydberg,
            rydberg_prime, identity, data_full, ancilla_full,
            forster_full, computational)
