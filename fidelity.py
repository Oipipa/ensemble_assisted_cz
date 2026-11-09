import numpy as np
from qutip import expect


def average_from_basis(basis_outputs, target_states, projector):
    identity = sum(basis_outputs[i, i] for i in range(4))
    term_one = expect(projector, identity)
    term_two = sum(
        target_states[i].overlap(basis_outputs[i, j] * target_states[j])
        for i in range(4) for j in range(4))
    return np.real((term_one + term_two) / 20)


def average_from_pure_outputs(outputs, target_states, projector):
    basis_map = {(i, j): (outputs[i] * outputs[j].dag()).ptrace([1, 2])
                 for i in range(4) for j in range(4)}
    return average_from_basis(basis_map, target_states, projector)
