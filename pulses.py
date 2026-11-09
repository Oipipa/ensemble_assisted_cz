import numpy as np
from qutip import qeye, sesolve, tensor


def smooth_envelope(t, start, duration, area):
    x = (t - start) / duration
    if x < 0 or x > 1:
        return 0.0
    return (np.exp(-0.5 * ((x - 0.5) / 0.25) ** 2) - np.exp(-2)) / area


def variable_drive_outputs(omega, smooth, base_omega, data_drive,
                           ancilla_drive, interaction, computational,
                           rydberg, ancilla_r, area, detuning=0,
                           amplitude_error=0, duration_error=0,
                           interaction_ratio=1):
    td = np.pi / omega * (1 + duration_error)
    ta = 2 * np.pi / omega * (1 + duration_error)
    t1 = np.linspace(0, td, 121)
    t2 = np.linspace(td, td + ta, 241)
    t3 = np.linspace(td + ta, 2 * td + ta, 121)
    h_data = omega / base_omega * (1 + amplitude_error) * data_drive
    h_ancilla = omega / base_omega * (1 + amplitude_error) * ancilla_drive
    identity = qeye(rydberg.shape[0])
    data_r = rydberg * rydberg.dag()
    detuning_term = -detuning * (
        tensor(qeye(3), data_r, identity) + tensor(qeye(3), identity, data_r)
        + tensor(ancilla_r * ancilla_r.dag(), identity, identity))
    idle = interaction_ratio * interaction + detuning_term
    if smooth:
        first_h = [idle, [h_data, lambda t: smooth_envelope(t, 0, td, area)]]
        last_h = [idle, [h_data, lambda t: smooth_envelope(t, td + ta, td, area)]]
    else:
        first_h = h_data + idle
        last_h = first_h
    outputs = []
    for initial in computational:
        first = sesolve(first_h, initial, t1)
        middle = sesolve(h_ancilla + idle, first.states[-1], t2)
        last = sesolve(last_h, middle.states[-1], t3)
        outputs.append(last.states[-1])
    return outputs
