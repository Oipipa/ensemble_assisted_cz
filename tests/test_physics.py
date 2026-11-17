import numpy as np
from qutip import basis, qeye, tensor

from fidelity import average_from_basis, average_from_pure_outputs
from models import forster_operators
from pulses import smooth_envelope, variable_drive_outputs
from reproduce import gate_point


OMEGA = 2 * np.pi * 10
INTERACTION = 2 * np.pi * 200


def test_forster_dark_state():
    coupling = np.array([[0, OMEGA / 2, 0],
                         [OMEGA / 2, 0, INTERACTION],
                         [0, INTERACTION, 0]])
    dark = np.array([2 * INTERACTION, 0, -OMEGA])
    dark /= np.linalg.norm(dark)
    assert np.linalg.norm(coupling @ dark) < 1e-11
    assert np.allclose(np.linalg.eigvalsh(coupling),
                       [-np.hypot(OMEGA, 2 * INTERACTION) / 2, 0,
                        np.hypot(OMEGA, 2 * INTERACTION) / 2])


def test_pulse_area_and_gate_fidelity():
    (ancilla_a, ancilla_r, _, zero, one, rydberg, _, identity,
     data_drive, ancilla_drive, interaction, states) = forster_operators(
        OMEGA, OMEGA, INTERACTION)
    assert data_drive.isherm and ancilla_drive.isherm and interaction.isherm
    x = np.linspace(0, 1, 2001)
    area = np.trapezoid(np.exp(-0.5 * ((x - 0.5) / 0.25) ** 2) - np.exp(-2), x)
    duration = np.pi / OMEGA
    times = np.linspace(0, duration, 2001)
    assert abs(OMEGA * np.trapezoid(
        [smooth_envelope(t, 0, duration, area) for t in times], times)
        - np.pi) < 1e-5

    data_states = [tensor(zero, zero), tensor(zero, one),
                   tensor(one, zero), tensor(one, one)]
    targets = [-data_states[0], -data_states[1],
               -data_states[2], data_states[3]]
    projector = sum(state * state.dag() for state in data_states)
    ideal_basis = {(i, j): targets[i] * targets[j].dag()
                   for i in range(4) for j in range(4)}
    assert abs(average_from_basis(ideal_basis, targets, projector) - 1) < 1e-12

    outputs = variable_drive_outputs(OMEGA, False, OMEGA, data_drive,
                                     ancilla_drive, interaction, states,
                                     rydberg, ancilla_r, area)
    fidelity = average_from_pure_outputs(outputs, targets, projector)
    assert 0.9997 < fidelity < 1


def test_detuning_and_timing_lower_fidelity():
    (_, ancilla_r, _, zero, one, rydberg, _, _, data_drive,
     ancilla_drive, interaction, states) = forster_operators(
        OMEGA, OMEGA, INTERACTION)
    x = np.linspace(0, 1, 2001)
    area = np.trapezoid(np.exp(-0.5 * ((x - 0.5) / 0.25) ** 2) - np.exp(-2), x)
    data_states = [tensor(zero, zero), tensor(zero, one),
                   tensor(one, zero), tensor(one, one)]
    targets = [-data_states[0], -data_states[1],
               -data_states[2], data_states[3]]
    projector = sum(state * state.dag() for state in data_states)
    baseline = variable_drive_outputs(OMEGA, False, OMEGA, data_drive,
                                      ancilla_drive, interaction, states,
                                      rydberg, ancilla_r, area)
    mistuned = variable_drive_outputs(OMEGA, False, OMEGA, data_drive,
                                      ancilla_drive, interaction, states,
                                      rydberg, ancilla_r, area,
                                      detuning=0.1 * OMEGA,
                                      duration_error=0.1)
    assert average_from_pure_outputs(mistuned, targets, projector) < 0.9
    assert average_from_pure_outputs(baseline, targets, projector) > 0.9997


def test_noisy_forster_reproduces_table_and_ancilla_shaping_suppresses_dip():
    square_10, _ = gate_point(10)
    assert abs(square_10 - 0.9980) < 0.0001
    noiseless_10, _ = gate_point(10, noisy=False)
    data_only, _ = gate_point(10, gamma_ancilla=0)
    ancilla_only, _ = gate_point(10, gamma_data=0)
    assert abs((noiseless_10 - data_only) - 0.0015) < 0.00002
    assert abs((noiseless_10 - ancilla_only) - 0.0003) < 0.00002
    square_30, square_leakage = gate_point(30, noisy=False)
    smooth_30, smooth_leakage = gate_point(30, shape="ancilla", noisy=False)
    assert square_30 < 0.99
    assert smooth_30 > 0.9999
    assert square_leakage > 0.01
    assert smooth_leakage < 0.0001
