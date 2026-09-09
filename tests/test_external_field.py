import numpy as np
import pytest

from analysis import (
    compress_with_layout,
    decode_structure,
    expand_bitstring,
    solve_exactly,
)
from enums import FieldMode
from particle import ExternalField
from tests.helpers import build_hamiltonian, diagonal_of

SEQUENCE = "APRLRFY"


def ground_state(**kwargs):
    _, hamiltonian = build_hamiltonian(SEQUENCE, **kwargs)
    compressed, _ = compress_with_layout(hamiltonian)
    return solve_exactly(compressed)


def test_uniform_field_is_a_multiple_of_the_identity():
    field = ExternalField.uniform(strength=2.0)

    operator = field.build_hamiltonian(chain_length=5, num_qubits=4)
    diagonal = diagonal_of(operator)

    assert np.allclose(diagonal, diagonal[0])
    assert np.isclose(diagonal[0], 10.0)


def test_gradient_field_is_not_a_multiple_of_the_identity():
    field = ExternalField.gradient(strength=1.0)

    diagonal = diagonal_of(field.build_hamiltonian(chain_length=5, num_qubits=8))

    assert not np.allclose(diagonal, diagonal[0])


def test_uniform_field_cannot_change_the_fold():
    baseline = ground_state()

    folds = {
        ground_state(external_field=ExternalField.uniform(strength)).bitstring
        for strength in (0.5, 2.0, 10.0, -10.0)
    }

    assert folds == {baseline.bitstring}


def test_uniform_field_shifts_every_energy_equally():
    baseline = ground_state()
    shifted = ground_state(external_field=ExternalField.uniform(1.0))

    assert np.allclose(shifted.spectrum - baseline.spectrum, len(SEQUENCE))


def test_gradient_field_changes_the_fold_at_sufficient_strength():
    baseline = ground_state()

    folds = {
        ground_state(external_field=ExternalField.gradient(strength)).bitstring
        for strength in (-2.0, -0.5, -0.1, 0.0, 0.5)
    }

    assert baseline.bitstring in folds
    assert len(folds) > 1


def test_gradient_field_stretches_the_chain_along_its_direction():
    spans = []
    for strength in (0.0, -2.0):
        _, hamiltonian = build_hamiltonian(
            SEQUENCE, external_field=ExternalField.gradient(strength)
        )
        compressed, kept = compress_with_layout(hamiltonian)
        solution = solve_exactly(compressed)
        structure = decode_structure(
            bitstring=expand_bitstring(
                solution.bitstring, kept, hamiltonian.num_qubits
            ),
            chain_length=len(SEQUENCE),
        )
        spans.append(np.ptp(structure.coordinates @ np.array([0.0, 0.0, 1.0])))

    assert spans[1] > spans[0]


def test_axis_couplings_project_onto_the_lattice_basis():
    field = ExternalField.gradient(1.0, direction=np.array([0.0, 0.0, 1.0]))

    couplings = field.axis_couplings()

    assert couplings.shape == (4,)
    assert np.allclose(np.abs(couplings), 1.0 / np.sqrt(3))


def test_direction_is_normalised():
    field = ExternalField.gradient(1.0, direction=np.array([0.0, 0.0, 5.0]))

    assert np.isclose(np.linalg.norm(field.direction), 1.0)


def test_factories_set_the_mode():
    assert ExternalField.uniform(1.0).mode is FieldMode.UNIFORM
    assert ExternalField.gradient(1.0).mode is FieldMode.GRADIENT


@pytest.mark.parametrize("strength", [float("nan"), float("inf")])
def test_non_finite_strength_is_rejected(strength):
    with pytest.raises(ValueError, match="finite"):
        ExternalField.uniform(strength)


def test_zero_direction_is_rejected():
    with pytest.raises(ValueError, match="non-zero"):
        ExternalField.gradient(1.0, direction=np.zeros(3))


def test_wrong_shaped_direction_is_rejected():
    with pytest.raises(ValueError, match="shape"):
        ExternalField.gradient(1.0, direction=np.array([1.0, 0.0]))
