import numpy as np
import pytest

from constants import QUBITS_PER_TURN
from particle import Ligand
from tests.helpers import diagonal_of
from utils.lattice_utils import build_squared_distance, build_walk_position


def test_walk_qubits_scale_with_steps():
    assert Ligand(num_steps=1).num_walk_qubits == QUBITS_PER_TURN
    assert Ligand(num_steps=3).num_walk_qubits == 3 * QUBITS_PER_TURN


@pytest.mark.parametrize(("steps", "parity"), [(1, 1), (2, 0), (3, 1), (4, 0)])
def test_sublattice_parity_matches_the_step_count(steps, parity):
    assert Ligand(num_steps=steps).sublattice_parity == parity


def test_even_ligand_can_only_reach_odd_residues():
    assert Ligand(num_steps=2).eligible_bead_indices(7) == [1, 3, 5]


def test_odd_ligand_can_only_reach_even_residues():
    assert Ligand(num_steps=3).eligible_bead_indices(7) == [0, 2, 4, 6]


def test_the_two_parities_partition_the_chain():
    even = set(Ligand(num_steps=2).eligible_bead_indices(7))
    odd = set(Ligand(num_steps=3).eligible_bead_indices(7))

    assert even.isdisjoint(odd)
    assert even | odd == set(range(7))


def test_position_operators_span_the_lattice_axes():
    ligand = Ligand(num_steps=2)

    position = ligand.position_operators(
        num_qubits=ligand.num_walk_qubits, walk_qubit_offset=0
    )

    assert len(position) == 4


def test_single_step_ligand_sits_one_site_from_the_origin():
    ligand = Ligand(num_steps=1)
    origin = build_walk_position(qubit_bases=[], num_qubits=ligand.num_walk_qubits)

    position = ligand.position_operators(
        num_qubits=ligand.num_walk_qubits, walk_qubit_offset=0
    )
    squared = diagonal_of(build_squared_distance(position, origin))

    assert np.allclose(squared, 1.0)


def test_two_step_ligand_stays_on_the_origin_sublattice():
    ligand = Ligand(num_steps=2)
    origin = build_walk_position(qubit_bases=[], num_qubits=ligand.num_walk_qubits)

    position = ligand.position_operators(
        num_qubits=ligand.num_walk_qubits, walk_qubit_offset=0
    )
    squared = diagonal_of(build_squared_distance(position, origin))

    assert np.all(np.isclose(np.round(squared) % 2, 0))


def test_walk_offset_moves_the_operators_onto_other_qubits():
    ligand = Ligand(num_steps=1)

    shifted = ligand.position_operators(num_qubits=8, walk_qubit_offset=4)

    acted_on = np.any([pauli.z for pauli in shifted[1].paulis], axis=0)
    assert not acted_on[:4].any()


@pytest.mark.parametrize("num_steps", [0, -1])
def test_non_positive_step_count_is_rejected(num_steps):
    with pytest.raises(ValueError, match="positive integer"):
        Ligand(num_steps=num_steps)


def test_repr_mentions_the_symbol_and_steps():
    assert "num_steps=2" in repr(Ligand(symbol="X", num_steps=2))
