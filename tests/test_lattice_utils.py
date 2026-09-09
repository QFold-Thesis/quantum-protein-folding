import numpy as np
import pytest

from constants import CONFORMATION_ENCODING, QUBITS_PER_TURN
from distance import DistanceMap
from exceptions import ConformationEncodingError
from protein import Protein
from tests.helpers import diagonal_of
from utils.lattice_utils import (
    build_chain_position,
    build_squared_distance,
    build_step_indicators,
    build_walk_position,
)
from utils.qubit_utils import fix_qubits


def test_step_indicators_are_one_hot():
    indicators = build_step_indicators(
        qubit_base=0, num_qubits=QUBITS_PER_TURN, encoding=CONFORMATION_ENCODING
    )

    diagonals = np.array([diagonal_of(indicator) for indicator in indicators])

    assert diagonals.shape[0] == 4
    assert np.allclose(diagonals.sum(axis=0), 1.0)
    assert np.all((np.isclose(diagonals, 0.0)) | (np.isclose(diagonals, 1.0)))


def test_step_indicators_select_distinct_directions():
    indicators = build_step_indicators(
        qubit_base=0, num_qubits=QUBITS_PER_TURN, encoding=CONFORMATION_ENCODING
    )

    active = [int(np.argmax(diagonal_of(indicator))) for indicator in indicators]

    assert len(set(active)) == 4


def test_walk_of_one_step_is_always_unit_distance_from_origin():
    origin = build_walk_position(qubit_bases=[], num_qubits=QUBITS_PER_TURN)
    stepped = build_walk_position(qubit_bases=[0], num_qubits=QUBITS_PER_TURN)

    squared = diagonal_of(build_squared_distance(stepped, origin))

    assert np.allclose(squared, 1.0)


def test_empty_walk_stays_at_the_origin():
    position = build_walk_position(qubit_bases=[], num_qubits=4)

    for component in position:
        assert np.allclose(diagonal_of(component), 0.0)


def test_sign_offset_flips_the_sublattice_step():
    forward = build_walk_position(qubit_bases=[0], num_qubits=QUBITS_PER_TURN)
    flipped = build_walk_position(
        qubit_bases=[0], num_qubits=QUBITS_PER_TURN, sign_offset=1
    )

    for ahead, behind in zip(forward, flipped, strict=True):
        assert np.allclose(diagonal_of(ahead), -diagonal_of(behind))


@pytest.mark.parametrize("sequence", ["HPPHH", "APRLRF"])
def test_chain_positions_reproduce_the_distance_map(sequence):
    protein = Protein(
        main_protein_sequence=sequence,
        side_protein_sequence="_" * len(sequence),
        valid_symbols=set(sequence),
    )
    distance_map = DistanceMap(protein=protein)
    num_qubits = (len(sequence) - 1) * QUBITS_PER_TURN

    positions = [
        build_chain_position(index, num_qubits) for index in range(len(sequence))
    ]

    for lower in range(len(sequence)):
        for upper in range(lower + 1, len(sequence)):
            mine = fix_qubits(
                build_squared_distance(positions[lower], positions[upper])
            )
            assert np.allclose(
                diagonal_of(mine.simplify()),
                diagonal_of(distance_map[lower][upper].simplify()),
            )


def test_consecutive_beads_are_lattice_neighbours():
    num_qubits = 4 * QUBITS_PER_TURN
    positions = [build_chain_position(index, num_qubits) for index in range(5)]

    for index in range(4):
        squared = diagonal_of(
            build_squared_distance(positions[index], positions[index + 1])
        )
        assert np.allclose(squared, 1.0)


def test_squared_distance_is_never_negative():
    num_qubits = 3 * QUBITS_PER_TURN
    positions = [build_chain_position(index, num_qubits) for index in range(4)]

    squared = diagonal_of(build_squared_distance(positions[0], positions[3]))

    assert np.all(squared >= -1e-9)


def test_same_sublattice_distances_are_even():
    num_qubits = 4 * QUBITS_PER_TURN
    positions = [build_chain_position(index, num_qubits) for index in range(5)]

    squared = diagonal_of(build_squared_distance(positions[0], positions[2]))

    assert np.allclose(squared, np.round(squared))
    assert np.all(np.isclose(np.round(squared) % 2, 0))


def test_unknown_encoding_is_rejected():
    with pytest.raises(ConformationEncodingError):
        build_step_indicators(qubit_base=0, num_qubits=2, encoding="dense")
