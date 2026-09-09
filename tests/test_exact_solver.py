import numpy as np
import pytest
from qiskit.quantum_info import SparsePauliOp

from analysis import compress_with_layout, expand_bitstring, solve_exactly, to_diagonal
from analysis.exact_solver import MAX_EXACT_QUBITS
from exceptions import InvalidOperatorError
from tests.helpers import build_hamiltonian


def test_identity_has_a_constant_diagonal():
    diagonal = to_diagonal(SparsePauliOp.from_list([("II", 2.5)]))

    assert np.allclose(diagonal, 2.5)


def test_single_z_splits_the_diagonal_by_parity():
    diagonal = to_diagonal(SparsePauliOp.from_list([("Z", 1.0)]))

    assert np.allclose(diagonal, [1.0, -1.0])


def test_diagonal_length_matches_the_register():
    diagonal = to_diagonal(SparsePauliOp.from_list([("ZZZ", 1.0)]))

    assert len(diagonal) == 8


def test_off_diagonal_operators_are_rejected():
    with pytest.raises(InvalidOperatorError, match="off-diagonal"):
        to_diagonal(SparsePauliOp.from_list([("X", 1.0)]))


def test_oversized_registers_are_refused():
    label = "Z" * (MAX_EXACT_QUBITS + 1)

    with pytest.raises(InvalidOperatorError, match="Refusing"):
        to_diagonal(SparsePauliOp.from_list([(label, 1.0)]))


def test_ground_state_is_the_smallest_diagonal_entry():
    operator = SparsePauliOp.from_list([("ZI", 1.0), ("IZ", 2.0)])

    solution = solve_exactly(operator)

    assert np.isclose(solution.energy, to_diagonal(operator).min())


def test_ground_state_bitstring_reproduces_the_energy():
    operator = SparsePauliOp.from_list([("ZI", 1.0), ("IZ", 2.0)])

    solution = solve_exactly(operator)

    assert np.isclose(
        to_diagonal(operator)[int(solution.bitstring, 2)], solution.energy
    )


def test_degeneracy_counts_the_minimisers():
    operator = SparsePauliOp.from_list([("ZZ", 1.0)])

    assert solve_exactly(operator).degeneracy == 2


def test_ground_state_has_rank_zero():
    operator = SparsePauliOp.from_list([("ZI", 1.0), ("IZ", 2.0)])
    solution = solve_exactly(operator)

    assert solution.rank_of(solution.bitstring) == 0


def test_compression_keeps_the_acting_qubits():
    operator = SparsePauliOp.from_list([("IZI", 1.0)])

    compressed, kept = compress_with_layout(operator)

    assert kept == [1]
    assert compressed.num_qubits == 1


def test_compression_preserves_the_spectrum_values():
    _, hamiltonian = build_hamiltonian("HPPHH")
    compressed, _ = compress_with_layout(hamiltonian)

    assert np.isclose(solve_exactly(compressed).energy, to_diagonal(compressed).min())


def test_expanding_restores_the_kept_qubits():
    expanded = expand_bitstring("1", kept_qubits=[1], num_qubits=3)

    assert expanded == "010"


def test_expanding_round_trips_through_compression():
    operator = SparsePauliOp.from_list([("IZI", 1.0), ("ZII", 2.0)])
    compressed, kept = compress_with_layout(operator)
    solution = solve_exactly(compressed)

    expanded = expand_bitstring(solution.bitstring, kept, operator.num_qubits)

    assert np.isclose(to_diagonal(operator)[int(expanded, 2)], solution.energy)


def test_length_mismatch_is_rejected():
    with pytest.raises(ValueError, match="does not match"):
        expand_bitstring("101", kept_qubits=[0], num_qubits=3)
