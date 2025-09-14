from __future__ import annotations

import contextlib

import numpy as np
from qiskit.quantum_info import (  # pyright: ignore[reportMissingTypeStubs]
    Pauli,
    SparsePauliOp,
)

from constants import (
    NORM_FACTOR,
    SIGN_FLIP_SECOND_QUBIT_INDEX,
    SIGN_FLIP_SIXTH_QUBIT_INDEX,
)


def build_full_identity(num_qubits: int) -> SparsePauliOp:
    """Builds a full identity Pauli operator for a given number of qubits."""
    identity_string: str = "I" * num_qubits
    return SparsePauliOp.from_list([(identity_string, 1.0)])


def build_turn_qubit(z_index: int, num_qubits: int) -> SparsePauliOp:
    """Builds a turn qubit Pauli operator with Z at the specified index."""
    z_operator: SparsePauliOp = SparsePauliOp.from_sparse_list(
        [("Z", [z_index], 1.0)], num_qubits=num_qubits
    )

    full_identity: SparsePauliOp = build_full_identity(num_qubits=num_qubits)

    return NORM_FACTOR * (full_identity - z_operator)


def build_pauli_z_operator(num_qubits: int, pauli_z_indices: set[int]) -> SparsePauliOp:
    if not pauli_z_indices:
        return SparsePauliOp.from_list([("I" * num_qubits, 1.0)])

    idx_sorted = sorted(pauli_z_indices)
    local_label = "Z" * len(idx_sorted)

    return SparsePauliOp.from_sparse_list(
        [(local_label, idx_sorted, 1.0)],
        num_qubits=num_qubits,
    )


def convert_to_qubits(pauli_op: SparsePauliOp) -> SparsePauliOp:
    num_qubits: int = int(pauli_op.num_qubits)  # pyright: ignore[reportArgumentType]
    full_id: SparsePauliOp = SparsePauliOp.from_list([("I" * num_qubits, 1.0)])

    return NORM_FACTOR * (full_id - pauli_op)


def fix_qubits(
    operator: int | SparsePauliOp | Pauli,
    *,
    has_side_chain_second_bead: bool = False,
) -> int | SparsePauliOp | Pauli:
    """
    Assigns predefined values for turns qubits on positions 0, 1, 2, 3, 5 in the main chain.
    Qubits on these position are considered fixed and not subject to optimization.
    """
    # return if operator is int (might be 0 because it is initialized as operator = 0)
    if not isinstance(operator, (SparsePauliOp, Pauli)):
        return operator

    # Normalize operators to SparsePauliOp
    coeffs: np.ndarray = np.array([1.0])
    if isinstance(operator, Pauli):
        operator = SparsePauliOp([operator], coeffs)

    # Check if only one Pauli term is present, if so don't negate coeffs
    if len(operator.paulis) == 1:
        table_z = np.copy(operator.paulis[0].z)
        table_x = np.copy(operator.paulis[0].x)
        _preset_binary_vals(
            table_z, has_side_chain_second_bead=has_side_chain_second_bead
        )
        return SparsePauliOp(
            [Pauli((table_z, table_x))], [operator.coeffs[0]]
        ).simplify()

    new_paulis: list[Pauli] = []
    new_coeffs: np.ndarray = np.array([])

    for idx, pauli in enumerate(operator.paulis):  # pyright: ignore[reportArgumentType]
        table_z = np.copy(pauli.z)
        table_x = np.copy(pauli.x)
        coeff = operator.coeffs[idx]

        coeff = _calc_updated_coeffs(
            table_z, coeff, has_side_chain_second_bead=has_side_chain_second_bead
        )
        _preset_binary_vals(
            table_z, has_side_chain_second_bead=has_side_chain_second_bead
        )
        # Create new Pauli from updated z, x
        new_pauli = Pauli((table_z, table_x))
        new_paulis.append(new_pauli)
        new_coeffs = np.append(new_coeffs, coeff)

    return SparsePauliOp(new_paulis, new_coeffs).simplify()


def _calc_updated_coeffs(
    table_z: np.ndarray, coeff: float, *, has_side_chain_second_bead: bool
) -> float:
    """
    Update coefficients based on fixed qubit positions. Negate if appropriate.
    """
    # Negate coeff if table_z[1] == True
    if (
        len(table_z) > SIGN_FLIP_SECOND_QUBIT_INDEX
        and table_z[SIGN_FLIP_SECOND_QUBIT_INDEX]
    ):
        coeff = -1 * coeff
    # Negate coeff if index 5 == True and no side_chain
    if (not has_side_chain_second_bead) and (
        len(table_z) > SIGN_FLIP_SIXTH_QUBIT_INDEX + 1
        and table_z[SIGN_FLIP_SIXTH_QUBIT_INDEX]
    ):
        coeff = -1 * coeff
    return coeff


def _preset_binary_vals(
    table_z: np.ndarray, *, has_side_chain_second_bead: bool
) -> None:
    main_beads_indices = [0, 1, 2, 3]
    if not has_side_chain_second_bead:
        main_beads_indices.append(5)
    for index in main_beads_indices:
        _preset_single_binary_val(table_z, index)


def _preset_single_binary_val(table_z: np.ndarray, index: int) -> None:
    with contextlib.suppress(IndexError):
        table_z[index] = False
