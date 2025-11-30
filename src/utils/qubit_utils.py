"""
Module providing tools to create and manipulate Pauli operators,
manage main chain bead states, and prepare operators for qubit-based
protein folding simulations.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
from qiskit.quantum_info import (
    Pauli,
    SparsePauliOp,
)

from constants import (
    IDENTITY_OP_COEFF,
    MAIN_CHAIN_FIFTH_FIXED_POSITION,
    MAIN_CHAIN_FIXED_POSITIONS,
    NORM_FACTOR,
    SIGN_FLIP_SECOND_QUBIT_INDEX,
    SIGN_FLIP_SIXTH_QUBIT_INDEX,
)
from exceptions import InvalidOperatorError
from logger import get_logger

logger = get_logger()

if TYPE_CHECKING:
    from numpy.typing import NDArray


def build_identity_op(
    num_qubits: int, coeff: float = IDENTITY_OP_COEFF
) -> SparsePauliOp:
    """
    Builds a full identity Pauli operator for a given number of qubits.

    Args:
        num_qubits (int): Number of qubits in the operator.
        coeff (float, optional): Coefficient for the operator. Defaults to IDENTITY_OP_COEFF.

    Returns:
        SparsePauliOp: Identity Pauli operator.

    """
    identity_string: str = "I" * num_qubits
    return SparsePauliOp.from_list([(identity_string, coeff)])


def build_turn_qubit(z_index: int, num_qubits: int) -> SparsePauliOp:
    """
    Builds a turn qubit Pauli operator with Z at the specified index.

    Args:
        z_index (int): Index of the qubit to place a Z operator.
        num_qubits (int): Total number of qubits.

    Returns:
        SparsePauliOp: Pauli operator representing the turn qubit.

    """
    z_operator: SparsePauliOp = SparsePauliOp.from_sparse_list(
        [("Z", [z_index], 1.0)], num_qubits=num_qubits
    )

    full_identity: SparsePauliOp = build_identity_op(num_qubits=num_qubits)

    return NORM_FACTOR * (full_identity - z_operator)


def build_pauli_z_operator(num_qubits: int, pauli_z_indices: set[int]) -> SparsePauliOp:
    """
    Build a Pauli operator with Z operators at specified positions and I elsewhere.

    Args:
        num_qubits (int): Total number of qubits.
        pauli_z_indices (set[int]): Indices where Z operators should be applied.

    Returns:
        SparsePauliOp: Constructed Pauli operator.

    """
    if not pauli_z_indices:
        return SparsePauliOp.from_list([("I" * num_qubits, 1.0)])

    idx_sorted = sorted(pauli_z_indices)
    local_label = "Z" * len(idx_sorted)

    return SparsePauliOp.from_sparse_list(
        [(local_label, idx_sorted, 1.0)],
        num_qubits=num_qubits,
    )


def convert_to_qubits(pauli_op: SparsePauliOp) -> SparsePauliOp:
    """
    Convert a Pauli operator to a qubit operator using the identity and normalization factor.

    Args:
        pauli_op (SparsePauliOp): Pauli operator to convert.

    Returns:
        SparsePauliOp: Converted operator.

    Raises:
        InvalidOperatorError: If the Pauli operator does not have a defined number of qubits.

    """
    if pauli_op.num_qubits is None:
        msg: str = "pauli_op.num_qubits is None, cannot convert to qubits."
        raise InvalidOperatorError(msg)

    num_qubits: int = int(pauli_op.num_qubits)
    full_id: SparsePauliOp = SparsePauliOp.from_list([("I" * num_qubits, 1.0)])

    return NORM_FACTOR * (full_id - pauli_op)


def fix_qubits(
    operator: SparsePauliOp,
    *,
    has_side_chain_second_bead: bool = False,
) -> SparsePauliOp:
    """
    Fixes specific qubits in a SparsePauliOp to predefined values for main chain turns.

    Qubits at positions 0, 1, 2, 3, and 5 correspond to fixed turn positions in the main chain
    and are not subject to optimization.

    Args:
        operator (SparsePauliOp): Operator to fix.
        has_side_chain_second_bead (bool, optional): Whether second bead of side chain exists. Defaults to False.

    Returns:
        SparsePauliOp: Operator with fixed qubits.

    """
    # Normalize operators to SparsePauliOp
    coeffs: NDArray[np.float64] = np.array([1.0])
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
            data=[Pauli((table_z, table_x))], coeffs=np.array([operator.coeffs[0]])
        ).simplify()

    new_paulis: list[Pauli] = []
    new_coeffs: NDArray[np.float64] = np.array([])

    for idx, pauli in enumerate(operator.paulis):
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
    table_z: NDArray[np.bool], coeff: float, *, has_side_chain_second_bead: bool
) -> float:
    """
    Update coefficients based on fixed qubit positions.

    Args:
        table_z (NDArray[np.bool]): Z values for each qubit.
        coeff (float): Original coefficient.
        has_side_chain_second_bead (bool): Whether second bead of side chain exists.

    Returns:
        float: Updated coefficient.

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
    table_z: NDArray[np.bool], *, has_side_chain_second_bead: bool
) -> None:
    """
    Set False for main bead indices in the Z table.

    Args:
        table_z (NDArray[np.bool]): Z values for each qubit.
        has_side_chain_second_bead (bool): Whether second bead of side chain exists.

    """
    main_beads_indices = MAIN_CHAIN_FIXED_POSITIONS.copy()

    if not has_side_chain_second_bead:
        main_beads_indices.append(MAIN_CHAIN_FIFTH_FIXED_POSITION)

    for index in main_beads_indices:
        _preset_single_binary_val(table_z, index)


def _preset_single_binary_val(table_z: NDArray[np.bool], index: int) -> None:
    """
    Set a single qubit value to False.

    Args:
        table_z (NDArray[np.bool]): Z values for each qubit.
        index (int): Qubit index to set.

    """
    if index < len(table_z):
        table_z[index] = False


def pad_to_n_qubits(op: SparsePauliOp, target: int) -> SparsePauliOp:
    """
    Extends a Pauli operator with identity qubits to reach the target size.

    Args:
        op (SparsePauliOp): Operator to pad.
        target (int): Target number of qubits.

    Returns:
        SparsePauliOp: Padded operator.

    Raises:
        InvalidOperatorError: If op.num_qubits is None.

    """
    if op.num_qubits is None:
        msg: str = "op.num_qubits is None, cannot pad operator."
        raise InvalidOperatorError(msg)

    if op.num_qubits == target:
        return op
    pad = target - op.num_qubits
    id_pad = build_identity_op(pad)
    logger.debug("Padding operator from %s to %s qubits.", op.num_qubits, target)
    return id_pad ^ op


def find_unused_qubits(op: SparsePauliOp) -> list[int]:
    """
    Return indices of qubits that are identity (I) in every term of the operator.

    Args:
        op (SparsePauliOp): Operator to check.

    Returns:
        list[int]: List of unused qubit indices.

    Raises:
        InvalidOperatorError: If op.num_qubits is None.

    """
    if op.num_qubits is None:
        msg: str = "op.num_qubits is None, cannot find unused qubits."
        raise InvalidOperatorError(msg)

    if op.num_qubits == 0 or len(op.paulis) == 0:
        return []
    used_mask = np.any(op.paulis.z, axis=0)
    return [i for i, used in enumerate(used_mask) if not used]


def remove_unused_qubits(
    op: SparsePauliOp,
) -> SparsePauliOp:
    """
    Remove qubits that are identity in all terms.

    Args:
        op (SparsePauliOp): Operator to remove unused qubits from.

    Returns:
        SparsePauliOp: Operator without unused qubits.

    Raises:
        InvalidOperatorError: If op.num_qubits is None.

    """
    if op.num_qubits is None:
        msg: str = "op.num_qubits is None, cannot remove unused qubits."
        raise InvalidOperatorError(msg)

    unused = find_unused_qubits(op)
    if not unused:
        return op.copy()

    mask = np.ones(op.num_qubits, dtype=bool)
    mask[unused] = False

    z_full = op.paulis.z
    x_full = op.paulis.x
    new_paulis: list[Pauli] = [
        Pauli((z_full[k][mask], x_full[k][mask])) for k in range(len(op.paulis))
    ]

    return SparsePauliOp(new_paulis, coeffs=op.coeffs).simplify()
