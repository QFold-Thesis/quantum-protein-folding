import numpy as np
from typing import Union
from qiskit.quantum_info import SparsePauliOp, Pauli  # pyright: ignore[reportMissingTypeStubs]

from constants import NORM_FACTOR


def build_full_identity(num_qubits: int) -> SparsePauliOp:
    """
    Creates a full identity operator of length num_qubits.
    """
    identity_str = "I" * num_qubits
    coeffs = np.array([1.0])
    return SparsePauliOp([identity_str], coeffs)


def build_turn_qubit(num_qubits: int, z_index: int) -> SparsePauliOp:
    """
    Builds a SparsePauliOp representing a turn qubit:
    0.5 * I^n - 0.5 * Z_on_index

    Args:
        num_qubits: total number of qubits
        z_index: index of qubit where Pauli Z acts

    Returns:
        SparsePauliOp representing the turn qubit operator

    """
    identity_op = build_full_identity(num_qubits)

    # Pauli Z on selected index
    pauli_list = ["I"] * num_qubits
    pauli_list[z_index] = "Z"
    pauli_list = pauli_list[::-1]  # Reverse for Qiskit ordering
    pauli_z_str = "".join(pauli_list)
    pauli_z_op = SparsePauliOp([pauli_z_str], np.array([1.0]))

    return NORM_FACTOR * identity_op - NORM_FACTOR * pauli_z_op

def fix_qubits(
    operator: Union[int, SparsePauliOp, Pauli],
    has_side_chain_second_bead: bool = False,
) -> Union[int, SparsePauliOp, Pauli]:
    """
    Assigns predefined values for turns qubits on positions 0, 1, 2, 3, 5 in the main chain.
    Qubits on these position are considered fixed and not subject to optimization.
    """
    # return if operator is int (might be 0 because it is initialized as operator = 0)
    if not isinstance(operator, (SparsePauliOp, Pauli)):
        return operator

    # Normalize operators to SparsePauliOp
    if isinstance(operator, Pauli):
        operator = SparsePauliOp([operator], [1.0])

    new_paulis = []
    new_coeffs = []

    for idx, pauli in enumerate(operator.paulis):
        table_z = np.copy(pauli.z)
        table_x = np.copy(pauli.x)
        coeff = operator.coeffs[idx]

        _preset_binary_vals(table_z, has_side_chain_second_bead)
        coeff = _calc_updated_coeffs(table_z, coeff, has_side_chain_second_bead)
        # Create new Pauli from updated z, x
        new_pauli = Pauli((table_z, table_x))
        new_paulis.append(new_pauli)
        new_coeffs.append(coeff)

    operator_updated = SparsePauliOp(new_paulis, new_coeffs).simplify()
    return operator_updated

def _calc_updated_coeffs(table_z, coeff, has_side_chain_second_bead: bool) -> float:
    """
    Update coefficients based on fixed qubit positions. Negate if appropriate.
    """
    # Negate coeff if table_z[3] == True
    if len(table_z) > 1 and table_z[3]:
        coeff = -1 * coeff
    # Negate coeff if index 5 == True and no side_chain
    if (not has_side_chain_second_bead) and (len(table_z) > 6 and table_z[4]):
        coeff = -1 * coeff
    return coeff

def _preset_binary_vals(table_z: np.ndarray, has_side_chain_second_bead: bool) -> None:
    main_beads_indices = [0, 1, 2, 3]
    if not has_side_chain_second_bead:
        main_beads_indices.append(5)
    for index in main_beads_indices:
        _preset_single_binary_val(table_z, index)

def _preset_single_binary_val(table_z: np.ndarray, index: int) -> None:
    try:
        table_z[index] = False
    except IndexError:
        pass
