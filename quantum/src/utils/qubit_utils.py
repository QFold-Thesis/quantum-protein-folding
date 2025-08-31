import numpy as np
from qiskit.quantum_info import SparsePauliOp  # pyright: ignore[reportMissingTypeStubs]

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
