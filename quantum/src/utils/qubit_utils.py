from qiskit.quantum_info import SparsePauliOp  # pyright: ignore[reportMissingTypeStubs]
from constants import NORM_FACTOR

def build_full_identity(num_qubits: int) -> SparsePauliOp:
    """Builds a full identity Pauli operator for a given number of qubits."""
    identity_string: str = "I" * num_qubits
    return SparsePauliOp.from_list([(identity_string, 1.0)])


def build_turn_qubit(z_index: int, num_qubits: int) -> SparsePauliOp:
    """Builds a turn qubit Pauli operator with Z at the specified index."""

    z_operator: SparsePauliOp = SparsePauliOp.from_sparse_list(
        [("Z", [z_index], 1.0)],
        num_qubits=num_qubits
    )

    full_identity = build_full_identity(num_qubits=num_qubits)

    return NORM_FACTOR * (full_identity - z_operator)