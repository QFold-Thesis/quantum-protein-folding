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
    num_qubits: int = int(pauli_op.num_qubits) # pyright: ignore[reportArgumentType]
    full_id: SparsePauliOp = SparsePauliOp.from_list([("I" * num_qubits, 1.0)])

    converted = NORM_FACTOR * (full_id - pauli_op)
    return converted
