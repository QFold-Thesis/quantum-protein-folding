from qiskit.quantum_info import SparsePauliOp  # pyright: ignore[reportMissingTypeStubs]


def build_full_identity(num_qubits: int) -> SparsePauliOp:
    """Builds a full identity Pauli operator for a given number of qubits."""
    identity_string = "I" * num_qubits
    return SparsePauliOp.from_list([(identity_string, 1.0)])
