from qiskit.quantum_info import Pauli # pyright: ignore[reportMissingTypeStubs]

def build_full_identity(num_qubits: int) -> Pauli:
    return Pauli(num_qubits * "I")