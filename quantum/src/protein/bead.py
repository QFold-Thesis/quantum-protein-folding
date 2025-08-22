from qiskit.quantum_info import SparsePauliOp


class Bead:
    def __init__(self, symbol: str, index: int) -> None:
        self.symbol: str = symbol
        self.index: int = index
        self.turn_qubits: (
            tuple[SparsePauliOp, SparsePauliOp]
            | tuple[SparsePauliOp, SparsePauliOp, SparsePauliOp, SparsePauliOp]
        )
