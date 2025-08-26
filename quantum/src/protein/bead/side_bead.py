from qiskit.quantum_info import SparsePauliOp  # pyright: ignore[reportMissingTypeStubs]

from protein.bead import Bead


class SideBead(Bead):
    def __init__(self, symbol: str, index: int, chain_length: int) -> None:
        super().__init__(symbol, index, chain_length)

    def turn_0(self) -> SparsePauliOp:
        raise NotImplementedError("Side beads do not have turn qubits.")

    def turn_1(self) -> SparsePauliOp:
        raise NotImplementedError("Side beads do not have turn qubits.")

    def turn_2(self) -> SparsePauliOp:
        raise NotImplementedError("Side beads do not have turn qubits.")

    def turn_3(self) -> SparsePauliOp:
        raise NotImplementedError("Side beads do not have turn qubits.")
