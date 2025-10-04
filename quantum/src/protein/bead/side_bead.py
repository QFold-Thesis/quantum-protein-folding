from qiskit.quantum_info import SparsePauliOp

from logger import get_logger
from protein.bead import Bead

logger = get_logger()


class SideBead(Bead):
    def __init__(self, symbol: str, index: int, parent_chain_len: int) -> None:
        super().__init__(symbol=symbol, index=index, parent_chain_len=parent_chain_len)

    def turn_0(self) -> SparsePauliOp:
        return (
            (
                (self._full_identity - self.turn_qubits[0])
                @ (self._full_identity - self.turn_qubits[1])
            )
            ^ self._full_identity
        ).simplify()

    def turn_1(self) -> SparsePauliOp:
        return (
            (self.turn_qubits[1] @ (self.turn_qubits[1] - self.turn_qubits[0]))
            ^ self._full_identity
        ).simplify()

    def turn_2(self) -> SparsePauliOp:
        return (
            (self.turn_qubits[0] @ (self.turn_qubits[0] - self.turn_qubits[1]))
            ^ self._full_identity
        ).simplify()

    def turn_3(self) -> SparsePauliOp:
        return (
            self.turn_qubits[0] @ self.turn_qubits[1] ^ self._full_identity
        ).simplify()
