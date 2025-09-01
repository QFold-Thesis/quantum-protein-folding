from qiskit.quantum_info import SparsePauliOp  # pyright: ignore[reportMissingTypeStubs]

from constants import CONFORMATION_ENCODING
from enums import ConformationEncoding
from exceptions import ConformationEncodingError
from protein.bead import Bead


class MainBead(Bead):
    def __init__(self, symbol: str, index: int, parent_chain_len: int) -> None:
        super().__init__(symbol=symbol, index=index, parent_chain_len=parent_chain_len)

    def turn_0(self) -> SparsePauliOp:
        if CONFORMATION_ENCODING == ConformationEncoding.SPARSE:
            return self.turn_qubits[0]
        if CONFORMATION_ENCODING == ConformationEncoding.DENSE:
            return (
                self._full_identity.tensor(
                    (self._full_identity - self.turn_qubits[0])
                    @ (self._full_identity - self.turn_qubits[1])
                )
            ).simplify()
        raise ConformationEncodingError

    def turn_1(self) -> SparsePauliOp:
        if CONFORMATION_ENCODING == ConformationEncoding.SPARSE:
            return self.turn_qubits[1]
        if CONFORMATION_ENCODING == ConformationEncoding.DENSE:
            return (
                self._full_identity.tensor(
                    (self.turn_qubits[1]) @ (self.turn_qubits[1] - self.turn_qubits[0])
                )
            ).simplify()
        raise ConformationEncodingError

    def turn_2(self) -> SparsePauliOp:
        if CONFORMATION_ENCODING == ConformationEncoding.SPARSE:
            return self.turn_qubits[2]
        if CONFORMATION_ENCODING == ConformationEncoding.DENSE:
            return (
                self._full_identity.tensor(
                    (self.turn_qubits[0]) @ (self.turn_qubits[0] - self.turn_qubits[1])
                )
            ).simplify()
        raise ConformationEncodingError

    def turn_3(self) -> SparsePauliOp:
        if CONFORMATION_ENCODING == ConformationEncoding.SPARSE:
            return self.turn_qubits[3]
        if CONFORMATION_ENCODING == ConformationEncoding.DENSE:
            return (
                self._full_identity.tensor(
                    (self.turn_qubits[0]) @ (self.turn_qubits[1])
                )
            ).simplify()
        raise ConformationEncodingError
