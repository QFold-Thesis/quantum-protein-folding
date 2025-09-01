from qiskit.quantum_info import SparsePauliOp  # pyright: ignore[reportMissingTypeStubs]
from quantum.src.constants import CONFORMATION_ENCODING

from enums import ConformationEncoding
from exceptions import ConformationEncodingError
from protein.bead import Bead


class MainBead(Bead):
    def __init__(self, symbol: str, index: int, full_identity: SparsePauliOp) -> None:
        super().__init__(symbol, index, full_identity)

    def turn_0(self) -> SparsePauliOp:
        if CONFORMATION_ENCODING == ConformationEncoding.SPARSE:
            return (
                self._full_identity.tensor(self.turn_qubits[0])
            )
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
            return (
                self._full_identity.tensor(self.turn_qubits[1])
            )
        if CONFORMATION_ENCODING == ConformationEncoding.DENSE:
            return (
                self._full_identity.tensor(
                    (self.turn_qubits[1]) @
                    (self.turn_qubits[1] - self.turn_qubits[0])
                )
            ).simplify()
        raise ConformationEncodingError

    def turn_2(self) -> SparsePauliOp:
        if CONFORMATION_ENCODING == ConformationEncoding.SPARSE:
            return (
                self._full_identity.tensor(self.turn_qubits[2])
            )
        if CONFORMATION_ENCODING == ConformationEncoding.DENSE:
            return (
                self._full_identity.tensor(
                    (self.turn_qubits[0]) @
                    (self.turn_qubits[0] - self.turn_qubits[1])
                )
            ).simplify()
        raise ConformationEncodingError

    def turn_3(self) -> SparsePauliOp:
        if CONFORMATION_ENCODING == ConformationEncoding.SPARSE:
            return (
                self._full_identity.tensor(self.turn_qubits[3])
            )
        if CONFORMATION_ENCODING == ConformationEncoding.DENSE:
            return (
                self._full_identity.tensor(
                    (self.turn_qubits[0]) @
                    (self.turn_qubits[1])
                )
            ).simplify()
        raise ConformationEncodingError
