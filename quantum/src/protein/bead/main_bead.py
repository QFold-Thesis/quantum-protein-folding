from qiskit.quantum_info import SparsePauliOp  # pyright: ignore[reportMissingTypeStubs]

from protein.bead import Bead
from quantum.src.constants import CONFORMATION_ENCODING
from enums import ConformationEncoding

class MainBead(Bead):
    def __init__(self, symbol: str, index: int, chain_length: int) -> None:
        super().__init__(symbol, index, chain_length)

    def turn_0(self) -> SparsePauliOp:
        if CONFORMATION_ENCODING == ConformationEncoding.SPARSE:
            return self.turn_qubits[0]
        if CONFORMATION_ENCODING == ConformationEncoding.DENSE:
            # full_identity = build_full_identity(2 * (self.chain_length - 1))
            # return full_identity ^ ((full_identity - self.turn_qubits[0]) @ (full_identity - self.turn_qubits[1]))
            pass
        raise NotImplementedError("Error.")

    def turn_1(self) -> SparsePauliOp:
        if CONFORMATION_ENCODING == ConformationEncoding.SPARSE:
            return self.turn_qubits[1]
        if CONFORMATION_ENCODING == ConformationEncoding.DENSE:
            pass
        raise NotImplementedError("Error.")

    def turn_2(self) -> SparsePauliOp:
        if CONFORMATION_ENCODING == ConformationEncoding.SPARSE:
            return self.turn_qubits[2]
        if CONFORMATION_ENCODING == ConformationEncoding.DENSE:
            pass
        raise NotImplementedError("Error.")

    def turn_3(self) -> SparsePauliOp:
        if CONFORMATION_ENCODING == ConformationEncoding.SPARSE:
            return self.turn_qubits[3]
        if CONFORMATION_ENCODING == ConformationEncoding.DENSE:
            pass
        raise NotImplementedError("Error.")
