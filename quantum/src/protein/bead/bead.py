from abc import ABC, abstractmethod

from qiskit.quantum_info import (  # pyright: ignore[reportMissingTypeStubs]
    SparsePauliOp,
)

from enums import SubLattice

class Bead(ABC):
    def __init__(self, symbol: str, index: int, full_identity: SparsePauliOp) -> None:
        self.symbol: str = symbol
        self.index: int = index
        self.turn_qubits: tuple[SparsePauliOp, ...] = tuple()
        self.sublattice: SubLattice = SubLattice.B if index % 2 == 1 else SubLattice.A
        
        self._full_identity: SparsePauliOp = full_identity

    @abstractmethod
    def turn_0(self) -> SparsePauliOp:
        pass

    @abstractmethod
    def turn_1(self) -> SparsePauliOp:
        pass

    @abstractmethod
    def turn_2(self) -> SparsePauliOp:
        pass

    @abstractmethod
    def turn_3(self) -> SparsePauliOp:
        pass
