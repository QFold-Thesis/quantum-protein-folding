from abc import ABC, abstractmethod

from qiskit.quantum_info import (  # pyright: ignore[reportMissingTypeStubs]
    SparsePauliOp,
)

from enums import SubLattice

DenseTurnQubits = tuple[SparsePauliOp, SparsePauliOp]
SparseTurnQubits = tuple[SparsePauliOp, SparsePauliOp, SparsePauliOp, SparsePauliOp]


class Bead(ABC):
    def __init__(
        self,
        symbol: str,
        index: int,
        chain_length: int,
        *,  # enforce keyword arguments
        is_turning: bool,
    ) -> None:
        self.symbol: str = symbol
        self.index: int = index
        self.is_turning: bool = is_turning
        self.turn_qubits: tuple[SparsePauliOp, ...] = ()
        self.sublattice: SubLattice = SubLattice.B if index % 2 == 1 else SubLattice.A
        self.chain_length: int = chain_length

    @abstractmethod
    def _initate_turn_qubits(self) -> None:
        pass

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
