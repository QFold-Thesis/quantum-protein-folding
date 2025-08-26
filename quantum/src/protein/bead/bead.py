from abc import ABC, abstractmethod

from qiskit.quantum_info import (  # pyright: ignore[reportMissingTypeStubs]
    Operator,
    SparsePauliOp,
)


class Bead(ABC):
    def __init__(self, symbol: str, index: int) -> None:
        self.symbol: str = symbol
        self.index: int = index
        self.turn_qubits: (
            tuple[SparsePauliOp, SparsePauliOp]
            | tuple[SparsePauliOp, SparsePauliOp, SparsePauliOp, SparsePauliOp]
        )

    @abstractmethod
    def turn_0(self) -> Operator:
        pass

    @abstractmethod
    def turn_1(self) -> Operator:
        pass

    @abstractmethod
    def turn_2(self) -> Operator:
        pass

    @abstractmethod
    def turn_3(self) -> Operator:
        pass
