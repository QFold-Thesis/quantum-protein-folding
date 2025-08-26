import numpy as np
from qiskit.quantum_info import Operator  # pyright: ignore[reportMissingTypeStubs]

from protein.bead import Bead


class MainBead(Bead):
    def __init__(self, symbol: str, index: int) -> None:
        super().__init__(symbol, index)

    def turn_0(self) -> Operator:
        return Operator(np.eye(4))

    def turn_1(self) -> Operator:
        return Operator(np.eye(4))

    def turn_2(self) -> Operator:
        return Operator(np.eye(4))

    def turn_3(self) -> Operator:
        return Operator(np.eye(4))
