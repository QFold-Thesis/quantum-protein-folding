from qiskit.quantum_info import Operator

from protein.bead import Bead


class SideBead(Bead):
    def __init__(self, symbol: str, index: int) -> None:
        super().__init__(symbol, index)

    def turn_0(self) -> Operator:
        pass

    def turn_1(self) -> Operator:
        pass

    def turn_2(self) -> Operator:
        pass

    def turn_3(self) -> Operator:
        pass
