from qiskit.quantum_info import SparsePauliOp

from logger import get_logger
from protein.bead import Bead

logger = get_logger()


class PlaceholderSideBead(Bead):
    def __init__(self, symbol: str, index: int, parent_chain_len: int) -> None:
        self.symbol = symbol
        self.index = index
        self.parent_chain_len = parent_chain_len

    def turn_0(self) -> SparsePauliOp:
        msg: str = "Placeholder side bead has no turn qubit operators."
        logger.error(msg)
        raise NotImplementedError(msg)

    def turn_1(self) -> SparsePauliOp:
        msg: str = "Placeholder side bead has no turn qubit operators."
        logger.error(msg)
        raise NotImplementedError(msg)

    def turn_2(self) -> SparsePauliOp:
        msg: str = "Placeholder side bead has no turn qubit operators."
        logger.error(msg)
        raise NotImplementedError(msg)

    def turn_3(self) -> SparsePauliOp:
        msg: str = "Placeholder side bead has no turn qubit operators."
        logger.error(msg)
        raise NotImplementedError(msg)
