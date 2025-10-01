from qiskit.quantum_info import SparsePauliOp

from logger import get_logger
from protein.bead import Bead

logger = get_logger()


class SideBead(Bead):
    """
    Represents a side bead attached to a protein's main chain.
    """

    def __init__(self, symbol: str, index: int, parent_chain_len: int) -> None:
        """
        Initializes a side bead with its symbol and position in the chain.
        """
        super().__init__(symbol=symbol, index=index, parent_chain_len=parent_chain_len)

    def turn_0(self) -> SparsePauliOp:
        """
        Turn operator for side bead in direction 0.
        """
        _msg: str = "Side beads are not yet implemented!"
        logger.error(_msg)
        raise NotImplementedError(_msg)

    def turn_1(self) -> SparsePauliOp:
        """
        Turn operator for side bead in direction 1.
        """
        _msg: str = "Side beads are not yet implemented!"
        logger.error(_msg)
        raise NotImplementedError(_msg)

    def turn_2(self) -> SparsePauliOp:
        """
        Turn operator for side bead in direction 2.
        """
        _msg: str = "Side beads are not yet implemented!"
        logger.error(_msg)
        raise NotImplementedError(_msg)

    def turn_3(self) -> SparsePauliOp:
        """
        Turn operator for side bead in direction 3.
        """
        _msg: str = "Side beads are not yet implemented!"
        logger.error(_msg)
        raise NotImplementedError(_msg)
