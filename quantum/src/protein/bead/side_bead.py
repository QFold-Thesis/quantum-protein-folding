from qiskit.quantum_info import SparsePauliOp  # pyright: ignore[reportMissingTypeStubs]

from protein.bead import Bead
from logger import get_logger


logger = get_logger()


class SideBead(Bead):
    def __init__(self, symbol: str, index: int, parent_chain_len: int) -> None:
        super().__init__(symbol=symbol, index=index, parent_chain_len=parent_chain_len)

    def turn_0(self) -> SparsePauliOp:
        _msg: str = "Side beads are not yet implemented!"
        logger.error(_msg)
        raise NotImplementedError(_msg)

    def turn_1(self) -> SparsePauliOp:
        _msg: str = "Side beads are not yet implemented!"
        logger.error(_msg)
        raise NotImplementedError(_msg)

    def turn_2(self) -> SparsePauliOp:
        _msg: str = "Side beads are not yet implemented!"
        logger.error(_msg)
        raise NotImplementedError(_msg)

    def turn_3(self) -> SparsePauliOp:
        _msg: str = "Side beads are not yet implemented!"
        logger.error(_msg)
        raise NotImplementedError(_msg)
