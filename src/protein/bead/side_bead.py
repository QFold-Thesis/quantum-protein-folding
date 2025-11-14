from qiskit.quantum_info import SparsePauliOp

from logger import get_logger
from protein.bead import Bead

logger = get_logger()


class _SideBead(Bead):
    """Represents a side bead attached to a protein's main chain."""

    def __init__(self, symbol: str, index: int, parent_chain_len: int) -> None:
        """
        Initialize a side bead with its symbol and position in the chain.

        Args:
            symbol (str): Amino acid symbol representing the side bead.
            index (int): Position of the bead in the chain.
            parent_chain_len (int): Total number of beads in the parent chain.

        """
        _msg: str = "Side beads are not yet implemented!"
        logger.error(_msg)
        raise NotImplementedError(_msg)

    def turn_0(self) -> SparsePauliOp:
        """
        Return the Pauli operator for a turn in direction 0.

        Raises:
            NotImplementedError: Always raised since side bead turns are not implemented.

        """
        _msg: str = "Side beads are not yet implemented!"
        logger.error(_msg)
        raise NotImplementedError(_msg)

    def turn_1(self) -> SparsePauliOp:
        """
        Return the Pauli operator for a turn in direction 1.

        Raises:
            NotImplementedError: Always raised since side bead turns are not implemented.

        """
        _msg: str = "Side beads are not yet implemented!"
        logger.error(_msg)
        raise NotImplementedError(_msg)

    def turn_2(self) -> SparsePauliOp:
        """
        Return the Pauli operator for a turn in direction 2.

        Raises:
            NotImplementedError: Always raised since side bead turns are not implemented.

        """
        _msg: str = "Side beads are not yet implemented!"
        logger.error(_msg)
        raise NotImplementedError(_msg)

    def turn_3(self) -> SparsePauliOp:
        """
        Return the Pauli operator for a turn in direction 3.

        Raises:
            NotImplementedError: Always raised since side bead turns are not implemented.

        """
        _msg: str = "Side beads are not yet implemented!"
        logger.error(_msg)
        raise NotImplementedError(_msg)
