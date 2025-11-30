from qiskit.quantum_info import SparsePauliOp

from logger import get_logger
from protein.bead import Bead

logger = get_logger()


class _PlaceholderSideBead(Bead):
    """Represents a empty side bead - not attached to a protein's main chain, used as a placeholder for storing empty symbols and padding.

    Attributes:
        symbol (str): Empty symbol.
        index (int): Position of the bead in the side chain.
        parent_chain_len (int): Total number of beads in the parent chain.

    """

    def __init__(self, symbol: str, index: int, parent_chain_len: int) -> None:
        """Initialize the placeholder side bead.

        Note:
            Placeholder side beads do not have turn qubit operators. Any attempt to access them will raise NotImplementedError.
            Symbol of the placeholder bead is set to EMPTY_SIDECHAIN_PLACEHOLDER constant, as initialized in SideChain.
            The "parent_chain_len" arg is used for consistency with other bead types, although it has no functional impact on placeholder beads.

        Args:
            symbol (str): Empty symbol.
            index (int): Position of the bead in the side chain.
            parent_chain_len (int): Total number of beads in the parent chain.

        """
        self.symbol = symbol
        self.index = index
        self.parent_chain_len = parent_chain_len

    def turn_0(self) -> SparsePauliOp:
        """Return the Pauli operator for a turn in direction 0.

        Raises:
            NotImplementedError: Always raised since placeholder side beads do not have turn qubits.

        """
        msg: str = "Placeholder side bead has no turn qubit operators."
        logger.error(msg)
        raise NotImplementedError(msg)

    def turn_1(self) -> SparsePauliOp:
        """Return the Pauli operator for a turn in direction 1.

        Raises:
            NotImplementedError: Always raised since placeholder side beads do not have turn qubits.

        """
        msg: str = "Placeholder side bead has no turn qubit operators."
        logger.error(msg)
        raise NotImplementedError(msg)

    def turn_2(self) -> SparsePauliOp:
        """Return the Pauli operator for a turn in direction 2.

        Raises:
            NotImplementedError: Always raised since placeholder side beads do not have turn qubits.

        """
        msg: str = "Placeholder side bead has no turn qubit operators."
        logger.error(msg)
        raise NotImplementedError(msg)

    def turn_3(self) -> SparsePauliOp:
        """Return the Pauli operator for a turn in direction 3.

        Raises:
            NotImplementedError: Always raised since placeholder side beads do not have turn qubits.

        """
        msg: str = "Placeholder side bead has no turn qubit operators."
        logger.error(msg)
        raise NotImplementedError(msg)
