"""
Defines the Protein class, representing proteins with main and side chains,
including sequence validation and chain initialization.
"""

from exceptions import ChainLengthError
from logger import get_logger
from protein.chain import MainChain, SideChain

logger = get_logger()


class Protein:
    """
    Represents a protein with main and side chains.

    The main chain is defined as a sequence of amino acid residues.
    Side chains are optional and can be empty; if present, they follow the same residue types.
    Side chains cannot be attached to the first or last main bead.

    Attributes:
        main_chain (MainChain): The main chain of the protein.
        side_chain (SideChain): The optional side chain of the protein.

    """

    def __init__(self, main_protein_sequence: str, side_protein_sequence: str) -> None:
        """
        Initialize a Protein instance with given main and side chain sequences.

        Validates that both sequences have the same length and initializes
        the corresponding chain objects.

        Args:
            main_protein_sequence (str): Sequence of residues for the main chain.
            side_protein_sequence (str): Sequence of residues for the side chain.

        Raises:
            ChainLengthError: If main and side chain sequences are not of the same length.

        """
        if len(main_protein_sequence) != len(side_protein_sequence):
            msg = "Main and side protein sequences must be of the same length."
            logger.error(msg)
            raise ChainLengthError(msg)

        self.main_chain: MainChain = MainChain(main_protein_sequence)
        self.side_chain: SideChain = SideChain(side_protein_sequence)
