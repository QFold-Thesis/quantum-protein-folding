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

    The main chain is defined as a sequence of amino acid residues, with valid types
    [A, C, D, E, F, G, H, I, K, L, M, N, P, Q, R, S, T, V, W, Y].
    Side chains are optional and can be empty; if present, they follow the same residue types.
    Side chains cannot be attached to the first or last main bead.
    """

    def __init__(self, main_protein_sequence: str, side_protein_sequence: str) -> None:
        """
        Initializes a Protein instance with given main and side chain sequences.
        """
        if len(main_protein_sequence) != len(side_protein_sequence):
            msg = "Main and side protein sequences must be of the same length."
            logger.error(msg)
            raise ChainLengthError(msg)

        self.main_chain: MainChain = MainChain(main_protein_sequence)
        self.side_chain: SideChain = SideChain(side_protein_sequence)
