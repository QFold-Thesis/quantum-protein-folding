from protein.chain import MainChain, SideChain
from exceptions import ChainLengthError
from logger import get_logger


logger = get_logger()


class Protein:
    def __init__(self, main_protein_sequence: str, side_protein_sequence: str) -> None:
        if len(main_protein_sequence) != len(side_protein_sequence):
            msg = "Main and side protein sequences must be of the same length."
            logger.error(msg)
            raise ChainLengthError(msg)

        self.main_chain: MainChain = MainChain(main_protein_sequence)
        self.side_chain: SideChain = SideChain(side_protein_sequence)
