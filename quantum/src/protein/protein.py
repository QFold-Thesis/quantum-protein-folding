from exceptions import ChainLengthError
from logger import get_logger
from protein.chain import MainChain

logger = get_logger()


class Protein:
    def __init__(self, main_protein_sequence: str, side_protein_sequences: list[str]) -> None:
        if len(main_protein_sequence) != len(side_protein_sequences):
            msg = "Main and side protein sequences must be of the same length."
            logger.error(msg)
            raise ChainLengthError(msg)

        self.main_chain: MainChain = MainChain(main_protein_sequence, side_protein_sequences)
