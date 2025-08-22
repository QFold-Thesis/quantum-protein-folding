from protein.chain import MainChain, SideChain


class Protein:
    def __init__(self, main_protein_sequence: str, side_protein_sequence: str) -> None:
        if len(main_protein_sequence) != len(side_protein_sequence):
            raise ValueError(
                "Main and side protein sequences must be of the same length."
            )

        self.main_chain: MainChain = MainChain(main_protein_sequence)
        self.side_chain: SideChain = SideChain(side_protein_sequence)
