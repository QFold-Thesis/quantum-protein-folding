from protein.chain import Chain


class Protein:
    def __init__(self, main_protein_sequence: str, side_protein_sequence: str) -> None:
        self.main_chain: Chain = Chain(main_protein_sequence)
        self.side_chain: Chain = Chain(side_protein_sequence)
