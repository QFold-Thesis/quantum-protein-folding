from interaction.mj_interaction import MJInteraction
from protein.protein import Protein
from constants import EMPTY_SIDECHAIN_PLACEHOLDER


def main() -> None:
    main_chain = "APRLRFY"
    side_chain = EMPTY_SIDECHAIN_PLACEHOLDER * len(main_chain)

    Protein(main_chain, side_chain)

    MJInteraction(protein_sequence=main_chain)


if __name__ == "__main__":
    main()
