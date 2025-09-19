from constants import EMPTY_SIDECHAIN_PLACEHOLDER
from contact import ContactMap
from distance import DistanceMap
from interaction import MJInteraction
from protein import Protein


def main() -> None:
    main_chain = "APRLRFY"
    side_chain = EMPTY_SIDECHAIN_PLACEHOLDER * len(main_chain)

    protein = Protein(
        main_protein_sequence=main_chain, side_protein_sequence=side_chain
    )

    _ = MJInteraction(protein=protein)

    _ = ContactMap(protein=protein)

    _ = DistanceMap(protein=protein)


if __name__ == "__main__":
    main()
