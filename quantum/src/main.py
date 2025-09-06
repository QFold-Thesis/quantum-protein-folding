from constants import EMPTY_SIDECHAIN_PLACEHOLDER
from interaction import MJInteraction
from protein import Protein
from contact import ContactMap


def main() -> None:
    main_chain = "APRLRFY"
    side_chain = EMPTY_SIDECHAIN_PLACEHOLDER * len(main_chain)

    protein = Protein(
        main_protein_sequence=main_chain, 
        side_protein_sequence=side_chain
    )
    
    mj_interaction = MJInteraction(
        protein=protein
    )

    contact_map = ContactMap(
        protein=protein
    )



if __name__ == "__main__":
    main()
