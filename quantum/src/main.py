from interaction.mj_interaction import MJInteraction


def main() -> None:
    main_chain = "APRLRFY"

    mj_interaction = MJInteraction(protein_sequence=main_chain)
    print(mj_interaction.energy_pairs)


if __name__ == "__main__":
    main()
