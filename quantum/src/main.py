from interaction.mj_interaction import MJInteraction


def main() -> None:
    main_chain = "APRLRFY"
    side_chain = "_______"
    print(main_chain, side_chain)

    mj_interaction = MJInteraction(protein_sequence=main_chain)
    print(mj_interaction.energy_pairs)


if __name__ == "__main__":
    main()
