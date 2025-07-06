from encoding import TetrahedralLattice
from protein import create_protein_sequence, all_turn_combinations


def main() -> None:
    main_chain = ["H", "P", "P", "H", "P", "H", "H", "P"]

    beads = create_protein_sequence(main_chain)

    print(f"\nProtein sequence: {main_chain}")
    print("\nBead information:")
    for bead in beads:
        print(f"  {bead}")

    lattice = TetrahedralLattice()
    lattice.generate_lattice(nx=3, ny=3, nz=3)

    all_turns = list(all_turn_combinations(len(main_chain)))
    print(f"\nAll possible conformations: {len(all_turns)}")

    folding_result = lattice.find_lowest_energy_conformation(beads, all_turns)

    if folding_result["best_turns"] is None:
        print("No valid conformation found.")
        return

    best_positions = folding_result["best_positions"]
    best_turns = folding_result["best_turns"]

    print("\nBest conformation positions:")
    for i, (pos, bead) in enumerate(zip(best_positions, beads)):
        print(f"  {bead.symbol}{i}: {pos}")

    print(f"\nBest turn sequence: {best_turns}")
    print(f"Best energy: {folding_result['best_energy']}")

    lattice.visualize_lattice(
        show_bonds=True,
        show_node_labels=False,
        protein_path=best_positions,
        protein_sequence=main_chain,
    )


if __name__ == "__main__":
    main()
