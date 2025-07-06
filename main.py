from encoding import TetrahedralLattice
from protein import create_protein_sequence, all_turn_combinations



def main() -> None:
    main_chain = [
        "H",
        "H", 
        "P", 
        "P", 
        "H", 
        "P", 
        "H", 
        "P", 
        "H",
    ]

    beads = create_protein_sequence(main_chain)
    
    print(f"\nProtein sequence: {main_chain}")
    print("\nBead information:")
    for bead in beads:
        print(f"  {bead}")
    
    
    lattice = TetrahedralLattice()
    lattice.generate_lattice(
        nx=2, 
        ny=2, 
        nz=2
    )
    
    all_turns = list(all_turn_combinations(len(main_chain)))
    print(f"\nAll possible conformations: {len(all_turns)}")
    
    example_turns = all_turns[0]
    print(f"\nExample turn sequence: {example_turns}")
    
    protein_path = lattice.generate_protein_path(beads, example_turns)
    print(f"Wygenerowana ścieżka białka: {len(protein_path)} pozycji")
    
    lattice.visualize_lattice(
        show_bonds=True, 
        show_node_labels=False, 
        protein_path=protein_path, 
        protein_sequence=main_chain
    )
    

if __name__ == "__main__":
    main()
