from encoding import TetrahedralLattice, all_turn_combinations
import matplotlib.pyplot as plt
from utils import draw_chain

def aa_to_hp(sequence: str) -> str:
    """
    Konwertuje sekwencję aminokwasów na model HP (H - hydrofobowy, P - polarny).
    Przyjmuje sekwencję jako string z jednoliterowym kodem aminokwasów.
    Zwraca string z literami H i P.
    """
    # Lista hydrofobowych aminokwasów wg najczęściej stosowanej konwencji HP
    hydrophobic = {'A', 'V', 'L', 'I', 'M', 'F', 'W', 'Y', 'C'}
    # Pozostałe uznaje się za polarne
    hp_sequence = ""
    for aa in sequence.upper():
        if aa in hydrophobic:
            hp_sequence += "H"
        else:
            hp_sequence += "P"
    return hp_sequence

def main() -> None:
    main_chain = ["H", "P", "P", "H", "P", "H", "H", "P", ]

    lattice = TetrahedralLattice(30, 30, 30)
    all_turns = all_turn_combinations(len(main_chain))

    min_energy = float('inf')
    min_turns = None
    for turns in all_turns:
        try:
            positions = lattice.generate_positions((1, 1, 1), turns)
            energy = lattice.compute_energy(positions, main_chain)
            if energy < min_energy:
                min_energy = energy
                min_turns = turns
        except ValueError as e:
            continue

    if min_turns is None:
        print("Nie znaleziono żadnej prawidłowej konformacji!")
        return
    
    print(f"Minimalna energia: {min_energy}")

    positions = lattice.generate_positions((1, 1, 1), min_turns)
    draw_chain(positions, main_chain, lattice)


if __name__ == "__main__":
    main()
