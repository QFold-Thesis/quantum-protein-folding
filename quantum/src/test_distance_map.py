from builder import HamiltonianBuilder
from protein import Protein


def test_first_neighbor():
    """
    Tests that Pauli operators for 1st neighbor interactions are created correctly.
    """
    main_chain_residue_seq = "SAASSS"
    side_chain_residue_seq = "_" * len(main_chain_residue_seq)

    protein = Protein(
        main_protein_sequence=main_chain_residue_seq,
        side_protein_sequence=side_chain_residue_seq,
    )
    lambda_1 = 2
    lower_main_bead_index = 0
    upper_main_bead_index = 3

    hamiltonian_builder = HamiltonianBuilder(protein)
    expr = hamiltonian_builder.get_first_neighbor_hamiltonian(
        lower_main_bead_index,
        upper_main_bead_index,
        lambda_1,
    )

    print(expr)  # noqa: T201


def test_second_neighbor():
    """
    Tests that Pauli operators for 2nd neighbor interactions are created correctly.
    """
    main_chain_residue_seq = "SAASS"
    side_chain_residue_sequences = ["", "", "", "", ""]

    protein = Protein(
        main_protein_sequence=main_chain_residue_seq,
        side_protein_sequence=side_chain_residue_sequences,
    )
    lambda_1 = 2
    lower_main_bead_index = 0
    upper_main_bead_index = 3

    hamiltonian_builder = HamiltonianBuilder(protein)
    expr = hamiltonian_builder.get_second_neighbor_hamiltonian(
        lower_main_bead_index,
        upper_main_bead_index,
        lambda_1,
    )
    print(f"result: {expr}")  # noqa: T201


if __name__ == "__main__":
    test_first_neighbor()
    test_second_neighbor()
