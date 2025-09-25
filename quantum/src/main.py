from builder import HamiltonianBuilder
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

    builder = HamiltonianBuilder(protein)
    hamiltonian, backbone, backtrack = builder.sum_hamiltonians()
    print(  # noqa: T201
        40 * "-",
        "\n backbone hamiltonian\n",
        backbone,
    )

    print(  # noqa: T201
        40 * "-",
        "\n backtrack hamiltonian\n",
        backtrack,
    )

    print(  # noqa: T201
        40 * "-",
        "\nfinal hamiltonian\n",
        hamiltonian,
    )

    from qiskit.circuit.library import RealAmplitudes
    from qiskit_algorithms.optimizers import COBYLA
    from qiskit_algorithms import SamplingVQE
    from qiskit.primitives import StatevectorSampler as Sampler

    # Classical optimizer
    optimizer = COBYLA(maxiter=50)

    # Variational ansatz
    ansatz = RealAmplitudes(reps=1)

    counts = []
    values = []

    def store_intermediate_result(eval_count, parameters, mean, std):
        counts.append(eval_count)
        values.append(mean)


    # Inicjalizacja VQE
    vqe = SamplingVQE(
        sampler=Sampler(),
        ansatz=ansatz,
        optimizer=optimizer,
        aggregation=0.1,  # CVaR dla alpha=0.1
        callback=store_intermediate_result
    )

    # Uruchomienie optymalizacji
    raw_result = vqe.compute_minimum_eigenvalue(hamiltonian)
    print(raw_result)



if __name__ == "__main__":
    main()
