from qiskit.circuit.library import RealAmplitudes
from qiskit.primitives import StatevectorSampler as Sampler
from qiskit_algorithms import SamplingVQE
from qiskit_algorithms.optimizers import COBYLA

from builder import HamiltonianBuilder
from constants import EMPTY_SIDECHAIN_PLACEHOLDER
from contact import ContactMap
from distance import DistanceMap
from interaction import HPInteraction, MJInteraction  # noqa: F401
from protein import Protein
from utils.qubit_utils import remove_unused_qubits


def main() -> None:
    main_chain = "APRLRFY"
    side_chain = EMPTY_SIDECHAIN_PLACEHOLDER * len(main_chain)

    protein = Protein(
        main_protein_sequence=main_chain, side_protein_sequence=side_chain
    )

    mj_interaction = MJInteraction()
    #hp_interaction = HPInteraction()  # noqa: ERA001

    contact_map = ContactMap(protein=protein)

    distance_map = DistanceMap(protein=protein)

    builder = HamiltonianBuilder(
        protein=protein,
        interaction=mj_interaction,
        distance_map=distance_map,
        contact_map=contact_map,
    )
    hamiltonian = builder.sum_hamiltonians()
    print(  # noqa: T201
        40 * "-",
        "\n hamiltonian\n",
        hamiltonian,
    )

    optimizer = COBYLA(maxiter=50)

    ansatz = RealAmplitudes(reps=1)

    counts = []
    values = []

    def store_intermediate_result(eval_count, parameters, mean, std):  # noqa: ARG001
        counts.append(eval_count)
        values.append(mean)

    vqe = SamplingVQE(
        sampler=Sampler(),
        ansatz=ansatz,
        optimizer=optimizer,
        aggregation=0.1,
        callback=store_intermediate_result,
    )

    compressed_h = remove_unused_qubits(hamiltonian)
    raw_result = vqe.compute_minimum_eigenvalue(compressed_h)
    print(raw_result)  # noqa: T201


if __name__ == "__main__":
    main()
