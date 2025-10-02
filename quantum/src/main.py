from qiskit.circuit.library import real_amplitudes
from qiskit.primitives import StatevectorSampler as Sampler
from qiskit_algorithms import SamplingVQE
from qiskit_algorithms.optimizers import COBYLA

from builder import HamiltonianBuilder
from constants import EMPTY_SIDECHAIN_PLACEHOLDER
from contact import ContactMap
from distance import DistanceMap
from interaction import HPInteraction, Interaction, MJInteraction
from logger import get_logger
from protein import Protein
from utils.qubit_utils import remove_unused_qubits

logger = get_logger()


def setup_folding_system(main_chain: str, side_chain: str) -> tuple[Protein, MJInteraction, ContactMap, DistanceMap]:
    """Setup the protein folding system components."""
    protein = Protein(
        main_protein_sequence=main_chain, side_protein_sequence=side_chain
    )
    mj_interaction = MJInteraction()
    hp_interaction = HPInteraction()  # noqa: F841

    contact_map = ContactMap(protein=protein)
    distance_map = DistanceMap(protein=protein)

    return protein, mj_interaction, contact_map, distance_map


def build_and_compress_hamiltonian(
    protein: Protein,
    interaction: Interaction,
    contact_map: ContactMap,
    distance_map: DistanceMap
) -> tuple[object, object]:
    """Build and compress the Hamiltonian."""
    h_builder = HamiltonianBuilder(
        protein=protein,
        interaction=interaction,
        distance_map=distance_map,
        contact_map=contact_map,
    )

    hamiltonian = h_builder.sum_hamiltonians()
    logger.debug("Original hamiltonian qubits: %d", hamiltonian.num_qubits)

    compressed_h = remove_unused_qubits(hamiltonian)
    logger.debug("Compressed hamiltonian qubits: %d", compressed_h.num_qubits)

    return hamiltonian, compressed_h


def setup_vqe_optimization(compressed_h: object) -> tuple[object, list, list]:
    """Setup VQE optimization components."""
    optimizer = COBYLA(maxiter=50)
    ansatz = real_amplitudes(num_qubits=compressed_h.num_qubits, reps=1)

    counts = []
    values = []

    def _store_intermediate_result(eval_count, parameters, mean, std):  # noqa: ARG001
        counts.append(eval_count)
        values.append(mean)

    vqe = SamplingVQE(
        sampler=Sampler(),
        ansatz=ansatz,
        optimizer=optimizer,
        aggregation=0.1,
        callback=_store_intermediate_result,
    )

    return vqe, counts, values


def main() -> None:
    main_chain: str = "APRLRFY"
    side_chain: str = EMPTY_SIDECHAIN_PLACEHOLDER * len(main_chain)

    protein, interaction, contact_map, distance_map = setup_folding_system(main_chain, side_chain)

    hamiltonian, compressed_h = build_and_compress_hamiltonian(
        protein=protein,
        interaction=interaction,
        contact_map=contact_map,
        distance_map=distance_map
    )

    vqe, counts, values = setup_vqe_optimization(compressed_h)

    logger.debug("Starting VQE optimization...")
    raw_result = vqe.compute_minimum_eigenvalue(compressed_h)

    logger.debug(raw_result)


if __name__ == "__main__":
    main()
