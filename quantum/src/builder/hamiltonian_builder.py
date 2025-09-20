import numpy as np

from logger import get_logger
from protein import Protein

logger = get_logger()


class HamiltonianBuilder:
    def __init__(self, protein: Protein):
        ## Add pair_energies and penalties to args.
        self.protein = protein

    def get_first_neighbor_hamiltonian(
        self,
        lower_bead_idx: int,
        upper_bead_idx: int,
        lambda_1: float,
        pair_energies: np.ndarray,
    ):
        pass

    def get_second_neighbor_hamiltonian(
        self,
        lower_bead_idx: int,
        upper_bead_idx: int,
        lambda_1: float,
        pair_energies: np.ndarray,
    ):
        pass
