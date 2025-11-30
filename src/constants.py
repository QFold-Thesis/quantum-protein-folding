"""Definitions of constants used throughout the protein folding simulations.

This module centralizes configuration parameters, file paths, numeric constants,
energy coefficients, and backend options required by different components of the
project.

ROOT_PROJECT_PATH : Path
    Root directory of the repository.
LOGS_DIRPATH : Path
    Directory where log files are stored.
MJ_INTERACTION_MATRIX_FILEPATH : Path
    Path to the Miyazawa-Jernigan interaction matrix file.
HP_INTERACTION_MATRIX_FILEPATH : Path
    Path to the hydrophobic-polar matrix file.
DEFAULT_TIMEZONE : tzinfo
    Timezone used for timestamping logs and outputs.
CONFORMATION_ENCODING : ConformationEncoding
    Selected encoding type for representing turns in the lattice.
QUBITS_PER_TURN : int
    Number of qubits required to encode a single turn.
EMPTY_SIDECHAIN_PLACEHOLDER : str
    Placeholder symbol for missing side-chain positions.
NORM_FACTOR : float
    Normalization factor applied during distance vector construction.
SIGN_FLIP_SECOND_QUBIT_INDEX : int
    Index of the second qubit whose sign may require inversion.
SIGN_FLIP_SIXTH_QUBIT_INDEX : int
    Index of the sixth qubit whose sign may require inversion.
DIST_VECTOR_AXES : int
    Number of axes in the distance operator vector.
BOUNDING_CONSTANT : int
    Limits how far the chain is allowed to spread in space, preventing exploration of physically impossible conformations.
MJ_ENERGY_MULTIPLIER : float
    Scaling factor applied to Miyazawa-Jernigan energies.
HP_HH_CONTACT_ENERGY : float
    Energy contribution for hydrophobic-hydrophobic contact in the HP model.
HP_NON_HH_CONTACT_ENERGY : float
    Energy contribution for non-hydrophobic contacts in the HP model.
GLOBAL_LOGGER_NAME : str
    Name of the root logger for the project.
LOGGER_DEFAULT_LEVEL : int
    Default logging verbosity.
MIN_DISTANCE_BETWEEN_CONTACTS : int
   Smallest allowed distance along the chain for beads to form a valid contact.
MAIN_CHAIN_FIXED_POSITIONS : list[int]
    Main-chain bead indices that must be assigned predefined values during encoding.
MAIN_CHAIN_FIFTH_FIXED_POSITION : int
    Index of the bead with a predetermined position in the chain.
INTERACTION_TYPE : InteractionType
    Selected interaction model (MJ or HP).
IDENTITY_OP_COEFF : float
    Coefficient applied to the identity operator when building the Hamiltonian.
EMPTY_OP_COEFF : float
    Coefficient used when an operator contributes no term.
SPARSE_TURN_INDICATORS : dict[TurnDirection, str]
    Sparse encoding for turn directions.
DENSE_TURN_INDICATORS : dict[TurnDirection, str]
    Dense encoding for turn directions.
XYZ_FILE_LINE_START_INDEX : int
    Number of header lines in .xyz coordinate files.
XYZ_FILE_PARTS_PER_LINE : int
    Expected number of fields per line in .xyz files.
RESULTS_DATA_DIRPATH : Path
    Directory where simulation's results are written.
RAW_VQE_RESULTS_FILENAME : str
    Filename for raw VQE results.
XYZ_FILENAME : str
    Name of the output XYZ coordinate file.
SPARSE_VQE_RESULTS_FILENAME : str
    Filename for sparse VQE result summaries.
VQE_ITERATIONS_FILENAME : str
    Filename recording energy per iteration.
GIF_VISUALIZATION_FILENAME : str
    Name of the generated 3D GIF visualization.
HTML_VISUALIZATION_FILENAME : str
    Name of the interactive HTML visualization.
FLAT_VISUALIZATION_FILENAME : str
    Name of the 2D conformational plot.
TETRAHEDRAL_LATTICE_PADDING : int
    Extra spacing added to bead coordinates to prevent index conflicts.
INDEX_COLNAME : str
    Column name for bead indices.
SYMBOL_COLNAME : str
    Column name for residue symbols.
ITERATION_COLNAME : str
    Column name for iteration counters.
COORDINATES_COLUMN_WIDTH : int
    Output formatting width for coordinate columns.
FCC_EDGE_LENGTH : float
    Edge length of fractional FCC unit vectors.
FCC_BASIS : ndarray
    Basis vectors of the face-centered cubic lattice.
SIDE_CHAIN_FIFTH_POSITION_INDEX : int
    Index of the side-chain bead associated with the fifth main-chain position.
BACKEND_TYPE : BackendType
    Selected quantum backend type.
IBM_QUANTUM_TOKEN : str or None
    IBM Quantum API token loaded from environment.
IBM_QUANTUM_BACKEND_NAME : str
    Name of the IBM Quantum device to execute circuits on.
IBM_QUANTUM_SHOTS : int
    Number of measurement shots for hardware execution.
MIN_CHAIN_LENGTH : int
    Minimum allowed protein chain length for simulations.
"""

from __future__ import annotations

import datetime
import logging
import os
from datetime import tzinfo
from pathlib import Path
from typing import TYPE_CHECKING

import numpy as np
from dotenv import load_dotenv

from enums import BackendType, ConformationEncoding, InteractionType, TurnDirection

if TYPE_CHECKING:
    from numpy.typing import NDArray

ROOT_PROJECT_PATH: Path = Path(__file__).parent.parent

load_dotenv(ROOT_PROJECT_PATH / ".env")

LOGS_DIRPATH: Path = ROOT_PROJECT_PATH / "output" / "logs"

MJ_INTERACTION_MATRIX_FILEPATH: Path = (
    ROOT_PROJECT_PATH / "src" / "resources" / "mj_matrix.txt"
)

HP_INTERACTION_MATRIX_FILEPATH: Path = (
    ROOT_PROJECT_PATH / "src" / "resources" / "hp_matrix.txt"
)

DEFAULT_TIMEZONE: tzinfo = datetime.UTC

CONFORMATION_ENCODING: ConformationEncoding = ConformationEncoding.DENSE

QUBITS_PER_TURN: int = CONFORMATION_ENCODING.value

EMPTY_SIDECHAIN_PLACEHOLDER: str = "_"

NORM_FACTOR: float = 0.5

SIGN_FLIP_SECOND_QUBIT_INDEX: int = 1

SIGN_FLIP_SIXTH_QUBIT_INDEX: int = 5

DIST_VECTOR_AXES: int = 4

BOUNDING_CONSTANT: int = 7

MJ_ENERGY_MULTIPLIER: float = 0.1

HP_HH_CONTACT_ENERGY: float = -1.0

HP_NON_HH_CONTACT_ENERGY: float = 0.0

GLOBAL_LOGGER_NAME: str = "global_logger"

LOGGER_DEFAULT_LEVEL: int = logging.DEBUG

MIN_DISTANCE_BETWEEN_CONTACTS: int = (
    5  # Minimum bonds between two beads to consider a contact
)

MAIN_CHAIN_FIXED_POSITIONS: list[int] = [0, 1, 2, 3, 5]
MAIN_CHAIN_FIFTH_FIXED_POSITION: int = 5

INTERACTION_TYPE: InteractionType = InteractionType.MJ

IDENTITY_OP_COEFF: float = 1.0

EMPTY_OP_COEFF: float = 0.0

SPARSE_TURN_INDICATORS: dict[TurnDirection, str] = {
    TurnDirection.DIR_0: "0001",
    TurnDirection.DIR_1: "0010",
    TurnDirection.DIR_2: "0100",
    TurnDirection.DIR_3: "1000",
}

DENSE_TURN_INDICATORS: dict[TurnDirection, str] = {
    TurnDirection.DIR_0: "00",
    TurnDirection.DIR_1: "01",
    TurnDirection.DIR_2: "10",
    TurnDirection.DIR_3: "11",
}

XYZ_FILE_LINE_START_INDEX: int = 2  # First two lines are header in .xyz files

XYZ_FILE_PARTS_PER_LINE: int = 4  # Each line has symbol, x, y, z

RESULTS_DATA_DIRPATH: Path = ROOT_PROJECT_PATH / "output" / "results"

RAW_VQE_RESULTS_FILENAME: str = "raw_vqe_results.json"

XYZ_FILENAME: str = "conformation.xyz"

SPARSE_VQE_RESULTS_FILENAME: str = "sparse_vqe_results.json"

VQE_ITERATIONS_FILENAME: str = "vqe_iterations.txt"

GIF_VISUALIZATION_FILENAME: str = "rotating_3d_visualization.gif"

HTML_VISUALIZATION_FILENAME: str = "interactive_3d_visualization.html"

FLAT_VISUALIZATION_FILENAME: str = "conformation_2d.png"

TETRAHEDRAL_LATTICE_PADDING: int = 1

INDEX_COLNAME: str = "Index"

SYMBOL_COLNAME: str = "Symbol"

ITERATION_COLNAME: str = "Iteration"

COORDINATES_COLUMN_WIDTH: int = (
    12  # Width for coordinate columns in output files (sign, integer part, decimals)
)

FCC_EDGE_LENGTH: float = 1.0 / np.sqrt(3)

FCC_BASIS: NDArray[np.float64] = FCC_EDGE_LENGTH * np.array(
    [[-1, 1, 1], [1, 1, -1], [-1, -1, -1], [1, -1, 1]]
)

SIDE_CHAIN_FIFTH_POSITION_INDEX: int = (
    4  # Index of the 5th bead in zero-indexed beads list
)

BACKEND_TYPE: BackendType = BackendType.LOCAL_STATEVECTOR

IBM_QUANTUM_TOKEN: str | None = os.environ.get("IBM_QUANTUM_TOKEN")

IBM_QUANTUM_BACKEND_NAME: str | None = "ibm_marrakesh"

IBM_QUANTUM_SHOTS: int = 1024

MIN_CHAIN_LENGTH: int = 5  # Minimum length of the protein chain to be analyzed
