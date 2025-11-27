"""Definitions of constants used throughout the protein folding simulations."""

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

OUTPUT_DATA_DIR: Path = ROOT_PROJECT_PATH / "output" / "results"

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
