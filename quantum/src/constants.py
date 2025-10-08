import logging
from pathlib import Path

from enums import ConformationEncoding, InteractionType, TurnDirection

ROOT_PROJECT_PATH: Path = Path(__file__).parent.parent

MJ_INTERACTION_MATRIX_FILEPATH: Path = (
    ROOT_PROJECT_PATH / "src" / "resources" / "mj_matrix.txt"
)

HP_INTERACTION_MATRIX_FILEPATH: Path = (
    ROOT_PROJECT_PATH / "src" / "resources" / "hp_matrix.txt"
)

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
    TurnDirection.DIR_1: "0001",
    TurnDirection.DIR_2: "0010",
    TurnDirection.DIR_3: "0100",
    TurnDirection.DIR_4: "1000",
}

DENSE_TURN_INDICATORS: dict[TurnDirection, str] = {
    TurnDirection.DIR_1: "00",
    TurnDirection.DIR_2: "01",
    TurnDirection.DIR_3: "10",
    TurnDirection.DIR_4: "11",
}

XYZ_FILE_LINE_START_INDEX: int = 2  # First two lines are header in .xyz files

XYZ_FILE_PARTS_PER_LINE: int = 4  # Each line has symbol, x, y, z

OUTPUT_DATA_DIR: Path = ROOT_PROJECT_PATH / "output"
