from pathlib import Path

from enums import ConformationEncoding

ROOT_PROJECT_PATH: Path = Path(__file__).parent.parent

MJ_INTERACTION_MATRIX_FILEPATH: Path = (
    ROOT_PROJECT_PATH / "src" / "resources" / "mj_matrix.txt"
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

GLOBAL_LOGGER_NAME: str = "global_logger"
