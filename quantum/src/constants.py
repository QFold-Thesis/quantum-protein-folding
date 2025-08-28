from pathlib import Path

from enums import ConformationEncoding

from exceptions import ConformationEncodingError

ROOT_PROJECT_PATH: Path = Path(__file__).parent.parent

MJ_INTERACTION_MATRIX_FILEPATH: Path = (
    ROOT_PROJECT_PATH / "src" / "resources" / "mj_matrix.txt"
)

CONFORMATION_ENCODING: ConformationEncoding = ConformationEncoding.DENSE

try:
    QUBITS_PER_TURN: int = CONFORMATION_ENCODING.value
except ValueError:
    raise ConformationEncodingError

EMPTY_SIDECHAIN_PLACEHOLDER = "_"
