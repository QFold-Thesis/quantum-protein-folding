from pathlib import Path
from enums import ConformationEncoding

ROOT_PROJECT_PATH: Path = Path(__file__).parent.parent

MJ_INTERACTION_MATRIX_FILEPATH: Path = (
    ROOT_PROJECT_PATH / "src" / "resources" / "mj_matrix.txt"
)

CONFORMATION_ENCODING = ConformationEncoding.DENSE

EMPTY_SIDECHAIN_PLACEHOLDER = "_"
