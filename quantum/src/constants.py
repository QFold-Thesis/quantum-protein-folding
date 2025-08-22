from pathlib import Path


ROOT_PROJECT_PATH: Path = Path(__file__).parent.parent

MJ_INTERACTION_MATRIX_FILEPATH: Path = (
    ROOT_PROJECT_PATH / "src" / "resources" / "mj_matrix.txt"
)
