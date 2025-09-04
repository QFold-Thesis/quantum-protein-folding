from collections import defaultdict
from qiskit.quantum_info import SparsePauliOp  # pyright: ignore[reportMissingTypeStubs]
from logger import get_logger
from utils.qubit_utils import build_full_identity
from protein import Protein

logger = get_logger()


class ContactMap:
    def __init__(self, protein: Protein):
        self.main_main_contacts: dict[int, dict[int, SparsePauliOp | None]] = defaultdict(
            lambda: defaultdict(lambda: None)
        )
        self.main_side_contacts: dict[int, dict[int, SparsePauliOp | None]] = defaultdict(
            lambda: defaultdict(lambda: None)
        )
        self.side_main_contacts: dict[int, dict[int, SparsePauliOp | None]] = defaultdict(
            lambda: defaultdict(lambda: None)
        )
        self.side_side_contacts: dict[int, dict[int, SparsePauliOp | None]] = defaultdict(
            lambda: defaultdict(lambda: None)
        )

        self._contacts_detected: int = 0


        self._num_contact_qubits: int = pow(len(protein.main_chain) - 1, 2)
        self._full_identity: SparsePauliOp = build_full_identity(num_qubits=self._num_contact_qubits)

        try:
            self.initialize_contact_map()
        except Exception as e:
            logger.error(f"Error in initializing contact map: {e}")
            raise e

        else:
            logger.debug("Contact map initialized successfully.")


    def initialize_contact_map(self):
        """Initializes all contact maps to empty dictionaries."""
        pass
    
