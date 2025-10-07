from qiskit.quantum_info import SparsePauliOp

from constants import CONFORMATION_ENCODING
from enums import ConformationEncoding
from exceptions import ConformationEncodingError
from protein.bead import Bead


class MainBead(Bead):
    """Represents a main bead in the protein's backbone."""

    def __init__(self, symbol: str, index: int, parent_chain_len: int) -> None:
        """
        Initialize the main bead and set up its turn qubits.

        Args:
            symbol (str): Amino acid symbol.
            index (int): Position of the bead in the chain.
            parent_chain_len (int): Total number of beads in the parent chain.

        """
        super().__init__(symbol=symbol, index=index, parent_chain_len=parent_chain_len)

    def turn_0(self) -> SparsePauliOp:
        """
        Return the Pauli operator for the turn in direction 0.

        Returns:
            SparsePauliOp: Pauli operator representing direction 0.

        Raises:
            ConformationEncodingError: If the conformation encoding is invalid.

        """
        if CONFORMATION_ENCODING == ConformationEncoding.SPARSE:
            return self.turn_qubits[0]
        if CONFORMATION_ENCODING == ConformationEncoding.DENSE:
            return self._dense_turn_fun_0()
        raise ConformationEncodingError

    def turn_1(self) -> SparsePauliOp:
        """
        Return the Pauli operator for the turn in direction 1.

        Returns:
            SparsePauliOp: Pauli operator representing direction 1.

        Raises:
            ConformationEncodingError: If the conformation encoding is invalid.

        """
        if CONFORMATION_ENCODING == ConformationEncoding.SPARSE:
            return self.turn_qubits[1]
        if CONFORMATION_ENCODING == ConformationEncoding.DENSE:
            return self._dense_turn_fun_1()
        raise ConformationEncodingError

    def turn_2(self) -> SparsePauliOp:
        """
        Return the Pauli operator for the turn in direction 2.

        Returns:
            SparsePauliOp: Pauli operator representing direction 2.

        Raises:
            ConformationEncodingError: If the conformation encoding is invalid.

        """
        if CONFORMATION_ENCODING == ConformationEncoding.SPARSE:
            return self.turn_qubits[2]
        if CONFORMATION_ENCODING == ConformationEncoding.DENSE:
            return self._dense_turn_fun_2()
        raise ConformationEncodingError

    def turn_3(self) -> SparsePauliOp:
        """
        Return the Pauli operator for the turn in direction 3.

        Returns:
            SparsePauliOp: Pauli operator representing direction 3.

        Raises:
            ConformationEncodingError: If the conformation encoding is invalid.

        """
        if CONFORMATION_ENCODING == ConformationEncoding.SPARSE:
            return self.turn_qubits[3]
        if CONFORMATION_ENCODING == ConformationEncoding.DENSE:
            return self._dense_turn_fun_3()
        raise ConformationEncodingError

    def _dense_turn_fun_0(self) -> SparsePauliOp:
        """
        Compute the dense encoding operator for direction 0.

        Returns:
            SparsePauliOp: Dense-encoded Pauli operator for direction 0.

        """
        return (
            (self._full_identity - self.turn_qubits[0])
            @ (self._full_identity - self.turn_qubits[1])
        ).simplify()

    def _dense_turn_fun_1(self) -> SparsePauliOp:
        """
        Compute the dense encoding operator for direction 1.

        Returns:
            SparsePauliOp: Dense-encoded Pauli operator for direction 1.

        """
        return (
            self.turn_qubits[1] @ (self.turn_qubits[1] - self.turn_qubits[0])
        ).simplify()

    def _dense_turn_fun_2(self) -> SparsePauliOp:
        """
        Compute the dense encoding operator for direction 2.

        Returns:
            SparsePauliOp: Dense-encoded Pauli operator for direction 2.

        """
        return (
            self.turn_qubits[0] @ (self.turn_qubits[0] - self.turn_qubits[1])
        ).simplify()

    def _dense_turn_fun_3(self) -> SparsePauliOp:
        """
        Compute the dense encoding operator for direction 3.

        Returns:
            SparsePauliOp: Dense-encoded Pauli operator for direction 3.

        """
        return (self.turn_qubits[0] @ self.turn_qubits[1]).simplify()
