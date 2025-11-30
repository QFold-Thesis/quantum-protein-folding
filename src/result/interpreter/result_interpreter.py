"""Utilities for interpreting protein folding results.

This module provides the `ResultInterpreter` class, which processes
the results of quantum simulations for protein folding, including detailed VQE results,
decoded turn sequences, and 3D coordinate mappings.
"""

from __future__ import annotations

import json
from typing import TYPE_CHECKING

import numpy as np

from constants import (
    CONFORMATION_ENCODING,
    COORDINATES_COLUMN_WIDTH,
    DENSE_TURN_INDICATORS,
    EMPTY_SIDECHAIN_PLACEHOLDER,
    FCC_BASIS,
    INDEX_COLNAME,
    ITERATION_COLNAME,
    QUBITS_PER_TURN,
    RAW_VQE_RESULTS_FILENAME,
    SIDE_CHAIN_FIFTH_POSITION_INDEX,
    SPARSE_TURN_INDICATORS,
    SPARSE_VQE_RESULTS_FILENAME,
    SYMBOL_COLNAME,
    VQE_ITERATIONS_FILENAME,
)
from enums import ConformationEncoding, TurnDirection
from exceptions import ConformationEncodingError
from logger import get_logger
from result.models import SparseVQEOutput
from utils.result_interpretation_utils import (
    BeadPosition,
    create_xyz_file,
    sanitize_for_json,
)

if TYPE_CHECKING:
    from pathlib import Path

    from numpy.typing import NDArray
    from qiskit_algorithms import SamplingMinimumEigensolverResult

    from protein import Protein

logger = get_logger()


class ResultInterpreter:
    """Interprets and processes the results of quantum protein folding simulations.

    Attributes:
        coordinates_3d (list[BeadPosition]): 3D coordinates of the protein beads.
        turn_sequence (list[TurnDirection]): Decoded sequence of turns for the protein chain.

    """

    def __init__(
        self,
        protein: Protein,
        dirpath: Path,
        raw_vqe_results: SamplingMinimumEigensolverResult,
        vqe_energies: list[float],
        vqe_iterations: list[int],
    ) -> None:
        """Initialize the ResultInterpreter with protein data and VQE results.

        Args:
            protein (Protein): The Protein object containing chain information.
            dirpath (Path): Directory path for saving output files.
            raw_vqe_results (SamplingMinimumEigensolverResult): Raw VQE results from quantum simulation.
            vqe_energies (list[float]): List of VQE energy values per iteration.
            vqe_iterations (list[int]): List of VQE iteration numbers.

        Raises:
            ConformationEncodingError: If the conformation encoding is not recognized.

        """
        self._dirpath: Path = dirpath
        self._raw_results: SamplingMinimumEigensolverResult = raw_vqe_results

        self._vqe_energies: list[float] = vqe_energies
        self._vqe_iterations: list[int] = vqe_iterations

        self._protein: Protein = protein
        self._fifth_bead_has_no_sidechain: bool = (
            len(self._protein.side_chain) >= SIDE_CHAIN_FIFTH_POSITION_INDEX + 1
            and self._protein.side_chain[SIDE_CHAIN_FIFTH_POSITION_INDEX].symbol
            == EMPTY_SIDECHAIN_PLACEHOLDER
        )

        if CONFORMATION_ENCODING == ConformationEncoding.DENSE:
            self._turn_encoding: dict[TurnDirection, str] = DENSE_TURN_INDICATORS
        elif CONFORMATION_ENCODING == ConformationEncoding.SPARSE:
            self._turn_encoding: dict[TurnDirection, str] = SPARSE_TURN_INDICATORS
        else:
            raise ConformationEncodingError

        self._vqe_output: SparseVQEOutput = self._interpret_raw_vqe_results()

        # Note - the sole reason for the bitstring here to be passed explicitly, is to ensure that we have a single
        # starting point for bitstring preprocessing. After preprocessing, all further methods use single, property.

        self._processed_bitstring: str = self._preprocess_bitstring(
            bitstring=self._vqe_output.bitstring
        )

        self._turn_sequence: list[TurnDirection] = self._generate_turn_sequence()
        self._log_turn_sequence()

        self._coordinates_3d: list[BeadPosition] = self._generate_3d_coordinates()
        self._log_coordinates_3d()

        self._main_main_contacts_detected: dict[int, int] = (
            self._find_main_main_contacts()
        )

    def _find_main_main_contacts(self) -> dict[int, int]:
        """Finds contacts between main chain beads based on the interaction bits.

        This reads the first part of the raw VQE bitstring, which contains
        the on/off flags for potential interactions.

        Returns:
            dict[int, int]: A dictionary mapping bead index `i` to bead index `j` for each detected interaction (we don't store contacts symmetrically, assume contacts are bidirectional).

        """
        logger.debug("Finding main-main contacts from interaction bits...")

        raw_bitstring: str = self._vqe_output.bitstring
        num_beads: int = len(self._protein.main_chain)

        num_shape_qubits: int = self._get_target_sequence_length_main_chain()
        interaction_bits: str = raw_bitstring[:-num_shape_qubits]

        contacts: dict[int, int] = {}
        qubit_index: int = 0

        for i in range(num_beads - 5):
            for j in range(i + 5, num_beads, 2):
                if qubit_index >= len(interaction_bits):
                    logger.warning(
                        "Ran out of interaction bits while checking pair (%s, %s). Expected more bits.",
                        i,
                        j,
                    )
                    break
                if interaction_bits[qubit_index] == "1":
                    contacts[i] = j

                qubit_index += 1
            if qubit_index >= len(interaction_bits):
                break

        if qubit_index < len(interaction_bits):
            logger.warning(
                "Finished checking all pairs, but %s interaction bits were left over.",
                len(interaction_bits) - qubit_index,
            )
        return contacts

    def _log_turn_sequence(self) -> None:
        """Logs the decoded turn sequence."""
        if not self._turn_sequence:
            logger.warning("No turn sequence to log.")
            return

        logger.info("Turn sequence decoded for %s turns.", len(self._turn_sequence))
        idx_width: int = len(str(len(self._turn_sequence)))
        val_width: int = max(len(str(t.value)) for t in self._turn_sequence)
        name_width: int = max(len(t.name) for t in self._turn_sequence)

        for i, turn in enumerate(self._turn_sequence, start=1):
            logger.info(
                "Turn %s - %s (%s)",
                f"{i:>{idx_width}}",
                f"{turn.value!s:>{val_width}}",
                f"{turn.name:<{name_width}}",
            )

    def _log_coordinates_3d(self) -> None:
        """Logs the generated 3D coordinates."""
        if not self._coordinates_3d:
            logger.warning("No 3D coordinates to log.")
            return

        logger.info("3D coordinates generated for %s beads.", len(self._coordinates_3d))

        idx_width: int = len(str(len(self._coordinates_3d) - 1))
        symbol_width: int = max(len(b.symbol) for b in self._coordinates_3d)
        coord_width: int = COORDINATES_COLUMN_WIDTH

        idx_col: int = max(len(INDEX_COLNAME), idx_width)
        sym_col: int = max(len(SYMBOL_COLNAME), symbol_width)
        c_col: int = coord_width

        header: str = (
            f"{INDEX_COLNAME:>{idx_col}}  "
            f"{SYMBOL_COLNAME:>{sym_col}}"
            f"{'X':>{c_col}}  "
            f"{'Y':>{c_col}}  "
            f"{'Z':>{c_col}}"
        )
        logger.info(header)

        for bead_pos in self._coordinates_3d:
            row = (
                f"{bead_pos.index:>{idx_col}}  "
                f"{bead_pos.symbol:>{sym_col}}"
                f"{bead_pos.x:>{c_col}.4f}  "
                f"{bead_pos.y:>{c_col}.4f}  "
                f"{bead_pos.z:>{c_col}.4f}"
            )
            logger.info(row)

    def _interpret_raw_vqe_results(self) -> SparseVQEOutput:
        """Interprets the raw VQE results into a structured SparseVQEOutput.

        Returns:
            SparseVQEOutput: The interpreted VQE results.

        Raises:
            ValueError: If the best measurement data is incomplete or missing.

        """
        logger.debug("Interpreting raw VQE results for %s...", self._dirpath)

        best_measurement = self._raw_results.best_measurement

        if not best_measurement:
            msg: str = "No best measurement found in VQE output."
            raise ValueError(msg)

        bitstring: str | None = best_measurement.get("bitstring")
        probability: float | None = best_measurement.get("probability")
        state: str | None = best_measurement.get("state")
        energy_value: np.complex128 | None = best_measurement.get("value")

        if (
            bitstring is None
            or probability is None
            or state is None
            or energy_value is None
        ):
            msg: str = "Incomplete best measurement data in VQE output."
            raise ValueError(msg)

        logger.info("VQE interpretation completed")
        logger.info("Best state found: %s (probability: %s)", state, probability)
        logger.info("Bitstring: %s", bitstring)
        logger.info("Energy value: %s", energy_value)

        return SparseVQEOutput(
            bitstring=bitstring,
            probability=probability,
            state=state,
            energy_value=energy_value,
        )

    def _preprocess_bitstring(self, bitstring: str) -> str:
        """Preprocesses the raw bitstring from VQE results to match the expected format.

        Note:
            Bitstring preprocessing includes:
            - Trimming to the target length.
            - Appending fixed turn indicators for the first two turns.
            - Handling special cases for the fifth bead's sidechain in dense encoding

            For details regarding processing, see _get_target_sequence_length_main_chain.

        Args:
            bitstring (str): The raw bitstring from VQE results.

        Returns:
            str: The preprocessed bitstring ready for turn sequence decoding.

        """
        logger.debug("Preprocessing bitstring for turn sequence decoding...")

        target_bitstring_length: int = self._get_target_sequence_length_main_chain()

        result_bitstring: str = bitstring[-target_bitstring_length:]

        result_bitstring: str = (
            result_bitstring
            + self._turn_encoding[TurnDirection.DIR_0][::-1]
            + self._turn_encoding[TurnDirection.DIR_1][::-1]
        )
        if CONFORMATION_ENCODING == ConformationEncoding.SPARSE:
            return result_bitstring[::-1]
        if CONFORMATION_ENCODING != ConformationEncoding.DENSE:
            raise ConformationEncodingError

        if self._fifth_bead_has_no_sidechain:
            result_bitstring = (
                result_bitstring[: -(SIDE_CHAIN_FIFTH_POSITION_INDEX + 1)]
                + "1"
                + result_bitstring[-(SIDE_CHAIN_FIFTH_POSITION_INDEX + 1) :]
            )
            logger.info(
                "Fifth bead has no sidechain. Turn #3 encoded as fixed '1' value."
            )

        logger.info(
            "Preprocessed bitstring to target length of %s bits: %s",
            target_bitstring_length,
            result_bitstring,
        )
        return result_bitstring[::-1]

    def _get_target_sequence_length_main_chain(self) -> int:
        """Calculates the target length (in bits) of the turn sequence for the main chain.

        Note:
            Each turn is represented by a fixed number of qubits (QUBITS_PER_TURN).
            The number of turns is N - 1, where N is the number of beads in the main chain.
            However, by the symmetry of tetrahedral lattice we can assume that:

            - The first two turns are fixed.
            - (Sparse encoding only) If we don't have a side bead attached to the second bead, we can encode the third turn as fixed value of "1".

        Returns:
            int: Target turn sequence length for the main chain in bits.

        Raises:
            ConformationEncodingError: If the conformation encoding is not recognized.

        """
        target_length: int = QUBITS_PER_TURN * (len(self._protein.main_chain) - 3)
        if CONFORMATION_ENCODING == ConformationEncoding.DENSE:
            target_length -= int(self._fifth_bead_has_no_sidechain)
        elif CONFORMATION_ENCODING != ConformationEncoding.SPARSE:
            raise ConformationEncodingError

        return target_length

    def _generate_turn_sequence(
        self,
    ) -> list[TurnDirection]:
        """Generates the turn sequence from the processed bitstring.

        Returns:
            list[TurnDirection]: The decoded sequence of turns.

        Raises:
            ConformationEncodingError: If the decoded turns contain unknown encodings.

        """
        logger.debug("Generating turn sequence from processed bitstring...")
        turns_length = len(self._processed_bitstring) // QUBITS_PER_TURN
        turns = [
            self._processed_bitstring[i * QUBITS_PER_TURN : (i + 1) * QUBITS_PER_TURN]
            for i in range(turns_length)
        ]

        bitstring_to_direction = {
            bitstring: direction for direction, bitstring in self._turn_encoding.items()
        }

        turn_sequence: list[TurnDirection] = []
        for turn in turns:
            if turn not in bitstring_to_direction:
                msg: str = f"Unknown turn encoding for: {turn}"
                raise ConformationEncodingError(msg)
            turn_sequence.append(bitstring_to_direction[turn])

        return turn_sequence

    def _generate_3d_coordinates(self) -> list[BeadPosition]:
        """Generates the coordinates for the beads in the main chain in 3D tetrahedral lattice.

        Note:
            The 3D coordinates are generated based on the tetrahedral lattice structure, which is used to represent the spatial arrangement of the beads in the protein.
            Calculations use the FCC basis vectors to determine the position of each bead based on the turn sequence and the sublattice they belong to (determined by the bead index - alternating signs).

        Returns:
            list[BeadPosition]: List of 3D coordinates for each bead in the main chain.

        Raises:
            ConformationEncodingError: If the decoded turns contain unknown encodings.

        """
        logger.debug("Generating 3D coordinates from processed bitstring...")
        tetrahedral_basis_vector: NDArray[np.float64] = FCC_BASIS.copy()

        tetrahedral_basis_vector /= np.linalg.norm(tetrahedral_basis_vector[0])

        bitstring_to_direction: dict[str, int] = {
            bitstring: direction.value
            for direction, bitstring in self._turn_encoding.items()
        }

        turns_length: int = len(self._processed_bitstring) // QUBITS_PER_TURN
        turns: list[str] = [
            self._processed_bitstring[i * QUBITS_PER_TURN : (i + 1) * QUBITS_PER_TURN]
            for i in range(turns_length)
        ]

        main_chain_symbols = [bead.symbol for bead in self._protein.main_chain.beads]

        current_pos = np.zeros(3)
        coords = [BeadPosition(0, main_chain_symbols[0], *current_pos)]

        for i, (turn, symbol) in enumerate(
            zip(turns, main_chain_symbols[1:], strict=True)
        ):
            if turn not in bitstring_to_direction:
                msg: str = f"Unknown turn encoding for: {turn}"
                raise ConformationEncodingError(msg)

            dir_idx = bitstring_to_direction[turn]
            direction_vec = ((-1) ** i) * tetrahedral_basis_vector[dir_idx]
            current_pos = current_pos + direction_vec
            coords.append(BeadPosition(len(coords), symbol, *current_pos))

        return coords

    def dump_results_to_files(self) -> None:
        """Dumps the interpreted results to output files in the specified directory."""
        create_xyz_file(self.coordinates_3d, self._dirpath)

        self._dump_result_dict_to_json(
            filename=RAW_VQE_RESULTS_FILENAME, results_dict=self._raw_results
        )
        self._dump_result_dict_to_json(
            filename=SPARSE_VQE_RESULTS_FILENAME, results_dict=self._vqe_output
        )

        self._dump_vqe_iterations_to_file(filename=VQE_ITERATIONS_FILENAME)

    def _dump_vqe_iterations_to_file(self, filename: str) -> None:
        """Dumps the VQE iterations and their corresponding energies to a text file.

        Args:
            filename (str): The name of the file to save the VQE iterations and energies.

        Raises:
            Exception: If there is an error writing to the file.

        """
        vqe_iterations_filepath: Path = self._dirpath / filename

        idx_width: int = len(ITERATION_COLNAME)
        try:
            with vqe_iterations_filepath.open("w", encoding="utf-8") as f:
                f.write(f"{ITERATION_COLNAME} | Energy\n")
                for iteration, energy in zip(
                    self._vqe_iterations, self._vqe_energies, strict=True
                ):
                    placeholder = " " if energy >= 0 else ""
                    f.write(f"{iteration:^{idx_width}} | {placeholder}{energy:.16f}\n")
        except Exception:
            logger.exception("Error dumping VQE iterations to file")
            raise
        else:
            logger.info("VQE iterations saved to %s", vqe_iterations_filepath)

    def _dump_result_dict_to_json(
        self,
        filename: str,
        results_dict: SparseVQEOutput | SamplingMinimumEigensolverResult,
    ) -> None:
        """Dumps the given results dictionary to a JSON file.

        Args:
            filename (str): The name of the file to save the results.
            results_dict (SparseVQEOutput | SamplingMinimumEigensolverResult): The results data to be saved.

        Raises:
            Exception: If there is an error during sanitization or file writing.

        """
        results_filepath: Path = self._dirpath / filename

        try:
            results_data = sanitize_for_json(results_dict)
            with results_filepath.open("w", encoding="utf-8") as f:
                json.dump(results_data, f, indent=2, ensure_ascii=False)
        except Exception:
            logger.exception("Error sanitizing results for JSON")
            raise
        else:
            logger.info("JSON results saved to %s", results_filepath)

    @property
    def coordinates_3d(self) -> list[BeadPosition]:
        """list[BeadPosition]: 3D coordinates of the protein beads."""
        return self._coordinates_3d

    @property
    def turn_sequence(self) -> list[TurnDirection]:
        """list[TurnDirection]: Sequence of turns in the protein folding."""
        return self._turn_sequence

    @property
    def main_main_contacts_detected(self) -> dict[int, int]:
        """dict[int, int]: Main-main contacts detected in the result sequence."""
        return self._main_main_contacts_detected
