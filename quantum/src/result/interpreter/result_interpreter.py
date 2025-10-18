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
    QUBITS_PER_TURN,
    RAW_VQE_RESULTS_FILENAME,
    SIDE_CHAIN_FIFTH_POSITION_INDEX,
    SPARSE_TURN_INDICATORS,
    SPARSE_VQE_RESULTS_FILENAME,
    SYMBOL_COLNAME,
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
    def __init__(
        self,
        protein: Protein,
        dirpath: Path,
        raw_vqe_results: SamplingMinimumEigensolverResult,
    ) -> None:
        self._dirpath: Path = dirpath
        self._raw_results: SamplingMinimumEigensolverResult = raw_vqe_results

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
        # starting point for bitstring preprocessing. After preprocessing, all further methods use single, property-bitstring.

        self._processed_bitstring: str = self._preprocess_bitstring(
            bitstring=self._vqe_output.bitstring
        )

        self._turn_sequence: list[TurnDirection] = self._generate_turn_sequence()
        self._log_turn_sequence()

        self._coordinates_3d: list[BeadPosition] = self._generate_3d_coordinates()
        self._log_coordinates_3d()

    def _log_turn_sequence(self) -> None:
        logger.info(f"Turn sequence decoded for {len(self._turn_sequence)} turns.")
        idx_width: int = len(str(len(self._turn_sequence)))
        val_width: int = max(len(str(t.value)) for t in self._turn_sequence)
        name_width: int = max(len(t.name) for t in self._turn_sequence)

        for i, turn in enumerate(self._turn_sequence, start=1):
            logger.info(
                f"Turn {i:>{idx_width}} - {turn.value!s:>{val_width}} ({turn.name:<{name_width}})"
            )

    def _log_coordinates_3d(self) -> None:
        logger.info(f"3D coordinates generated for {len(self._coordinates_3d)} beads.")

        if not self._coordinates_3d:
            return

        idx_width: int = len(str(len(self._coordinates_3d) - 1))
        symbol_width: int = max(len(b.symbol) for b in self._coordinates_3d)
        coord_width: int = COORDINATES_COLUMN_WIDTH

        idx_col: int = max(len(INDEX_COLNAME), idx_width)
        sym_col: int = max(len(SYMBOL_COLNAME), symbol_width)
        c_col: int = coord_width

        header = (
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
        logger.info(f"Interpreting raw VQE results for {self._dirpath}")

        best_measurement = self._raw_results.best_measurement

        if not best_measurement:
            msg: str = "No best measurement found in VQE output."
            raise ValueError(msg)

        bitstring: str | None = best_measurement.get("bitstring")
        probability: float | None = best_measurement.get("probability")
        state: str | None = best_measurement.get("state")
        energy_value: np.complex128 | None = best_measurement.get("value")

        if None in (bitstring, probability, state, energy_value):
            msg: str = "Incomplete best measurement data in VQE output."
            raise ValueError(msg)

        logger.info("VQE interpretation complete")
        logger.info(f"Best state found: {state} (probability: {probability})")
        logger.info(f"Bitstring: {bitstring}")
        logger.info(f"Energy value: {energy_value}")

        return SparseVQEOutput(
            bitstring=bitstring,
            probability=probability,
            state=state,
            energy_value=energy_value,
        )

    def _preprocess_bitstring(self, bitstring: str) -> str:
        """Preprocesses the bitstring by appending initial turns and reversing it."""
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
                "Fifth bead has no sidechain. Turn 3 encoded as fixed '1' value."
            )

        logger.info(
            f"Preprocessed bitstring to target length of {target_bitstring_length} bits: {result_bitstring}"
        )
        return result_bitstring[::-1]

    def _get_target_sequence_length_main_chain(self) -> int:
        """
        Calculates the target length (in bits) of the turn sequence for the main chain.

        Notes:
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
                logger.warning(f"Unknown turn encoding: {turn}")
                continue

            dir_idx = bitstring_to_direction[turn]
            direction_vec = ((-1) ** i) * tetrahedral_basis_vector[dir_idx]
            current_pos = current_pos + direction_vec
            coords.append(BeadPosition(len(coords), symbol, *current_pos))

        return coords

    def dump_results_to_files(self) -> None:
        create_xyz_file(self.coordinates_3d, self._dirpath)

        self._dump_result_dict_to_json(
            filename=RAW_VQE_RESULTS_FILENAME, results_dict=self._raw_results
        )
        self._dump_result_dict_to_json(
            filename=SPARSE_VQE_RESULTS_FILENAME, results_dict=self._vqe_output
        )

    def _dump_result_dict_to_json(
        self,
        filename: str,
        results_dict: SparseVQEOutput | SamplingMinimumEigensolverResult,
    ) -> None:
        results_filepath: Path = self._dirpath / filename

        try:
            results_data = sanitize_for_json(results_dict)
            with results_filepath.open("w", encoding="utf-8") as f:
                json.dump(results_data, f, indent=2, ensure_ascii=False)
        except Exception:
            logger.exception("Error sanitizing results for JSON")
            raise
        else:
            logger.info(f"JSON results saved to {results_filepath}")

    @property
    def coordinates_3d(self) -> list[BeadPosition]:
        return self._coordinates_3d

    @property
    def turn_sequence(self) -> list[TurnDirection]:
        return self._turn_sequence
