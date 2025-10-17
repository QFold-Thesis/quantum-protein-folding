from __future__ import annotations

import json
from typing import TYPE_CHECKING

import numpy as np

from constants import (
    CONFORMATION_ENCODING,
    COORDINATES_COLUMN_WIDTH,
    DENSE_TURN_INDICATORS,
    INDEX_COLNAME,
    QUBITS_PER_TURN,
    RAW_VQE_RESULTS_FILENAME,
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

        self._main_chain_symbols: list[str] = [
            bead.symbol for bead in protein.main_chain.beads
        ]
        self._side_chain_symbols: list[str] = [
            bead.symbol for bead in protein.side_chain.beads
        ]

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
        return "".join(
            reversed(
                bitstring
                + self._turn_encoding[TurnDirection.DIR_1]
                + self._turn_encoding[TurnDirection.DIR_2]
            )
        )

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

    def _generate_3d_coordinates(
        self,
    ) -> list[BeadPosition]:
        tetra_dirs = np.array(
            [[1, 1, 1], [1, -1, -1], [-1, 1, -1], [-1, -1, 1]], dtype=float
        )
        tetra_dirs /= np.linalg.norm(tetra_dirs[0])

        bitstring_to_direction = {
            bitstring: direction.value
            for direction, bitstring in self._turn_encoding.items()
        }

        turns_length = len(self._processed_bitstring) // QUBITS_PER_TURN
        turns = [
            self._processed_bitstring[i * QUBITS_PER_TURN : (i + 1) * QUBITS_PER_TURN]
            for i in range(turns_length)
        ]

        # initialize the starting position (first bead)
        current_pos = np.array([0.0, 0.0, 0.0])
        coords = [
            BeadPosition(
                index=0,
                symbol=self._main_chain_symbols[0],
                x=current_pos[0],
                y=current_pos[1],
                z=current_pos[2],
            )
        ]

        for turn, symbol in zip(turns, self._main_chain_symbols[1::], strict=True):
            if turn not in bitstring_to_direction:
                logger.warning(f"Unknown turn encoding: {turn}")
                continue
            direction_idx = bitstring_to_direction[turn]
            direction = tetra_dirs[direction_idx]
            current_pos = current_pos + direction
            coords.append(
                BeadPosition(
                    index=len(coords),
                    symbol=symbol,
                    x=current_pos[0],
                    y=current_pos[1],
                    z=current_pos[2],
                )
            )

        return coords

    def save_to_files(self) -> None:
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
