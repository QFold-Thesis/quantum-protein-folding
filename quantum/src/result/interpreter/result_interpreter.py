from __future__ import annotations

import json
from typing import TYPE_CHECKING

import numpy as np

from constants import (
    CONFORMATION_ENCODING,
    DENSE_TURN_INDICATORS,
    QUBITS_PER_TURN,
    RAW_VQE_RESULTS_FILENAME,
    SPARSE_TURN_INDICATORS,
    SPARSE_VQE_RESULTS_FILENAME,
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
        self.dirpath: Path = dirpath
        self.raw_results: SamplingMinimumEigensolverResult = raw_vqe_results

        self.main_chain_symbols: list[str] = [
            bead.symbol for bead in protein.main_chain.beads
        ]
        self.side_chain_symbols: list[str] = [
            bead.symbol for bead in protein.side_chain.beads
        ]

        if CONFORMATION_ENCODING == ConformationEncoding.DENSE:
            self.turn_encoding: dict[TurnDirection, str] = DENSE_TURN_INDICATORS
        elif CONFORMATION_ENCODING == ConformationEncoding.SPARSE:
            self.turn_encoding: dict[TurnDirection, str] = SPARSE_TURN_INDICATORS
        else:
            raise ConformationEncodingError

        self.vqe_output: SparseVQEOutput = self._interpret_raw_vqe_results()
        self.formatted_bitstring: str = self._preprocess_bitstring(
            self.vqe_output.bitstring
        )

        self.coordinates_3d: list[BeadPosition] = self._generate_3d_coordinates(
            bitstring=self.formatted_bitstring,
        )

    def save_to_files(self) -> None:
        create_xyz_file(self.coordinates_3d, self.dirpath)

        self._dump_result_dict_to_json(
            filename=RAW_VQE_RESULTS_FILENAME, results_dict=self.raw_results
        )
        self._dump_result_dict_to_json(
            filename=SPARSE_VQE_RESULTS_FILENAME, results_dict=self.vqe_output
        )

    def _interpret_raw_vqe_results(self) -> SparseVQEOutput:
        logger.info(f"Interpreting raw VQE results for {self.dirpath}")

        best_measurement = self.raw_results.best_measurement

        if not best_measurement:
            msg = "No best measurement found in VQE output."
            raise ValueError(msg)

        bitstring: str | None = best_measurement.get("bitstring")
        probability: float | None = best_measurement.get("probability")
        state: str | None = best_measurement.get("state")
        energy_value: np.complex128 | None = best_measurement.get("value")

        if None in (bitstring, probability, state, energy_value):
            msg = "Incomplete best measurement data in VQE output."
            raise ValueError(msg)

        return SparseVQEOutput(
            bitstring=bitstring,
            probability=probability,
            state=state,
            energy_value=energy_value,
        )

    def _generate_3d_coordinates(
        self,
        bitstring: str,
    ) -> list[BeadPosition]:
        tetra_dirs = np.array(
            [[1, 1, 1], [1, -1, -1], [-1, 1, -1], [-1, -1, 1]], dtype=float
        )
        tetra_dirs /= np.linalg.norm(tetra_dirs[0])

        bitstring_to_direction = {
            bitstring: direction.value
            for direction, bitstring in self.turn_encoding.items()
        }

        turns_length = len(bitstring) // QUBITS_PER_TURN
        turns = [
            bitstring[i * QUBITS_PER_TURN : (i + 1) * QUBITS_PER_TURN]
            for i in range(turns_length)
        ]

        # initialize the starting position (first bead)
        current_pos = np.array([0.0, 0.0, 0.0])
        coords = [
            BeadPosition(
                index=0,
                symbol=self.main_chain_symbols[0],
                x=current_pos[0],
                y=current_pos[1],
                z=current_pos[2],
            )
        ]

        for turn, symbol in zip(turns, self.main_chain_symbols[1::], strict=True):
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

    def _preprocess_bitstring(self, bitstring: str) -> str:
        """Preprocesses the bitstring by appending initial turns and reversing it."""
        return "".join(
            reversed(
                bitstring
                + self.turn_encoding[TurnDirection.DIR_1]
                + self.turn_encoding[TurnDirection.DIR_2]
            )
        )

    def _dump_result_dict_to_json(
        self,
        filename: str,
        results_dict: SparseVQEOutput | SamplingMinimumEigensolverResult,
    ) -> None:
        results_filepath: Path = self.dirpath / filename

        try:
            results_data = sanitize_for_json(results_dict)
            with results_filepath.open("w", encoding="utf-8") as f:
                json.dump(results_data, f, indent=2, ensure_ascii=False)
        except Exception:
            logger.exception("Error sanitizing results for JSON")
            raise
