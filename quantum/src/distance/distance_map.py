from collections import defaultdict
from pathlib import Path

from qiskit.quantum_info import SparsePauliOp  # pyright: ignore[reportMissingTypeStubs]

from logger import get_logger
from protein import Protein
from protein.bead import Bead
from utils.qubit_utils import fix_qubits

logger = get_logger()

class DistanceMap:

    def __init__(self, protein: Protein):
        self.protein = protein
        self._distance_map_axes: list[defaultdict] = [
            defaultdict(lambda: defaultdict(int)),
            defaultdict(lambda: defaultdict(int)),
            defaultdict(lambda: defaultdict(int)),
            defaultdict(lambda: defaultdict(int)),
        ]
        self._calc_distances_main_chain()
        print(self._pretty_print_distance_map_axes())
        logger.debug(self._distance_map_axes)
        # also save to file for inspection
        try:
            self.save_distance_map_axes()
        except Exception:
            logger.exception("Failed to save distance_map_axes to file")

    def _format_operator_brief(self, op: object) -> str:
        """Return a short, human-readable representation for an operator or value."""
        # Simple cases first
        try:
            if op == 0:
                return "0"
        except Exception:
            pass

        # Prefer short repr but truncate long strings
        try:
            s = repr(op)
        except Exception:
            try:
                s = str(type(op).__name__)
            except Exception:
                s = "<unprintable>"

        # Normalize whitespace and truncate
        s = " ".join(s.split())
        max_len = 120
        if len(s) > max_len:
            s = s[: max_len - 3] + "..."
        return s

    def _pretty_print_distance_map_axes(self) -> str:
        """Format `self._distance_map_axes` into a readable multi-line string."""
        lines = ["_distance_map_axes:"]
        for axis_idx, dist_map_ax in enumerate(self._distance_map_axes):
            lines.append(f"  Axis {axis_idx}:")
            if not dist_map_ax:
                lines.append("    (no entries)")
                continue
            for lower, inner in dist_map_ax.items():
                for upper, op in inner.items():
                    # Render `Bead` keys as `main_chain_{index}` when they belong
                    # to this DistanceMap's protein main_chain for readability.
                    def bead_label(b):
                        try:
                            for i, mb in enumerate(self.protein.main_chain):
                                if mb is b:
                                    return f"main_chain_{i+1}"
                        except Exception:
                            pass
                        return str(b)

                    lower_s = bead_label(lower)
                    upper_s = bead_label(upper)
                    op_s = self._format_operator_brief(op)
                    lines.append(f"    {lower_s} -> {upper_s}: {op_s}")
        return "\n".join(lines)

    def _calc_distances_main_chain(self):
        main_chain_len = len(self.protein.main_chain)

        for lower_bead_idx in range(1, main_chain_len):
            for upper_bead_idx in range(
                lower_bead_idx + 1, main_chain_len + 1
            ):
                lower_bead: Bead = self.protein.main_chain[lower_bead_idx - 1]
                upper_bead: Bead = self.protein.main_chain[upper_bead_idx - 1]

                for k in range(lower_bead_idx, upper_bead_idx):
                    indic_funcs = self.protein.main_chain[k - 1].turn_funcs()
                    for indic_fun_x, dist_map_ax in zip(
                        indic_funcs, self._distance_map_axes
                    ):
                        dist_map_ax[lower_bead][upper_bead] += (
                            -1
                        ) ** k * indic_fun_x

                for dist_map_ax in self._distance_map_axes:
                    dist_map_ax[lower_bead][upper_bead] = fix_qubits(
                        dist_map_ax[lower_bead][upper_bead]
                    )

    def save_distance_map_axes(self, path: str | Path = "distance_map_axes.txt") -> Path:
        """Save pretty-printed `_distance_map_axes` to a text file and return the Path.

        The file will contain the exact output of `_pretty_print_distance_map_axes()`.
        """
        p = Path(path)
        p.write_text(self._pretty_print_distance_map_axes(), encoding="utf-8")
        return p
