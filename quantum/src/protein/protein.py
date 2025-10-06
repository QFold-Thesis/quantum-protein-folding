from __future__ import annotations

from datetime import datetime, timezone
from pathlib import Path

from exceptions import ChainLengthError
from logger import get_logger
from protein.chain import MainChain

logger = get_logger()


class Protein:
    def __init__(self, main_protein_sequence: str, side_protein_sequences: list[str]) -> None:
        if len(main_protein_sequence) != len(side_protein_sequences):
            msg = "Main and side protein sequences must be of the same length."
            logger.error(msg)
            raise ChainLengthError(msg)

        self.main_chain: MainChain = MainChain(main_protein_sequence, side_protein_sequences)

        # Automatically pretty-print after construction (debug/inspection aid)
        try:
            self.pretty_print_to_file()
        except Exception:  # pragma: no cover - defensive
            logger.exception("Failed to pretty print protein structure.")

    # ------------------------------------------------------------------
    # Pretty printing utilities
    # ------------------------------------------------------------------
    def pretty_print(self) -> str:
        """
        Build a human-readable multi-line representation of the protein.

        Format:
            MAIN : <concatenated main chain symbols>
            BEAD i (symbol): side_chain=<side symbols or ->

        Empty / missing side chain yields '-'.
        """
        lines: list[str] = []
        main_symbols = "".join(bead.symbol for bead in self.main_chain)
        lines.append(f"MAIN : {main_symbols}")
        for idx, main_bead in enumerate(self.main_chain):
            side_chain = getattr(main_bead, "side_chain", None)
            if side_chain and getattr(side_chain, "beads", None):
                side_symbols = "".join(b.symbol for b in side_chain.beads)
            else:
                side_symbols = "-"
            lines.append(f"BEAD {idx:02d} ({main_bead.symbol}): side_chain={side_symbols}")
        return "\n".join(lines)

    def pretty_print_to_file(self, filename: str | None = None) -> Path:
        """
        Write the pretty print output to a timestamped (or given) file.

        Returns the path to the written file.
        """
        if filename is None:
            timestamp = datetime.now(tz=timezone.utc).strftime("%Y%m%d_%H%M%S")
            filename = f"protein_structure_{timestamp}.txt"
        path = Path(filename)
        content = self.pretty_print() + "\n"
        path.write_text(content, encoding="utf-8")
        logger.info("Protein structure written to %s", path.resolve())
        return path
