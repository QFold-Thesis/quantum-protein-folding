"""Analysis of external-field influence on the quantum protein folding energy landscape.

This module provides the :class:`FieldInfluenceAnalysis` class, which orchestrates
comparative VQE runs with different external-field configurations and produces
publication-ready plots of:

* Minimum energy vs. field strength λ (uniform field sweep)
* Per-bitstring probability distributions for each field variant
* Energy comparison bar chart across all field scenarios

Usage example::

    from src.analysis.field_influence_analysis import FieldInfluenceAnalysis
    from src.enums import InteractionType

    analysis = FieldInfluenceAnalysis(
        main_chain="HPPHH",
        interaction_type=InteractionType.HP,
        uniform_lambdas=[0.1, 0.5, 1.0, 2.0],
        vqe_max_iter=100,
    )
    analysis.run()
    analysis.plot(output_dir=Path("output/analysis"))
"""

from __future__ import annotations

import dataclasses
from pathlib import Path
from typing import TYPE_CHECKING, Any

import numpy as np
from qiskit.circuit.library import real_amplitudes
from qiskit.quantum_info import SparsePauliOp
from qiskit_algorithms import SamplingVQE
from qiskit_algorithms.optimizers import COBYLA

from backend import get_sampler
from builder import HamiltonianBuilder
from constants import EMPTY_SIDECHAIN_PLACEHOLDER
from contact import ContactMap
from distance import DistanceMap
from enums import InteractionType
from exceptions import InvalidInteractionTypeError
from interaction import HPInteraction, MJInteraction
from logger import get_logger
from particle.external_field import ExternalField
from protein import Protein
from utils.qubit_utils import remove_unused_qubits

if TYPE_CHECKING:
    from qiskit_algorithms import SamplingMinimumEigensolverResult

logger = get_logger()


# ---------------------------------------------------------------------------
# Data containers
# ---------------------------------------------------------------------------


@dataclasses.dataclass
class ScenarioResult:
    """Stores the result of a single VQE run for one field scenario.

    Attributes:
        label (str): Human-readable scenario label (e.g. ``"λ=0.5 (uniform)"``).
        field (ExternalField | None): The external field used, or ``None`` for
            the baseline (field-free) run.
        minimum_energy (float): Lowest energy found by the VQE.
        best_bitstring (str): Bitstring of the optimal quantum state.
        state_probabilities (dict[str, float]): Mapping from bitstring to its
            probability in the best measurement distribution.
        vqe_iterations (list[int]): Evaluation count at each VQE callback step.
        vqe_energies (list[float]): Energy value at each VQE callback step.

    """

    label: str
    field: ExternalField | None
    minimum_energy: float
    best_bitstring: str
    state_probabilities: dict[str, float]
    vqe_iterations: list[int]
    vqe_energies: list[float]


# ---------------------------------------------------------------------------
# Main analysis class
# ---------------------------------------------------------------------------


class FieldInfluenceAnalysis:
    """Comparative VQE analysis of external-field influence on protein folding.

    Runs VQE for each of the following field configurations and collects results:

    1. **Baseline** - no external field (``external_field=None``).
    2. **Uniform sweep** - uniform field at each lambda in *uniform_lambdas*.
    3. **Non-uniform** - stronger field at the central bead(s) of the chain,
       falling off toward the termini.

    After calling :meth:`run`, the collected results can be visualised with
    :meth:`plot`.

    Attributes:
        main_chain (str): Amino-acid sequence of the protein main chain.
        side_chain (str): Side-chain sequence (all ``_`` by default).
        interaction_type (InteractionType): HP or MJ interaction model.
        uniform_lambdas (list[float]): Uniform-field strengths to sweep over.
        vqe_max_iter (int): Maximum COBYLA iterations per VQE run.
        results (list[ScenarioResult]): Populated after :meth:`run` is called.

    """

    def __init__(
        self,
        main_chain: str,
        *,
        interaction_type: InteractionType = InteractionType.HP,
        uniform_lambdas: list[float] | None = None,
        side_chain: str | None = None,
        vqe_max_iter: int = 100,
        vqe_aggregation: float = 0.1,
    ) -> None:
        """Initialise the analysis with protein and field configuration.

        Args:
            main_chain (str): One-letter amino-acid sequence of the main chain.
                Must be compatible with the chosen *interaction_type*.
            interaction_type (InteractionType, optional): Interaction model to
                use.  Defaults to HP.
            uniform_lambdas (list[float] | None, optional): Field strengths for
                the uniform sweep.  Defaults to ``[0.1, 0.5, 1.0, 2.0]``.
            side_chain (str | None, optional): Side-chain sequence.  Defaults to
                all ``_`` (no side chains).
            vqe_max_iter (int, optional): Maximum COBYLA iterations per run.
                Defaults to 100.
            vqe_aggregation (float, optional): CVaR aggregation fraction for
                SamplingVQE.  Defaults to 0.1.

        Raises:
            InvalidInteractionTypeError: If *interaction_type* is unrecognised.

        """
        self.main_chain: str = main_chain
        self.side_chain: str = (
            side_chain
            if side_chain is not None
            else EMPTY_SIDECHAIN_PLACEHOLDER * len(main_chain)
        )
        self.interaction_type: InteractionType = interaction_type
        self.uniform_lambdas: list[float] = (
            uniform_lambdas if uniform_lambdas is not None else [0.1, 0.5, 1.0, 2.0]
        )
        self.vqe_max_iter: int = vqe_max_iter
        self.vqe_aggregation: float = vqe_aggregation

        self.results: list[ScenarioResult] = []

        # Build shared protein components (reused across all runs).
        self._interaction = self._build_interaction()
        self._protein = Protein(
            main_protein_sequence=self.main_chain,
            side_protein_sequence=self.side_chain,
            valid_symbols=self._interaction.valid_symbols,
        )
        self._contact_map = ContactMap(protein=self._protein)
        self._distance_map = DistanceMap(protein=self._protein)

        chain_len = len(self._protein.main_chain)
        logger.info(
            "FieldInfluenceAnalysis initialised: chain=%s (len=%d), model=%s",
            self.main_chain,
            chain_len,
            self.interaction_type.name,
        )

    # ------------------------------------------------------------------
    # Public interface
    # ------------------------------------------------------------------

    def run(self) -> None:
        """Execute all field scenarios and populate :attr:`results`.

        Runs VQE for:
        - baseline (no field)
        - one uniform run per λ in :attr:`uniform_lambdas`
        - one non-uniform run (Gaussian-like profile centred on the chain)

        All results are appended to :attr:`results` in order.

        """
        logger.info("Starting FieldInfluenceAnalysis.run() ...")
        self.results.clear()

        # 1. Baseline (no field)
        logger.info("--- Scenario: baseline (no field) ---")
        self.results.append(self._run_scenario(label="baseline (no field)", field=None))

        # 2. Uniform sweep
        for lam in self.uniform_lambdas:
            label = f"lambda={lam} (uniform)"
            logger.info("--- Scenario: %s ---", label)
            field = ExternalField.uniform(strength=lam)
            self.results.append(self._run_scenario(label=label, field=field))

        # 3. Non-uniform: Gaussian-like profile centred on the chain midpoint
        nu_field = self._build_non_uniform_field()
        label = "non-uniform (centre boost)"
        logger.info("--- Scenario: %s ---", label)
        self.results.append(self._run_scenario(label=label, field=nu_field))

        logger.info(
            "FieldInfluenceAnalysis.run() complete: %d scenarios finished.",
            len(self.results),
        )

    def plot(self, output_dir: Path | None = None) -> None:
        """Generate and save (or display) all comparison plots.

        Produces three figures:

        1. **Energy vs lambda** - minimum VQE energy on the y-axis, uniform-field lambda
           on the x-axis.  The baseline and non-uniform scenarios are drawn as
           horizontal reference lines.
        2. **Energy comparison bar chart** - one bar per scenario.
        3. **Probability distributions** - stacked subplots, one per scenario,
           showing the top-k bitstring probabilities.

        Args:
            output_dir (Path | None, optional): Directory to save the figures.
                If ``None`` the figures are displayed interactively.  The
                directory is created if it does not exist.

        Raises:
            RuntimeError: If :meth:`run` has not been called yet.

        """
        if not self.results:
            msg = "No results to plot. Call run() first."
            raise RuntimeError(msg)

        try:
            import matplotlib.pyplot as plt  # noqa: PLC0415
            from matplotlib import ticker  # noqa: PLC0415
        except ImportError as e:
            msg = "matplotlib is required for plotting. Install it with: pip install matplotlib"
            raise ImportError(msg) from e

        if output_dir is not None:
            output_dir = Path(output_dir)
            output_dir.mkdir(parents=True, exist_ok=True)

        # Consistent colour palette
        _palette = [
            "#1f77b4",
            "#ff7f0e",
            "#2ca02c",
            "#d62728",
            "#9467bd",
            "#8c564b",
            "#e377c2",
            "#7f7f7f",
        ]

        # ---- Figure 1: energy vs lambda ----------------------------------
        self._plot_energy_vs_lambda(plt, ticker, _palette, output_dir)

        # ---- Figure 2: bar chart of all scenarios ------------------------
        self._plot_energy_bar(plt, _palette, output_dir)

        # ---- Figure 3: probability distributions -------------------------
        self._plot_probability_distributions(plt, _palette, output_dir)

        if output_dir is None:
            plt.show()
        else:
            logger.info("All plots saved to %s", output_dir)

    # ------------------------------------------------------------------
    # Internal helpers - scenario execution
    # ------------------------------------------------------------------

    def _build_interaction(self) -> HPInteraction | MJInteraction:
        """Build the interaction model matching :attr:`interaction_type`."""
        if self.interaction_type == InteractionType.HP:
            return HPInteraction()
        if self.interaction_type == InteractionType.MJ:
            return MJInteraction()
        msg = f"Unknown interaction type: {self.interaction_type}"
        raise InvalidInteractionTypeError(msg)

    def _build_non_uniform_field(self) -> ExternalField:
        """Build a non-uniform field with a Gaussian-like profile centred on the chain.

        The field energy at bead ``i`` is:

            E_field((i,)) = lambda_max * exp(-((i - mid) / sigma)^2)

        where ``mid`` is the fractional midpoint of the chain and sigma = N/4.
        This produces a smooth peak at the centre that falls to ~2 % of the
        maximum at the termini for a typical chain length.

        Returns:
            ExternalField: Non-uniform field instance.

        """
        chain_len: int = len(self._protein.main_chain)
        mid: float = (chain_len - 1) / 2.0
        sigma: float = max(chain_len / 4.0, 1.0)
        lambda_max: float = -1.0  # attractive (negative = stabilising)

        energy_map: dict[tuple[int, ...], float] = {
            (i,): lambda_max * float(np.exp(-((i - mid) ** 2) / sigma**2))
            for i in range(chain_len)
        }
        logger.debug(
            "Non-uniform field profile (Gaussian, lambda_max=%s, sigma=%.2f): %s",
            lambda_max,
            sigma,
            {k: round(v, 4) for k, v in energy_map.items()},
        )
        return ExternalField.non_uniform(energy_map, default_energy=0.0)

    def _run_scenario(
        self,
        label: str,
        field: ExternalField | None,
    ) -> ScenarioResult:
        """Run a single VQE scenario and return its :class:`ScenarioResult`.

        Args:
            label (str): Human-readable label for this scenario.
            field (ExternalField | None): External field to use, or ``None``
                for the baseline run.

        Returns:
            ScenarioResult: Populated result dataclass.

        """
        # Build Hamiltonian
        h_builder = HamiltonianBuilder(
            protein=self._protein,
            interaction=self._interaction,
            distance_map=self._distance_map,
            contact_map=self._contact_map,
            external_field=field,
        )
        full_h: SparsePauliOp = h_builder.sum_hamiltonians()
        compressed_h: SparsePauliOp = remove_unused_qubits(full_h)

        logger.info(
            "Scenario '%s': Hamiltonian has %d qubits (%d after compression).",
            label,
            full_h.num_qubits,
            compressed_h.num_qubits,
        )

        # VQE setup
        iterations: list[int] = []
        energies: list[float] = []

        def _callback(
            eval_count: int,
            _params: np.ndarray,
            mean: float,
            _std: dict[str, Any],
        ) -> None:
            iterations.append(eval_count)
            energies.append(mean)

        sampler, _ = get_sampler()
        ansatz = real_amplitudes(num_qubits=int(compressed_h.num_qubits), reps=1)
        vqe = SamplingVQE(
            sampler=sampler,
            ansatz=ansatz,
            optimizer=COBYLA(maxiter=self.vqe_max_iter),
            aggregation=self.vqe_aggregation,
            callback=_callback,
        )

        raw: SamplingMinimumEigensolverResult = vqe.compute_minimum_eigenvalue(
            compressed_h
        )

        # Extract results
        best = raw.best_measurement or {}
        minimum_energy: float = float(np.real(best.get("value", float("nan"))))
        best_bitstring: str = best.get("bitstring", "")

        # Build probability distribution from quasi-distribution.
        # qiskit-algorithms may return keys as int (older) or str (newer).
        state_probs: dict[str, float] = {}
        if raw.eigenstate is not None:
            quasi_dist = raw.eigenstate
            n_bits = int(compressed_h.num_qubits)
            for key, prob in quasi_dist.items():
                if prob > 0:
                    if isinstance(key, int):
                        bs = format(key, f"0{n_bits}b")
                    else:
                        # Already a bitstring; normalise to n_bits width
                        bs = str(key).zfill(n_bits)
                    state_probs[bs] = float(prob)

        logger.info(
            "Scenario '%s' complete: E_min=%.6f, best_state=%s",
            label,
            minimum_energy,
            best_bitstring,
        )

        return ScenarioResult(
            label=label,
            field=field,
            minimum_energy=minimum_energy,
            best_bitstring=best_bitstring,
            state_probabilities=state_probs,
            vqe_iterations=iterations,
            vqe_energies=energies,
        )

    # ------------------------------------------------------------------
    # Internal helpers - plotting
    # ------------------------------------------------------------------

    def _plot_energy_vs_lambda(
        self,
        plt: Any,
        ticker: Any,
        palette: list[str],
        output_dir: Path | None,
    ) -> None:
        """Plot minimum VQE energy vs. uniform-field strength λ.

        Args:
            plt: matplotlib.pyplot module.
            ticker: matplotlib.ticker module.
            palette (list[str]): Hex colour palette.
            output_dir (Path | None): Save directory, or None to display.

        """
        # Collect uniform-sweep results (skip baseline and non-uniform)
        uniform_results = [r for r in self.results if r.label.startswith("lambda=")]
        baseline = next(
            (r for r in self.results if r.label.startswith("baseline")), None
        )
        nu_result = next((r for r in self.results if "non-uniform" in r.label), None)

        if not uniform_results:
            logger.warning(
                "No uniform-sweep results found; skipping energy-vs-lambda plot."
            )
            return

        lambdas = [float(r.label.split("=")[1].split(" ")[0]) for r in uniform_results]
        energies = [r.minimum_energy for r in uniform_results]

        fig, ax = plt.subplots(figsize=(8, 5))
        fig.patch.set_facecolor("#0f1117")
        ax.set_facecolor("#1a1d27")

        ax.plot(
            lambdas,
            energies,
            marker="o",
            linewidth=2.5,
            markersize=8,
            color=palette[0],
            label="Uniform field",
            zorder=3,
        )
        ax.fill_between(lambdas, energies, alpha=0.15, color=palette[0])

        # Reference lines
        if baseline is not None:
            ax.axhline(
                baseline.minimum_energy,
                linestyle="--",
                linewidth=1.8,
                color=palette[1],
                label=f"Baseline (no field)  E={baseline.minimum_energy:.4f}",
                alpha=0.85,
            )
        if nu_result is not None:
            ax.axhline(
                nu_result.minimum_energy,
                linestyle=":",
                linewidth=1.8,
                color=palette[2],
                label=f"Non-uniform (centre)  E={nu_result.minimum_energy:.4f}",
                alpha=0.85,
            )

        # Styling
        _style_axes(ax, ticker)
        ax.set_xlabel("Field strength λ", fontsize=13, color="#e0e0e0")
        ax.set_ylabel("Minimum VQE energy", fontsize=13, color="#e0e0e0")
        ax.set_title(
            f"Energy vs. uniform field strength — {self.main_chain} "
            f"({self.interaction_type.name})",
            fontsize=14,
            color="#ffffff",
            pad=12,
        )
        ax.legend(
            fontsize=10, facecolor="#1a1d27", edgecolor="#444", labelcolor="#e0e0e0"
        )

        _save_or_show(fig, plt, output_dir, "energy_vs_lambda.png")

    def _plot_energy_bar(
        self,
        plt: Any,
        palette: list[str],
        output_dir: Path | None,
    ) -> None:
        """Plot a bar chart comparing minimum energy across all scenarios.

        Args:
            plt: matplotlib.pyplot module.
            palette (list[str]): Hex colour palette.
            output_dir (Path | None): Save directory, or None to display.

        """
        labels = [r.label for r in self.results]
        energies = [r.minimum_energy for r in self.results]
        colours = [palette[i % len(palette)] for i in range(len(self.results))]

        fig, ax = plt.subplots(figsize=(max(8.0, len(labels) * 1.4), 5))
        fig.patch.set_facecolor("#0f1117")
        ax.set_facecolor("#1a1d27")

        bars = ax.bar(
            range(len(labels)),
            energies,
            color=colours,
            edgecolor="#2a2d3a",
            linewidth=0.8,
            width=0.65,
            zorder=3,
        )

        # Annotate bars with value
        for bar, e in zip(bars, energies, strict=False):
            ax.text(
                bar.get_x() + bar.get_width() / 2,
                bar.get_height() + (max(energies) - min(energies)) * 0.01,
                f"{e:.4f}",
                ha="center",
                va="bottom",
                fontsize=8,
                color="#e0e0e0",
            )

        ax.set_xticks(range(len(labels)))
        ax.set_xticklabels(labels, rotation=30, ha="right", fontsize=9, color="#c0c0c0")
        ax.yaxis.set_tick_params(labelcolor="#c0c0c0")
        ax.set_ylabel("Minimum VQE energy", fontsize=13, color="#e0e0e0")
        ax.set_title(
            f"Energy comparison across field scenarios — {self.main_chain} "
            f"({self.interaction_type.name})",
            fontsize=13,
            color="#ffffff",
            pad=12,
        )
        ax.set_facecolor("#1a1d27")
        ax.grid(axis="y", color="#333", linewidth=0.6, zorder=0)
        ax.spines[:].set_edgecolor("#333")

        fig.tight_layout()
        _save_or_show(fig, plt, output_dir, "energy_comparison_bar.png")

    def _plot_probability_distributions(
        self,
        plt: Any,
        palette: list[str],
        output_dir: Path | None,
        top_k: int = 8,
    ) -> None:
        """Plot top-k bitstring probabilities for each scenario.

        Args:
            plt: matplotlib.pyplot module.
            palette (list[str]): Hex colour palette.
            output_dir (Path | None): Save directory, or None to display.
            top_k (int, optional): Number of highest-probability states to show
                per scenario.  Defaults to 8.

        """
        n = len(self.results)
        fig, axes = plt.subplots(
            n, 1, figsize=(10, 3.5 * n), squeeze=False, sharex=False
        )
        fig.patch.set_facecolor("#0f1117")

        for ax_row, result, colour in zip(
            axes, self.results, palette * (n // len(palette) + 1), strict=False
        ):
            ax = ax_row[0]
            ax.set_facecolor("#1a1d27")

            probs = result.state_probabilities
            if not probs:
                ax.text(
                    0.5,
                    0.5,
                    "No distribution data",
                    ha="center",
                    va="center",
                    color="#888",
                    fontsize=11,
                    transform=ax.transAxes,
                )
                ax.set_title(result.label, color="#ffffff", fontsize=11)
                continue

            # Sort by probability descending and keep top-k
            sorted_states = sorted(probs.items(), key=lambda x: x[1], reverse=True)[
                :top_k
            ]
            states, probs_vals = zip(*sorted_states, strict=False)

            bars = ax.bar(
                range(len(states)),
                probs_vals,
                color=colour,
                edgecolor="#2a2d3a",
                linewidth=0.6,
                width=0.7,
                zorder=3,
            )

            # Highlight the best bitstring
            for i, state in enumerate(states):
                if state == result.best_bitstring:
                    bars[i].set_edgecolor("#ffffff")
                    bars[i].set_linewidth(2.0)

            ax.set_xticks(range(len(states)))
            ax.set_xticklabels(
                [f"|{s}⟩" for s in states],
                rotation=45,
                ha="right",
                fontsize=8,
                color="#c0c0c0",
            )
            ax.yaxis.set_tick_params(labelcolor="#c0c0c0")
            ax.set_ylabel("Probability", fontsize=10, color="#e0e0e0")
            ax.set_title(
                f"{result.label}  |  E_min={result.minimum_energy:.4f}  |  "
                f"best: |{result.best_bitstring}⟩",
                fontsize=10,
                color="#ffffff",
                pad=6,
            )
            ax.grid(axis="y", color="#333", linewidth=0.5, zorder=0)
            ax.spines[:].set_edgecolor("#333")

        fig.suptitle(
            f"State probability distributions — {self.main_chain} ({self.interaction_type.name})",
            fontsize=14,
            color="#ffffff",
            y=1.01,
        )
        fig.tight_layout()
        _save_or_show(fig, plt, output_dir, "probability_distributions.png")

    # ------------------------------------------------------------------
    # Convenience: summary table
    # ------------------------------------------------------------------

    def summary(self) -> str:
        """Return a formatted text summary of all scenario results.

        Returns:
            str: Multi-line table with scenario label, minimum energy, and best
            bitstring for each completed run.

        Raises:
            RuntimeError: If :meth:`run` has not been called yet.

        """
        if not self.results:
            msg = "No results available. Call run() first."
            raise RuntimeError(msg)

        col_label = max(len(r.label) for r in self.results) + 2
        header = f"{'Scenario':<{col_label}}  {'E_min':>12}  {'Best bitstring'}"
        separator = "-" * (col_label + 30)
        rows = [header, separator]
        for r in self.results:
            rows.append(
                f"{r.label:<{col_label}}  {r.minimum_energy:>12.6f}  {r.best_bitstring}"
            )
        return "\n".join(rows)


# ---------------------------------------------------------------------------
# Private plot helpers
# ---------------------------------------------------------------------------


def _style_axes(ax: Any, ticker: Any) -> None:
    """Apply dark-theme styling to a matplotlib Axes object.

    Args:
        ax: The matplotlib Axes to style.
        ticker: The matplotlib.ticker module.

    """
    ax.tick_params(colors="#c0c0c0", which="both")
    ax.xaxis.set_minor_locator(ticker.AutoMinorLocator())
    ax.yaxis.set_minor_locator(ticker.AutoMinorLocator())
    ax.grid(True, which="major", color="#2a2d3a", linewidth=0.8, zorder=0)
    ax.grid(True, which="minor", color="#222530", linewidth=0.4, zorder=0)
    ax.spines[:].set_edgecolor("#333")


def _save_or_show(
    fig: Any,
    plt: Any,
    output_dir: Path | None,
    filename: str,
) -> None:
    """Save figure to *output_dir/filename* or display if *output_dir* is None.

    Args:
        fig: matplotlib Figure object.
        plt: matplotlib.pyplot module.
        output_dir (Path | None): Target directory, or None to show.
        filename (str): Output file name (PNG).

    """
    if output_dir is not None:
        filepath = output_dir / filename
        fig.savefig(
            filepath, dpi=150, bbox_inches="tight", facecolor=fig.get_facecolor()
        )
        logger.info("Saved plot: %s", filepath)
        plt.close(fig)
    else:
        plt.tight_layout()
