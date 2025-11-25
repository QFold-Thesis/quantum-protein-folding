"""
Benchmark script for comparing different optimizers for VQE protein folding.

Tests the following optimizers with 100 iterations:
- COBYLA (Constrained Optimization BY Linear Approximations)
- SLSQP (Sequential Least Squares Programming)
- L_BFGS_B (Limited-memory BFGS with bounds)
- SPSA (Simultaneous Perturbation Stochastic Approximation)
- NFT (Nakanishi-Fujii-Todo optimizer)

Metrics collected:
- Execution time (wall-clock)
- Final eigenvalue (energy)
- Convergence history (iterations vs energy)
- Number of function evaluations
"""

from __future__ import annotations

import csv
import json
import time
from typing import TYPE_CHECKING, Any

from constants import EMPTY_SIDECHAIN_PLACEHOLDER, OUTPUT_DATA_DIR
from logger import get_logger

if TYPE_CHECKING:
    from pathlib import Path

    from qiskit_algorithms.optimizers import Optimizer

logger = get_logger()


def create_optimizer(name: str, maxiter: int = 100) -> Optimizer:
    """
    Create an optimizer instance by name.

    Args:
        name: Optimizer name (COBYLA, SLSQP, L_BFGS_B, SPSA, NFT)
        maxiter: Maximum number of iterations

    Returns:
        Optimizer instance configured for VQE

    """
    from qiskit_algorithms.optimizers import COBYLA, L_BFGS_B, NFT, SLSQP, SPSA

    if name == "COBYLA":
        return COBYLA(maxiter=maxiter)
    if name == "SLSQP":
        return SLSQP(maxiter=maxiter)
    if name == "L_BFGS_B":
        return L_BFGS_B(maxiter=maxiter)
    if name == "SPSA":
        return SPSA(maxiter=maxiter)
    if name == "NFT":
        return NFT(maxiter=maxiter)

    msg = f"Unknown optimizer: {name}"
    raise ValueError(msg)


def run_single_benchmark(
    optimizer_name: str,
    run_number: int,
    protein_sequence: str,
    maxiter: int = 100,
) -> dict[str, Any]:
    """
    Run a single benchmark for one optimizer.

    Args:
        optimizer_name: Name of the optimizer to test
        run_number: Run number (for repeated experiments)
        protein_sequence: Main chain protein sequence
        maxiter: Maximum iterations for optimizer

    Returns:
        Dictionary with benchmark results

    """
    from utils.setup_utils import (
        build_and_compress_hamiltonian,
        run_vqe_optimization,
        setup_folding_system,
        setup_vqe_optimization,
    )

    logger.info("Starting %s run %d", optimizer_name, run_number + 1)

    side_chain = EMPTY_SIDECHAIN_PLACEHOLDER * len(protein_sequence)

    protein, interaction, contact_map, distance_map = setup_folding_system(
        main_chain=protein_sequence, side_chain=side_chain
    )

    _, compressed_h = build_and_compress_hamiltonian(
        protein=protein,
        interaction=interaction,
        contact_map=contact_map,
        distance_map=distance_map,
    )

    optimizer = create_optimizer(optimizer_name, maxiter=maxiter)
    vqe, counts, values = setup_vqe_optimization(
        num_qubits=compressed_h.num_qubits,
        optimizer=optimizer,
        maxiter=maxiter,
    )

    start_time = time.perf_counter()
    try:
        raw_results = run_vqe_optimization(vqe=vqe, hamiltonian=compressed_h)
        elapsed_time = time.perf_counter() - start_time
        eigenvalue = float(raw_results.eigenvalue)

        # Extract best measurement information
        best_measurement = raw_results.best_measurement
        if best_measurement is not None:
            bitstring = best_measurement.get("bitstring")
            state = best_measurement.get("state")
            probability = best_measurement.get("probability")
            value = str(best_measurement.get("value"))
        else:
            bitstring = None
            state = None
            probability = None
            value = None

        success = True
        error_msg = None
    except Exception as e:
        elapsed_time = time.perf_counter() - start_time
        eigenvalue = None
        bitstring = None
        state = None
        probability = None
        value = None
        success = False
        error_msg = str(e)
        logger.exception("Optimizer %s run %d failed", optimizer_name, run_number + 1)

    result = {
        "optimizer": optimizer_name,
        "run": run_number,
        "protein_sequence": protein_sequence,
        "maxiter": maxiter,
        "time_seconds": elapsed_time,
        "eigenvalue": eigenvalue,
        "best_bitstring": bitstring,
        "best_state": state,
        "best_probability": probability,
        "best_value": value,
        "success": success,
        "error": error_msg,
        "num_evaluations": counts[-1] if counts else 0,
        "convergence_history": {
            "iterations": counts,
            "energies": values,
        },
    }

    logger.info(
        "Completed %s run %d: time=%.2fs, eigenvalue=%s",
        optimizer_name, run_number + 1, elapsed_time, eigenvalue
    )

    return result


def run_benchmark_suite(
    protein_sequence: str = "APRLRFY",
    optimizers: list[str] | None = None,
    num_runs: int = 10,
    maxiter: int = 100,
) -> list[dict[str, Any]]:
    """
    Run complete benchmark suite for multiple optimizers.

    Args:
        protein_sequence: Main chain protein sequence to test
        optimizers: List of optimizer names to test (default: all 5)
        num_runs: Number of repeated runs per optimizer
        maxiter: Maximum iterations per optimizer

    Returns:
        List of all benchmark results

    """
    if optimizers is None:
        optimizers = ["COBYLA", "SLSQP", "L_BFGS_B", "SPSA", "NFT"]

    logger.info("=" * 80)
    logger.info("Starting Optimizer Benchmark Suite")
    logger.info("Protein sequence: %s", protein_sequence)
    logger.info("Optimizers to test: %s", optimizers)
    logger.info("Runs per optimizer: %d", num_runs)
    logger.info("Max iterations: %d", maxiter)
    logger.info("=" * 80)

    all_results = []

    for optimizer_name in optimizers:
        logger.info("\n--- Testing optimizer: %s ---", optimizer_name)
        for run_num in range(num_runs):
            result = run_single_benchmark(
                optimizer_name=optimizer_name,
                run_number=run_num,
                protein_sequence=protein_sequence,
                maxiter=maxiter,
            )
            all_results.append(result)

    return all_results


def save_results(results: list[dict[str, Any]], output_dir: Path | None = None) -> None:
    """
    Save benchmark results to JSON and CSV files.

    Args:
        results: List of benchmark results
        output_dir: Output directory (default: OUTPUT_DATA_DIR/benchmarks)

    """
    if output_dir is None:
        output_dir = OUTPUT_DATA_DIR / "benchmarks"

    output_dir.mkdir(parents=True, exist_ok=True)

    timestamp = time.strftime("%Y%m%d_%H%M%S")

    # Save full results as JSON (including convergence history)
    json_path = output_dir / f"benchmark_results_{timestamp}.json"
    with json_path.open("w") as f:
        json.dump(results, f, indent=2)
    logger.info("Saved full results to: %s", json_path)

    # Save summary as CSV (without convergence history)
    csv_path = output_dir / f"benchmark_summary_{timestamp}.csv"
    if results:
        fieldnames = [
            "optimizer",
            "run",
            "protein_sequence",
            "maxiter",
            "time_seconds",
            "eigenvalue",
            "best_bitstring",
            "best_state",
            "best_probability",
            "best_value",
            "num_evaluations",
            "success",
            "error",
        ]
        with csv_path.open("w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            for result in results:
                row = {k: result[k] for k in fieldnames if k in result}
                writer.writerow(row)
        logger.info("Saved summary to: %s", csv_path)


def print_summary(results: list[dict[str, Any]]) -> None:
    """
    Print statistical summary of benchmark results.

    Args:
        results: List of benchmark results

    """
    import statistics

    separator = "=" * 80
    logger.info("\n%s", separator)
    logger.info("BENCHMARK SUMMARY")
    logger.info("%s", separator)

    by_optimizer: dict[str, list[dict[str, Any]]] = {}
    for result in results:
        opt_name = result["optimizer"]
        if opt_name not in by_optimizer:
            by_optimizer[opt_name] = []
        by_optimizer[opt_name].append(result)

    for opt_name, opt_results in by_optimizer.items():
        successful = [r for r in opt_results if r["success"]]

        if not successful:
            logger.info("\n%s: All runs failed", opt_name)
            continue

        times = [r["time_seconds"] for r in successful]
        eigenvalues = [r["eigenvalue"] for r in successful]
        evaluations = [r["num_evaluations"] for r in successful]

        logger.info("\n%s:", opt_name)
        logger.info("  Successful runs: %d/%d", len(successful), len(opt_results))
        logger.info("  Time (s):       mean=%.2f, median=%.2f, std=%.2f",
                   statistics.mean(times), statistics.median(times),
                   statistics.stdev(times) if len(times) > 1 else 0)
        logger.info("  Eigenvalue:     mean=%.6f, best=%.6f, std=%.6f",
                   statistics.mean(eigenvalues), min(eigenvalues),
                   statistics.stdev(eigenvalues) if len(eigenvalues) > 1 else 0)
        logger.info("  Evaluations:    mean=%.0f, median=%.0f",
                   statistics.mean(evaluations), statistics.median(evaluations))

    separator_end = "=" * 80
    logger.info("\n%s", separator_end)


def main() -> None:
    """Main entry point for optimizer benchmark."""
    protein_sequence = "APRLRFY"
    optimizers = ["COBYLA", "SLSQP", "L_BFGS_B", "SPSA", "NFT"]
    optimizers = ["COBYLA"]
    num_runs = 2
    max_iterations = 100

    results = run_benchmark_suite(
        protein_sequence=protein_sequence,
        optimizers=optimizers,
        num_runs=num_runs,
        maxiter=max_iterations,
    )

    save_results(results)
    print_summary(results)

    logger.info("\nBenchmark suite completed!")


if __name__ == "__main__":
    main()
