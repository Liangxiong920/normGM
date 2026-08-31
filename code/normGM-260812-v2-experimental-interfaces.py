"""Reproduce the published-data post-processings added in manuscript v2.

The inputs below are copied from the cited articles.  Values obtained by
subtracting a displayed error bar are screening inputs only; they are not
reinterpreted as distribution-free confidence bounds.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path


def boundary_from_threshold(fidelity: float, threshold: float) -> float:
    """Evaluate R_A(F) for A=threshold and B=1-A."""
    if not 0.0 <= fidelity <= 1.0:
        raise ValueError("fidelity must lie in [0,1]")
    if not 0.0 <= threshold <= 1.0:
        raise ValueError("threshold must lie in [0,1]")
    if fidelity <= threshold or threshold == 1.0:
        return 0.0
    return (
        math.sqrt((1.0 - threshold) * fidelity)
        - math.sqrt(threshold * (1.0 - fidelity))
    ) ** 2


def reference_spectrum_boundary(
    fidelity: float, spectrum: list[float], rank_threshold: int
) -> tuple[float, float]:
    """Return (A_r, R_{r,phi}(F)) for an ordered reference spectrum."""
    if not 2 <= rank_threshold <= len(spectrum):
        raise ValueError("rank_threshold must lie in {2,...,m}")
    if any(spectrum[i] < spectrum[i + 1] for i in range(len(spectrum) - 1)):
        raise ValueError("spectrum must be nonincreasing")
    if not math.isclose(sum(spectrum), 1.0, rel_tol=0.0, abs_tol=1e-14):
        raise ValueError("spectrum must sum to one")
    threshold = sum(spectrum[: rank_threshold - 1])
    return threshold, boundary_from_threshold(fidelity, threshold)


def largest_coefficient_relaxation(
    fidelity: float, spectrum: list[float], rank_threshold: int
) -> float:
    """Evaluate the nu_1 relaxation in Corollary 5 of the manuscript."""
    dimension = len(spectrum)
    lambda_eff = max(fidelity / (dimension * spectrum[0]), 1.0 / dimension)
    uniform_threshold = (rank_threshold - 1) / dimension
    return boundary_from_threshold(lambda_eff, uniform_threshold)


def uniform_case(
    dimension: int, fidelity: float, rank_thresholds: list[int]
) -> dict[str, object]:
    spectrum = [1.0 / dimension] * dimension
    rows = {}
    for rank_threshold in rank_thresholds:
        threshold, value = reference_spectrum_boundary(
            fidelity, spectrum, rank_threshold
        )
        rows[str(rank_threshold)] = {"A_r": threshold, "R": value}
    return {"dimension": dimension, "fidelity_input": fidelity, "tails": rows}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--output",
        type=Path,
        help="write the JSON result to this UTF-8 file instead of stdout",
    )
    args = parser.parse_args()

    guo_spectrum = [3.0 / 4.0, 1.0 / 12.0, 1.0 / 12.0, 1.0 / 12.0]
    guo_fidelity = 0.988
    guo_rows = {}
    for rank_threshold in (2, 3, 4):
        threshold, value = reference_spectrum_boundary(
            guo_fidelity, guo_spectrum, rank_threshold
        )
        guo_rows[str(rank_threshold)] = {
            "A_r": threshold,
            "R": value,
            "nu1_relaxation": largest_coefficient_relaxation(
                guo_fidelity, guo_spectrum, rank_threshold
            ),
        }

    output = {
        "sources": {
            "Guo2018GMLE": "10.1103/PhysRevA.97.062309",
            "Hu2020Multipath": "10.1103/PhysRevLett.125.090503",
            "Chang2026TimeFrequency": "10.1126/sciadv.aee1333",
        },
        "guo_nonuniform_target": {
            "reported_fidelity": 0.991,
            "displayed_error": 0.003,
            "screening_input": guo_fidelity,
            "reference_spectrum": guo_spectrum,
            "tails": guo_rows,
        },
        "hu_uniform_d32": {
            "reported_fidelity": 0.933,
            "displayed_error": 0.001,
            "screening": uniform_case(32, 0.932, [2, 16, 30]),
        },
        "chang_uniform_d1021_central": uniform_case(
            1021, 0.654, [2, 512, 668]
        ),
        "chang_uniform_d1021_displayed_lower_endpoint": uniform_case(
            1021, 0.650, [2, 512, 664, 668]
        ),
    }
    serialized = json.dumps(output, indent=2, sort_keys=True) + "\n"
    if args.output is None:
        print(serialized, end="")
    else:
        args.output.write_text(serialized, encoding="utf-8")


if __name__ == "__main__":
    main()
