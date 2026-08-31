"""Reproduce the published-data screening calculation in Sec. VI.

Source: W.-H. Zhang et al., npj Quantum Information 5, 4 (2019), Fig. 4,
doi:10.1038/s41534-018-0120-0.

The Schmidt-amplitude bar heights were digitized as integer pixels from a
300-dpi rendering of Fig. 4. The script stores those heights explicitly,
normalizes each row in Euclidean norm, and evaluates the fixed-fidelity
Vidal-tail bounds. The source figure provides no calibrated
pixel-to-amplitude uncertainty model, so the digitization is a deterministic
screening input rather than a confidence region. The manuscript rounds
probabilities and derived bounds to three decimal places. Full-precision JSON
values are computational records, not claims of experimental precision. The
script does not reconstruct raw counts or turn the reported standard
deviations into confidence limits. Only the Python standard library is
required.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path


SOURCE_DOI = "10.1038/s41534-018-0120-0"
RENDER_DPI = 300
MANUSCRIPT_DECIMAL_PLACES = 3

# Central teal-bar heights in pixels, ordered as c_0, c_1, c_2, c_3.
DIGITIZED_HEIGHTS = {
    "psi_0": [325, 320, 305, 293],
    "psi_1": [321, 311, 309, 303],
    "psi_2": [323, 309, 309, 303],
    "psi_3": [314, 313, 311, 307],
    "psi_4": [317, 316, 310, 303],
    "psi_5": [324, 239, 179, 151],
    "psi_6": [329, 249, 199, 85],
    "psi_7": [296, 293, 145, 144],
    "psi_8": [314, 142, 141, 2],
    "psi_9": [298, 14, 10, 0],
}

TOMOGRAPHY_FIDELITIES = {
    "psi_0": (0.974, 0.003),
    "psi_1": (0.960, 0.001),
    "psi_2": (0.969, 0.001),
    "psi_3": (0.968, 0.002),
    "psi_4": (0.951, 0.004),
    "psi_5": (0.966, 0.003),
    "psi_6": (0.957, 0.004),
    "psi_7": (0.955, 0.002),
    "psi_8": (0.950, 0.006),
    "psi_9": (0.981, 0.004),
}


def normalized_probabilities(heights: list[int]) -> list[float]:
    squared_norm = sum(height * height for height in heights)
    if squared_norm == 0:
        raise ValueError("at least one digitized height must be positive")
    probabilities = [height * height / squared_norm for height in heights]
    return sorted(probabilities, reverse=True)


def reference_spectrum_boundary(
    probabilities: list[float], r: int, fidelity: float
) -> float:
    head = sum(probabilities[: r - 1])
    if fidelity <= head:
        return 0.0
    tail = 1.0 - head
    return (
        math.sqrt(tail * fidelity)
        - math.sqrt(head * (1.0 - fidelity))
    ) ** 2


def fixed_lambda_boundary(m: int, r: int, spectral_parameter: float) -> float:
    threshold = (r - 1) / m
    if spectral_parameter <= threshold:
        return 0.0
    return (
        math.sqrt((m - r + 1) * spectral_parameter)
        - math.sqrt((r - 1) * (1.0 - spectral_parameter))
    ) ** 2 / m


def largest_coefficient_relaxation(
    probabilities: list[float], r: int, fidelity: float
) -> float:
    m = len(probabilities)
    spectral_parameter = max(fidelity / (m * probabilities[0]), 1.0 / m)
    return fixed_lambda_boundary(m, r, min(spectral_parameter, 1.0))


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output",
        type=Path,
        help="write the JSON result to this UTF-8 file instead of stdout",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    rows = []
    for state, heights in DIGITIZED_HEIGHTS.items():
        probabilities = normalized_probabilities(heights)
        fidelity, standard_deviation = TOMOGRAPHY_FIDELITIES[state]
        screening_fidelity = fidelity - standard_deviation

        adapted = {
            r: reference_spectrum_boundary(
                probabilities, r, screening_fidelity
            )
            for r in range(2, len(probabilities) + 1)
        }
        largest_coefficient = {
            r: largest_coefficient_relaxation(
                probabilities, r, screening_fidelity
            )
            for r in range(2, len(probabilities) + 1)
        }
        activated = [r for r, value in adapted.items() if value > 0.0]
        largest_activated_r = max(activated) if activated else None

        rows.append(
            {
                "state": state,
                "digitized_amplitude_heights": heights,
                "reference_probabilities": probabilities,
                "tomography_fidelity": fidelity,
                "reported_standard_deviation": standard_deviation,
                "screening_fidelity": screening_fidelity,
                "adapted_bounds": adapted,
                "largest_coefficient_bounds": largest_coefficient,
                "largest_activated_r": largest_activated_r,
            }
        )

    output = {
        "source_doi": SOURCE_DOI,
        "render_dpi": RENDER_DPI,
        "manuscript_decimal_places": MANUSCRIPT_DECIMAL_PLACES,
        "digitization_scope": (
            "integer central bar heights from a 300-dpi rendering; no "
            "calibrated pixel-to-amplitude uncertainty model"
        ),
        "statistical_scope": (
            "F_T minus one reported standard deviation is a screening input, "
            "not a one-sided confidence bound"
        ),
        "rows": rows,
    }
    serialized = json.dumps(output, indent=2, sort_keys=True) + "\n"
    if args.output is None:
        print(serialized, end="")
    else:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(serialized, encoding="utf-8")


if __name__ == "__main__":
    main()
