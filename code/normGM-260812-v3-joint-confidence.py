#!/usr/bin/env python3
"""Conservative joint one-sided Vidal-tail certificates from count data.

The script implements the two statistical interfaces stated in
normGM-260812-v3.tex:

1. held_out_selected_reference: a target is frozen after calibration and its
   fidelity is tested on fresh Bernoulli rounds;
2. calibrated_unknown_spectrum: multinomial mode counts upper-bound the
   cumulative target weights while Bernoulli rounds lower-bound fidelity.

It deliberately rejects adaptive target selection and fidelity testing on
overlapping data. Such an analysis needs a selection-uniform joint confidence
region supplied by the experiment-specific model.
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path
from typing import Any


SCHEMA_VERSION = "1.0"
ALLOWED_MODES = {
    "held_out_selected_reference",
    "calibrated_unknown_spectrum",
    "precomputed_joint_region",
}


class InputError(ValueError):
    """Raised when an input would not support the claimed certificate."""


def require_mapping(value: Any, name: str) -> dict[str, Any]:
    if not isinstance(value, dict):
        raise InputError(f"{name} must be a JSON object")
    return value


def require_string(value: Any, name: str) -> str:
    if not isinstance(value, str) or not value.strip():
        raise InputError(f"{name} must be a nonempty string")
    return value.strip()


def require_probability(value: Any, name: str, *, allow_zero: bool = True) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise InputError(f"{name} must be numeric")
    result = float(value)
    lower_ok = result >= 0.0 if allow_zero else result > 0.0
    if not lower_ok or result >= 1.0 or not math.isfinite(result):
        interval = "[0,1)" if allow_zero else "(0,1)"
        raise InputError(f"{name} must lie in {interval}")
    return result


def require_unit_interval(value: Any, name: str) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise InputError(f"{name} must be numeric")
    result = float(value)
    if not 0.0 <= result <= 1.0 or not math.isfinite(result):
        raise InputError(f"{name} must lie in [0,1]")
    return result


def require_count(value: Any, name: str) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or value < 0:
        raise InputError(f"{name} must be a nonnegative integer")
    return value


def require_counts(value: Any, name: str) -> list[int]:
    if not isinstance(value, list) or len(value) < 2:
        raise InputError(f"{name} must contain at least two integer counts")
    counts = [require_count(item, f"{name}[{index}]") for index, item in enumerate(value)]
    if sum(counts) <= 0:
        raise InputError(f"{name} must have a positive total")
    return counts


def require_spectrum(value: Any, name: str) -> list[float]:
    if not isinstance(value, list) or len(value) < 2:
        raise InputError(f"{name} must contain at least two probabilities")
    spectrum = [require_unit_interval(item, f"{name}[{index}]") for index, item in enumerate(value)]
    total = math.fsum(spectrum)
    if not math.isclose(total, 1.0, rel_tol=0.0, abs_tol=1e-10):
        raise InputError(f"{name} must sum to one; received {total:.16g}")
    return sorted(spectrum, reverse=True)


def require_thresholds(value: Any, dimension: int) -> list[int]:
    if not isinstance(value, list) or not value:
        raise InputError("thresholds must be a nonempty list")
    thresholds: list[int] = []
    for index, item in enumerate(value):
        if isinstance(item, bool) or not isinstance(item, int):
            raise InputError(f"thresholds[{index}] must be an integer")
        if not 2 <= item <= dimension:
            raise InputError(f"thresholds[{index}] must lie in [2,{dimension}]")
        thresholds.append(item)
    if len(set(thresholds)) != len(thresholds):
        raise InputError("thresholds must not contain duplicates")
    return sorted(thresholds)


def require_true(value: Any, name: str) -> None:
    if value is not True:
        raise InputError(f"{name} must be explicitly set to true after verification")


def validate_assumptions(data: dict[str, Any], mode: str) -> dict[str, Any]:
    checks = require_mapping(data.get("assumption_checks"), "assumption_checks")
    require_true(checks.get("sampling_model_justified"), "assumption_checks.sampling_model_justified")
    require_true(checks.get("stopping_rule_documented"), "assumption_checks.stopping_rule_documented")
    require_true(checks.get("measurement_operator_validated"), "assumption_checks.measurement_operator_validated")
    scope = checks.get("dimension_scope")
    if scope not in {"certified_full_state", "normalized_in_subspace_state"}:
        raise InputError(
            "assumption_checks.dimension_scope must be 'certified_full_state' "
            "or 'normalized_in_subspace_state'"
        )
    if mode == "held_out_selected_reference":
        require_true(checks.get("target_frozen_before_test"), "assumption_checks.target_frozen_before_test")
        require_true(
            checks.get("calibration_test_independence_justified"),
            "assumption_checks.calibration_test_independence_justified",
        )
    elif mode == "calibrated_unknown_spectrum":
        require_true(checks.get("reference_fixed_before_data"), "assumption_checks.reference_fixed_before_data")
        require_true(checks.get("schmidt_mode_mapping_validated"), "assumption_checks.schmidt_mode_mapping_validated")
    return checks


def log_binomial_coefficient(n: int, k: int) -> float:
    if not 0 <= k <= n:
        raise InputError("invalid binomial coefficient indices")
    return math.lgamma(n + 1) - math.lgamma(k + 1) - math.lgamma(n - k + 1)


def hoeffding_lower(successes: int, trials: int, epsilon: float) -> tuple[float, float, float]:
    if trials <= 0 or successes > trials:
        raise InputError("fidelity counts require 0 <= successes <= trials and trials > 0")
    estimate = successes / trials
    radius = math.sqrt(math.log(1.0 / epsilon) / (2.0 * trials))
    return estimate, radius, max(0.0, estimate - radius)


def tail_boundary(a_value: float, fidelity: float) -> float:
    a_value = min(1.0, max(0.0, a_value))
    fidelity = min(1.0, max(0.0, fidelity))
    if fidelity <= a_value or a_value >= 1.0:
        return 0.0
    b_value = 1.0 - a_value
    return (
        math.sqrt(b_value * fidelity)
        - math.sqrt(a_value * (1.0 - fidelity))
    ) ** 2


def parse_fidelity_counts(data: dict[str, Any]) -> tuple[int, int, float, float, float]:
    counts = require_mapping(data.get("fidelity_counts"), "fidelity_counts")
    successes = require_count(counts.get("successes"), "fidelity_counts.successes")
    trials = require_count(counts.get("trials"), "fidelity_counts.trials")
    epsilon = require_probability(
        counts.get("failure_probability"),
        "fidelity_counts.failure_probability",
        allow_zero=False,
    )
    estimate, radius, lower = hoeffding_lower(successes, trials, epsilon)
    operator_error = require_unit_interval(data.get("projector_operator_error", 0.0), "projector_operator_error")
    return successes, trials, epsilon, estimate, max(0.0, lower - operator_error)


def base_output(data: dict[str, Any], mode: str) -> dict[str, Any]:
    provenance = require_mapping(data.get("data_provenance"), "data_provenance")
    source = require_string(provenance.get("source"), "data_provenance.source")
    calibration_record_id = require_string(
        provenance.get("calibration_record_id"),
        "data_provenance.calibration_record_id",
    )
    test_record_id = require_string(
        provenance.get("test_record_id"),
        "data_provenance.test_record_id",
    )
    raw_file_sha256 = require_string(
        provenance.get("raw_file_sha256"),
        "data_provenance.raw_file_sha256",
    )
    if len(raw_file_sha256) != 64 or any(character not in "0123456789abcdefABCDEF" for character in raw_file_sha256):
        raise InputError("data_provenance.raw_file_sha256 must contain 64 hexadecimal characters")
    return {
        "schema_version": SCHEMA_VERSION,
        "analysis_id": require_string(data.get("analysis_id"), "analysis_id"),
        "analysis_mode": mode,
        "data_provenance": {
            "source": source,
            "calibration_record_id": calibration_record_id,
            "test_record_id": test_record_id,
            "raw_file_sha256": raw_file_sha256.lower(),
        },
        "warnings": [],
    }


def held_out_analysis(data: dict[str, Any]) -> dict[str, Any]:
    protocol = data.get("selection_protocol")
    if protocol not in {"calibration_then_independent_holdout", "prespecified_reference"}:
        raise InputError(
            "held_out_selected_reference requires selection_protocol equal to "
            "'calibration_then_independent_holdout' or 'prespecified_reference'"
        )
    spectrum = require_spectrum(data.get("selected_target_probabilities"), "selected_target_probabilities")
    thresholds = require_thresholds(data.get("thresholds"), len(spectrum))
    checks = validate_assumptions(data, "held_out_selected_reference")
    successes, trials, epsilon_f, estimate, fidelity_lower = parse_fidelity_counts(data)

    output = base_output(data, "held_out_selected_reference")
    output.update(
        {
            "selection_protocol": protocol,
            "coverage_lower_bound": 1.0 - epsilon_f,
            "simultaneous_over_thresholds": True,
            "fidelity": {
                "successes": successes,
                "trials": trials,
                "estimate": estimate,
                "failure_probability": epsilon_f,
                "operator_error": float(data.get("projector_operator_error", 0.0)),
                "lower_bound": fidelity_lower,
                "method": "one-sided Hoeffding after target freeze",
            },
            "threshold_results": [],
            "assumptions": [
                "The selected target was frozen before the test record was acquired.",
                "Conditional on calibration, the test outcomes are IID Bernoulli trials for one fixed test state.",
                "The supplied target probabilities exactly define the declared post-calibration target.",
                "Any stated projector operator error is a valid uniform operator-norm bound.",
            ],
            "dimension_scope": checks["dimension_scope"],
        }
    )
    for threshold in thresholds:
        a_value = math.fsum(spectrum[: threshold - 1])
        output["threshold_results"].append(
            {
                "r": threshold,
                "A": a_value,
                "A_upper": a_value,
                "fidelity_lower": fidelity_lower,
                "vidal_tail_lower": tail_boundary(a_value, fidelity_lower),
                "activated": fidelity_lower > a_value,
            }
        )
    return output


def calibrated_spectrum_analysis(data: dict[str, Any]) -> dict[str, Any]:
    protocol = data.get("selection_protocol")
    if protocol != "fixed_reference_no_adaptive_selection":
        raise InputError(
            "calibrated_unknown_spectrum requires selection_protocol equal to "
            "'fixed_reference_no_adaptive_selection'"
        )
    spectrum_counts = require_counts(data.get("spectrum_counts"), "spectrum_counts")
    dimension = len(spectrum_counts)
    thresholds = require_thresholds(data.get("thresholds"), dimension)
    checks = validate_assumptions(data, "calibrated_unknown_spectrum")
    epsilon_a = require_probability(
        data.get("spectrum_failure_probability"),
        "spectrum_failure_probability",
        allow_zero=False,
    )
    successes, trials, epsilon_f, estimate, fidelity_lower = parse_fidelity_counts(data)
    if epsilon_a + epsilon_f >= 1.0:
        raise InputError("the spectrum and fidelity failure probabilities must sum to less than one")

    total = sum(spectrum_counts)
    empirical = sorted((count / total for count in spectrum_counts), reverse=True)
    threshold_count = len(thresholds)
    output = base_output(data, "calibrated_unknown_spectrum")
    output.update(
        {
            "selection_protocol": protocol,
            "coverage_lower_bound": 1.0 - epsilon_a - epsilon_f,
            "simultaneous_over_thresholds": True,
            "fidelity": {
                "successes": successes,
                "trials": trials,
                "estimate": estimate,
                "failure_probability": epsilon_f,
                "operator_error": float(data.get("projector_operator_error", 0.0)),
                "lower_bound": fidelity_lower,
                "method": "one-sided Hoeffding",
            },
            "spectrum_calibration": {
                "counts": spectrum_counts,
                "trials": total,
                "failure_probability": epsilon_a,
                "method": "Hoeffding union over mode subsets and requested thresholds",
            },
            "threshold_results": [],
            "assumptions": [
                "The reference is fixed independently of the reported observations; it was not selected adaptively from them.",
                "Spectrum counts are IID multinomial trials whose labels correspond to calibrated Schmidt modes.",
                "Mode counts identify Schmidt probabilities under a validated phase, basis-alignment, and leakage model.",
                "Fidelity counts support a valid IID Bernoulli model for the fixed reference projector.",
                "Any stated projector operator error is a valid uniform operator-norm bound.",
            ],
            "dimension_scope": checks["dimension_scope"],
        }
    )
    for threshold in thresholds:
        subset_size = threshold - 1
        a_hat = math.fsum(empirical[:subset_size])
        radius = math.sqrt(
            (
                log_binomial_coefficient(dimension, subset_size)
                + math.log(threshold_count / epsilon_a)
            )
            / (2.0 * total)
        )
        a_upper = min(1.0, a_hat + radius)
        output["threshold_results"].append(
            {
                "r": threshold,
                "A_estimate": a_hat,
                "A_radius": radius,
                "A_upper": a_upper,
                "fidelity_lower": fidelity_lower,
                "vidal_tail_lower": tail_boundary(a_upper, fidelity_lower),
                "activated": fidelity_lower > a_upper,
            }
        )
    if dimension > 100:
        output["warnings"].append(
            "The subset-union spectrum region can be very conservative at high dimension."
        )
    return output


def precomputed_joint_analysis(data: dict[str, Any]) -> dict[str, Any]:
    region = require_mapping(data.get("joint_region"), "joint_region")
    epsilon = require_probability(
        region.get("failure_probability"),
        "joint_region.failure_probability",
        allow_zero=False,
    )
    fidelity_lower = require_unit_interval(region.get("fidelity_lower"), "joint_region.fidelity_lower")
    a_upper_raw = require_mapping(region.get("A_upper_by_r"), "joint_region.A_upper_by_r")
    if not a_upper_raw:
        raise InputError("joint_region.A_upper_by_r must not be empty")
    results = []
    for key, value in a_upper_raw.items():
        try:
            threshold = int(key)
        except (TypeError, ValueError) as error:
            raise InputError(f"invalid threshold key {key!r}") from error
        if threshold < 2:
            raise InputError("joint-region threshold keys must be at least two")
        a_upper = require_unit_interval(value, f"joint_region.A_upper_by_r[{key!r}]")
        results.append(
            {
                "r": threshold,
                "A_upper": a_upper,
                "fidelity_lower": fidelity_lower,
                "vidal_tail_lower": tail_boundary(a_upper, fidelity_lower),
                "activated": fidelity_lower > a_upper,
            }
        )
    output = base_output(data, "precomputed_joint_region")
    output.update(
        {
            "selection_protocol": "experiment_specific_uniform_joint_region",
            "coverage_lower_bound": 1.0 - epsilon,
            "simultaneous_over_thresholds": True,
            "joint_region_method": require_string(region.get("method"), "joint_region.method"),
            "threshold_results": sorted(results, key=lambda item: item["r"]),
            "assumptions": [
                "The supplied region has simultaneous coverage for fidelity, target selection, and every reported A_r.",
                "The experiment-specific method includes all stopping, calibration, and dependence effects.",
            ],
        }
    )
    output["warnings"].append(
        "This script post-processes the supplied joint region; it does not validate that region's coverage proof."
    )
    return output


def analyze(data: dict[str, Any]) -> dict[str, Any]:
    mode = data.get("analysis_mode")
    if mode not in ALLOWED_MODES:
        raise InputError(f"analysis_mode must be one of {sorted(ALLOWED_MODES)}")
    if (
        mode != "precomputed_joint_region"
        and data.get("selection_protocol") == "overlapping_adaptive_data"
    ):
        raise InputError(
            "ordinary pointwise intervals are invalid after overlapping adaptive selection; "
            "supply a precomputed_joint_region from a selection-uniform analysis"
        )
    if mode == "held_out_selected_reference":
        return held_out_analysis(data)
    if mode == "calibrated_unknown_spectrum":
        return calibrated_spectrum_analysis(data)
    return precomputed_joint_analysis(data)


def self_test() -> None:
    expected = 0.00728390480771747
    actual = tail_boundary(0.75, 0.82)
    if not math.isclose(actual, expected, rel_tol=0.0, abs_tol=1e-12):
        raise AssertionError((actual, expected))
    _, _, lower = hoeffding_lower(410, 500, 0.01)
    if not math.isclose(lower, 0.7521385957558488, rel_tol=0.0, abs_tol=1e-12):
        raise AssertionError(lower)
    payload = {
        "analysis_id": "internal-self-test",
        "analysis_mode": "held_out_selected_reference",
        "selection_protocol": "calibration_then_independent_holdout",
        "data_provenance": {
            "source": "synthetic internal self-test",
            "calibration_record_id": "self-test-calibration",
            "test_record_id": "self-test-test",
            "raw_file_sha256": "0" * 64,
        },
        "assumption_checks": {
            "sampling_model_justified": True,
            "stopping_rule_documented": True,
            "measurement_operator_validated": True,
            "target_frozen_before_test": True,
            "calibration_test_independence_justified": True,
            "dimension_scope": "certified_full_state",
        },
        "selected_target_probabilities": [0.25, 0.25, 0.25, 0.25],
        "thresholds": [2, 3, 4],
        "fidelity_counts": {"successes": 410, "trials": 500, "failure_probability": 0.01},
        "projector_operator_error": 0.0,
    }
    result = analyze(payload)
    if result["coverage_lower_bound"] != 0.99 or len(result["threshold_results"]) != 3:
        raise AssertionError(result)

    spectrum_payload = {
        "analysis_id": "internal-spectrum-self-test",
        "analysis_mode": "calibrated_unknown_spectrum",
        "selection_protocol": "fixed_reference_no_adaptive_selection",
        "data_provenance": {
            "source": "synthetic internal self-test",
            "calibration_record_id": "self-test-spectrum-calibration",
            "test_record_id": "self-test-spectrum-test",
            "raw_file_sha256": "1" * 64,
        },
        "assumption_checks": {
            "sampling_model_justified": True,
            "stopping_rule_documented": True,
            "measurement_operator_validated": True,
            "reference_fixed_before_data": True,
            "schmidt_mode_mapping_validated": True,
            "dimension_scope": "certified_full_state",
        },
        "spectrum_counts": [400, 300, 200, 100],
        "spectrum_failure_probability": 0.01,
        "thresholds": [2, 3, 4],
        "fidelity_counts": {
            "successes": 900,
            "trials": 1000,
            "failure_probability": 0.01,
        },
        "projector_operator_error": 0.0,
    }
    spectrum_result = analyze(spectrum_payload)
    if not math.isclose(spectrum_result["coverage_lower_bound"], 0.98, abs_tol=1e-15):
        raise AssertionError(spectrum_result)
    if any(
        item["A_upper"] < item["A_estimate"]
        for item in spectrum_result["threshold_results"]
    ):
        raise AssertionError(spectrum_result)

    overlapping_payload = dict(payload)
    overlapping_payload["selection_protocol"] = "overlapping_adaptive_data"
    try:
        analyze(overlapping_payload)
    except InputError:
        pass
    else:
        raise AssertionError("overlapping adaptive data were not rejected")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", nargs="?", type=Path, help="JSON count-data input")
    parser.add_argument("--output", type=Path, help="write JSON output to this path")
    parser.add_argument("--self-test", action="store_true", help="run deterministic internal checks")
    args = parser.parse_args()

    try:
        if args.self_test:
            self_test()
            print("self-test passed")
            return 0
        if args.input is None:
            parser.error("input is required unless --self-test is used")
        with args.input.open("r", encoding="utf-8") as handle:
            data = require_mapping(json.load(handle), "root")
        output = analyze(data)
        rendered = json.dumps(output, indent=2, sort_keys=True) + "\n"
        if args.output:
            args.output.write_text(rendered, encoding="utf-8")
        else:
            sys.stdout.write(rendered)
        return 0
    except (InputError, json.JSONDecodeError, OSError) as error:
        print(f"input rejected: {error}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
