"""Reproducible numerical checks for normGM-260814-v1.

The pure-state test solves the original coefficient-matrix problem with SLSQP.
It does not provide the optimizer with the analytic two-block spectrum.

The mixed-state benchmark uses a colored-noise family and searches rank-rho HJW
ensembles by generic complex Givens rotations.  Every returned ensemble is an
exact decomposition up to the reported floating-point reconstruction residual,
so its average Vidal tail is a certified numerical upper bound (subject only to
the printed floating-point residual).  It is not claimed to be the convex roof.
"""

from __future__ import annotations

import csv
import json
import math
import sys
from pathlib import Path

import numpy as np
import scipy
from scipy.optimize import minimize
from scipy.stats import beta, binom


ROOT = Path(__file__).resolve().parent
DATA_ROOT = ROOT.parent / "data"
SEED = 20260806
THETA = 0.25
KAPPA = 0.45


def vidal_tail(vector: np.ndarray, dimension: int, threshold: int) -> float:
    singular_values = np.linalg.svd(
        vector.reshape(dimension, dimension), compute_uv=False
    )
    probabilities = np.sort(singular_values**2)[::-1]
    return float(probabilities[threshold - 1 :].sum())


def reference_boundary(fidelity: float, head_weight: float) -> float:
    if fidelity <= head_weight:
        return 0.0
    return (
        math.sqrt((1.0 - head_weight) * fidelity)
        - math.sqrt(head_weight * (1.0 - fidelity))
    ) ** 2


def isotropic_boundary(dimension: int, threshold: int, parameter: float) -> float:
    head_size = threshold - 1
    tail_size = dimension - threshold + 1
    if parameter <= head_size / dimension:
        return 0.0
    return (
        math.sqrt(tail_size * parameter)
        - math.sqrt(head_size * (1.0 - parameter))
    ) ** 2 / dimension


def reference_state(
    dimension: int, nonuniformity: float, angle: float = THETA
) -> tuple[np.ndarray, np.ndarray]:
    probabilities = np.exp(-nonuniformity * np.arange(dimension))
    probabilities /= probabilities.sum()
    rotation = np.eye(dimension)
    cosine, sine = math.cos(angle), math.sin(angle)
    rotation[:2, :2] = np.array([[cosine, -sine], [sine, cosine]])
    vector = np.zeros(dimension * dimension, dtype=complex)
    for index in range(dimension):
        for output in range(dimension):
            vector[index * dimension + output] += (
                math.sqrt(probabilities[index]) * rotation[output, index]
            )
    return vector, probabilities


def colored_noise_state(
    dimension: int, visibility: float, color: float = KAPPA
) -> tuple[np.ndarray, np.ndarray]:
    probabilities = np.exp(-color * np.arange(dimension))
    probabilities /= probabilities.sum()
    maximally_entangled = np.eye(dimension).reshape(-1) / math.sqrt(dimension)
    state = visibility * np.outer(maximally_entangled, maximally_entangled.conj())
    for index, probability in enumerate(probabilities):
        basis_vector = np.zeros(dimension * dimension)
        basis_vector[index * dimension + index] = 1.0
        state += (1.0 - visibility) * probability * np.outer(
            basis_vector, basis_vector
        )
    return state, probabilities


def pure_fixed_fidelity_search(
    reference: np.ndarray,
    dimension: int,
    threshold: int,
    fidelity: float,
    starts: int = 10,
) -> dict[str, float]:
    """Minimize E_r over a general complex coefficient matrix."""

    size = dimension * dimension
    rng = np.random.default_rng(SEED + 100 * dimension + threshold)

    def unpack(parameters: np.ndarray) -> np.ndarray:
        return parameters[:size] + 1j * parameters[size:]

    def objective(parameters: np.ndarray) -> float:
        return vidal_tail(unpack(parameters), dimension, threshold)

    constraints = [
        {
            "type": "eq",
            "fun": lambda parameters: float(
                np.vdot(unpack(parameters), unpack(parameters)).real - 1.0
            ),
        },
        {
            "type": "eq",
            "fun": lambda parameters: float(
                np.vdot(reference, unpack(parameters)).real - math.sqrt(fidelity)
            ),
        },
        {
            "type": "eq",
            "fun": lambda parameters: float(
                np.vdot(reference, unpack(parameters)).imag
            ),
        },
    ]
    best_value = math.inf
    best_residual = math.inf
    successful_starts = 0
    for _ in range(starts):
        orthogonal = rng.normal(size=size) + 1j * rng.normal(size=size)
        orthogonal -= reference * np.vdot(reference, orthogonal)
        orthogonal /= np.linalg.norm(orthogonal)
        vector = (
            math.sqrt(fidelity) * reference
            + math.sqrt(1.0 - fidelity) * orthogonal
        )
        initial = np.concatenate([vector.real, vector.imag])
        result = minimize(
            objective,
            initial,
            method="SLSQP",
            constraints=constraints,
            options={"maxiter": 1800, "ftol": 1.0e-12},
        )
        if not result.success:
            continue
        successful_starts += 1
        candidate = unpack(result.x)
        residual = max(
            abs(np.vdot(candidate, candidate).real - 1.0),
            abs(np.vdot(reference, candidate).real - math.sqrt(fidelity)),
            abs(np.vdot(reference, candidate).imag),
        )
        if result.fun < best_value:
            best_value = float(result.fun)
            best_residual = float(residual)
    if not math.isfinite(best_value):
        raise RuntimeError("all SLSQP starts failed")
    return {
        "numerical_value": best_value,
        "constraint_residual": best_residual,
        "successful_starts": successful_starts,
    }


def haar_unitary(order: int, rng: np.random.Generator) -> np.ndarray:
    matrix = rng.normal(size=(order, order)) + 1j * rng.normal(
        size=(order, order)
    )
    unitary, triangular = np.linalg.qr(matrix)
    diagonal = np.diag(triangular)
    return unitary * (diagonal / np.abs(diagonal)).conj()


def hjw_objective(
    unitary: np.ndarray,
    weighted_eigenvectors: np.ndarray,
    dimension: int,
    threshold: int,
) -> float:
    subnormalized = weighted_eigenvectors @ unitary.T
    total = 0.0
    for column in range(unitary.shape[0]):
        vector = subnormalized[:, column]
        probability = float(np.vdot(vector, vector).real)
        if probability > 1.0e-14:
            total += probability * vidal_tail(
                vector / math.sqrt(probability), dimension, threshold
            )
    return total


def hjw_upper_search(
    state: np.ndarray,
    dimension: int,
    threshold: int,
    seed: int,
    starts: int = 8,
    iterations: int = 2500,
) -> dict[str, float]:
    """Search rank-rho HJW ensembles by generic complex row rotations."""

    eigenvalues, eigenvectors = np.linalg.eigh(state)
    support = eigenvalues > 1.0e-12
    eigenvalues = eigenvalues[support]
    eigenvectors = eigenvectors[:, support]
    rank = len(eigenvalues)
    weighted_eigenvectors = eigenvectors * np.sqrt(eigenvalues)
    rng = np.random.default_rng(seed)
    best_value = math.inf
    best_unitary = np.eye(rank, dtype=complex)
    for start in range(starts):
        unitary = (
            np.eye(rank, dtype=complex)
            if start == 0
            else haar_unitary(rank, rng)
        )
        current = hjw_objective(
            unitary, weighted_eigenvectors, dimension, threshold
        )
        scale = 0.45
        for iteration in range(iterations):
            first, second = rng.choice(rank, 2, replace=False)
            angle = rng.normal(scale=scale)
            phase = rng.uniform(0.0, 2.0 * math.pi)
            cosine, sine = math.cos(angle), math.sin(angle)
            rotation = np.array(
                [
                    [cosine, -np.exp(1j * phase) * sine],
                    [np.exp(-1j * phase) * sine, cosine],
                ]
            )
            candidate = unitary.copy()
            candidate[[first, second], :] = (
                rotation @ unitary[[first, second], :]
            )
            value = hjw_objective(
                candidate, weighted_eigenvectors, dimension, threshold
            )
            if value < current:
                unitary, current = candidate, value
            if (iteration + 1) % 500 == 0:
                scale *= 0.62
        if current < best_value:
            best_value = current
            best_unitary = unitary

    subnormalized = weighted_eigenvectors @ best_unitary.T
    reconstruction = subnormalized @ subnormalized.conj().T
    residual = float(np.linalg.norm(reconstruction - state, ord="fro"))
    return {
        "upper_bound": float(best_value),
        "rank": rank,
        "reconstruction_residual": residual,
    }


def relaxed_weighted_envelope(
    parameter: float,
    costs: list[float],
    dimensions: list[float],
    ambient_dimension: int,
    tolerance: float = 1.0e-13,
) -> float:
    zero_dimension = sum(
        dimension
        for cost, dimension in zip(costs, dimensions)
        if cost == 0.0
    )
    if parameter <= zero_dimension / ambient_dimension:
        return 0.0

    def values(alpha: float) -> tuple[float, float]:
        z1 = sum(
            dimension / (cost + alpha)
            for cost, dimension in zip(costs, dimensions)
        )
        z2 = sum(
            dimension / (cost + alpha) ** 2
            for cost, dimension in zip(costs, dimensions)
        )
        spectral_parameter = z1 * z1 / (ambient_dimension * z2)
        objective = z1 / z2 - alpha
        return spectral_parameter, objective

    low, high = 1.0e-15, 1.0
    while values(high)[0] < parameter:
        high *= 2.0
    while high - low > tolerance * max(1.0, high):
        midpoint = (low + high) / 2.0
        if values(midpoint)[0] < parameter:
            low = midpoint
        else:
            high = midpoint
    return values((low + high) / 2.0)[1]


def write_csv(name: str, rows: list[dict[str, float]]) -> None:
    DATA_ROOT.mkdir(parents=True, exist_ok=True)
    with (DATA_ROOT / name).open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def run_pure_checks() -> list[dict[str, float]]:
    rows: list[dict[str, float]] = []
    for dimension in (3, 4):
        for sample in range(3):
            rng = np.random.default_rng(100 + 10 * dimension + sample)
            probabilities = np.sort(
                rng.dirichlet(np.ones(dimension))
            )[::-1]
            reference = np.diag(np.sqrt(probabilities)).reshape(-1)
            threshold = min(sample + 2, dimension)
            head = float(probabilities[: threshold - 1].sum())
            fidelity = head + 0.55 * (1.0 - head)
            analytic = reference_boundary(fidelity, head)
            result = pure_fixed_fidelity_search(
                reference, dimension, threshold, fidelity
            )
            rows.append(
                {
                    "dimension": dimension,
                    "threshold": threshold,
                    "head_weight": head,
                    "fidelity": fidelity,
                    "analytic_value": analytic,
                    "numerical_value": result["numerical_value"],
                    "absolute_gap": abs(result["numerical_value"] - analytic),
                    "constraint_residual": result["constraint_residual"],
                    "successful_starts": result["successful_starts"],
                }
            )
    write_csv("normGM-260814-v1-pure-check.csv", rows)
    return rows


def mixed_row(
    dimension: int,
    threshold: int,
    nonuniformity: float,
    visibility: float,
    upper_cache: dict[tuple[int, float], dict[str, float]],
) -> dict[str, float]:
    reference, probabilities = reference_state(dimension, nonuniformity)
    state, _ = colored_noise_state(dimension, visibility)
    fidelity = float(np.vdot(reference, state @ reference).real)
    head = float(probabilities[: threshold - 1].sum())
    main_bound = reference_boundary(fidelity, head)
    old_parameter = max(
        fidelity / (dimension * probabilities[0]), 1.0 / dimension
    )
    independent_lower = isotropic_boundary(
        dimension, threshold, old_parameter
    )
    key = (dimension, visibility)
    if key not in upper_cache:
        upper_cache[key] = hjw_upper_search(
            state,
            dimension,
            threshold,
            SEED + 1000 * dimension + int(round(100 * visibility)),
        )
    upper = upper_cache[key]
    return {
        "dimension": dimension,
        "threshold": threshold,
        "nonuniformity": nonuniformity,
        "basis_angle": THETA,
        "color": KAPPA,
        "visibility": visibility,
        "fidelity": fidelity,
        "main_lower_bound": main_bound,
        "nu1_lower_bound": independent_lower,
        "hjw_upper_bound": upper["upper_bound"],
        "lower_upper_gap": upper["upper_bound"] - main_bound,
        "state_rank": upper["rank"],
        "reconstruction_residual": upper["reconstruction_residual"],
    }


def run_mixed_checks() -> tuple[list[dict[str, float]], list[dict[str, float]]]:
    cache: dict[tuple[int, float], dict[str, float]] = {}
    noise_rows = [
        mixed_row(4, 3, 0.2, visibility, cache)
        for visibility in (0.60, 0.70, 0.80, 0.90, 0.95)
    ]
    dimension_rows = [
        mixed_row(dimension, 3, nonuniformity, 0.90, cache)
        for dimension in (3, 4, 5)
        for nonuniformity in (0.2, 0.6, 1.0)
    ]
    write_csv("normGM-260814-v1-colored-noise.csv", noise_rows)
    write_csv("normGM-260814-v1-dimension-scan.csv", dimension_rows)
    return noise_rows, dimension_rows


def run_weighted_comparison() -> dict[str, object]:
    probabilities = [0.4, 0.3, 0.2, 0.1]
    weights = [1.0 / 3.0] * 3
    costs = [0.0, 1.0 / 3.0, 2.0 / 3.0, 1.0]
    grid_size = 20001
    differences: list[float] = []
    rows: list[dict[str, float]] = []
    for index in range(grid_size):
        fidelity = index / (grid_size - 1)
        full = sum(
            weight
            * reference_boundary(
                fidelity, sum(probabilities[: threshold - 1])
            )
            for weight, threshold in zip(weights, (2, 3, 4))
        )
        parameter = max(fidelity / (4.0 * probabilities[0]), 0.25)
        auxiliary = relaxed_weighted_envelope(
            parameter, costs, [1.0] * 4, 4
        )
        difference = full - auxiliary
        differences.append(difference)
        if index % 100 == 0:
            rows.append(
                {
                    "fidelity": fidelity,
                    "full_spectrum": full,
                    "auxiliary": auxiliary,
                    "best": max(full, auxiliary),
                }
            )
    maximum_index = int(np.argmax(differences))
    positive_integral = float(
        np.trapezoid(np.maximum(differences, 0.0), dx=1.0 / (grid_size - 1))
    )
    write_csv("normGM-260814-v1-weighted-comparison.csv", rows)
    fidelity = 0.9256
    full = sum(
        weight
        * reference_boundary(
            fidelity, sum(probabilities[: threshold - 1])
        )
        for weight, threshold in zip(weights, (2, 3, 4))
    )
    auxiliary = relaxed_weighted_envelope(
        max(fidelity / 1.6, 0.25), costs, [1.0] * 4, 4
    )
    secondary_probabilities = [0.28, 0.27, 0.24, 0.21]
    secondary_fidelity = 0.47
    secondary_full = sum(
        weight
        * reference_boundary(
            secondary_fidelity,
            sum(secondary_probabilities[: threshold - 1]),
        )
        for weight, threshold in zip(weights, (2, 3, 4))
    )
    secondary_auxiliary = relaxed_weighted_envelope(
        max(
            secondary_fidelity / (4.0 * secondary_probabilities[0]),
            0.25,
        ),
        costs,
        [1.0] * 4,
        4,
    )
    delayed_weights = [0.0, 0.5, 0.5]
    delayed_costs = [0.0, 0.0, 0.5, 1.0]
    delayed_fidelity = 0.75
    delayed_full = sum(
        weight
        * reference_boundary(
            delayed_fidelity, sum(probabilities[: threshold - 1])
        )
        for weight, threshold in zip(delayed_weights, (2, 3, 4))
    )
    delayed_auxiliary = relaxed_weighted_envelope(
        max(delayed_fidelity / 1.6, 0.25),
        delayed_costs,
        [1.0] * 4,
        4,
    )
    return {
        "display_fidelity": fidelity,
        "display_full_spectrum": full,
        "display_auxiliary": auxiliary,
        "maximum_positive_difference": differences[maximum_index],
        "maximum_difference_fidelity": maximum_index / (grid_size - 1),
        "positive_difference_integral": positive_integral,
        "grid_size": grid_size,
        "auxiliary_advantage_example": {
            "probabilities": secondary_probabilities,
            "fidelity": secondary_fidelity,
            "full_spectrum": secondary_full,
            "auxiliary": secondary_auxiliary,
        },
        "strict_activation_example": {
            "probabilities": probabilities,
            "weights": delayed_weights,
            "fidelity": delayed_fidelity,
            "full_spectrum": delayed_full,
            "auxiliary": delayed_auxiliary,
            "best": max(delayed_full, delayed_auxiliary),
        },
    }


def two_qubit_concurrence(state: np.ndarray) -> float:
    sigma_y = np.array([[0.0, -1.0j], [1.0j, 0.0]])
    spin_flip = np.kron(sigma_y, sigma_y)
    product = state @ spin_flip @ state.conj() @ spin_flip
    roots = np.sort(
        np.sqrt(np.maximum(np.linalg.eigvals(product).real, 0.0))
    )[::-1]
    return max(0.0, float(roots[0] - roots[1:].sum()))


def run_exact_two_qubit_benchmark() -> list[dict[str, float]]:
    reference_probability = 0.55
    second_probability = 0.65
    angle = 0.5
    reference = np.array(
        [math.sqrt(reference_probability), 0.0, 0.0,
         math.sqrt(1.0 - reference_probability)]
    )
    second = np.array(
        [math.sqrt(second_probability), 0.0, 0.0,
         math.sqrt(1.0 - second_probability)]
    )
    rotation = np.array(
        [[math.cos(angle), -math.sin(angle)],
         [math.sin(angle), math.cos(angle)]]
    )
    second = np.kron(np.eye(2), rotation) @ second
    maximally_entangled = np.array([1.0, 0.0, 0.0, 1.0]) / math.sqrt(2.0)
    rows: list[dict[str, float]] = []
    for mixture_probability in (0.4, 0.6, 0.8, 0.9):
        state = (
            mixture_probability * np.outer(reference, reference)
            + (1.0 - mixture_probability) * np.outer(second, second)
        )
        fidelity = float(np.vdot(reference, state @ reference).real)
        concurrence = two_qubit_concurrence(state)
        exact_tail = (1.0 - math.sqrt(1.0 - concurrence**2)) / 2.0
        analytic = reference_boundary(fidelity, reference_probability)
        max_fidelity = float(
            np.vdot(maximally_entangled, state @ maximally_entangled).real
        )
        max_reference = reference_boundary(max_fidelity, 0.5)
        rows.append(
            {
                "mixture_probability": mixture_probability,
                "reference_fidelity": fidelity,
                "concurrence": concurrence,
                "exact_E2": exact_tail,
                "reference_bound": analytic,
                "relative_tightness": analytic / exact_tail,
                "maximally_entangled_fidelity": max_fidelity,
                "maximally_entangled_bound": max_reference,
            }
        )
    write_csv("normGM-260814-v1-two-qubit-exact.csv", rows)
    return rows


def clopper_pearson_lower(successes: int, trials: int, alpha: float) -> float:
    if successes == 0:
        return 0.0
    return float(beta.ppf(alpha, successes, trials - successes + 1))


def run_finite_sample_analysis() -> dict[str, object]:
    successes, trials, alpha = 410, 500, 0.01
    estimate = successes / trials
    hoeffding_lower = max(
        0.0, estimate - math.sqrt(math.log(1.0 / alpha) / (2.0 * trials))
    )
    exact_lower = clopper_pearson_lower(successes, trials, alpha)
    threshold = 0.75
    true_fidelity = 0.82
    hoeffding_power = sum(
        float(binom.pmf(count, trials, true_fidelity))
        for count in range(trials + 1)
        if count / trials
        - math.sqrt(math.log(1.0 / alpha) / (2.0 * trials))
        > threshold
    )
    exact_power = sum(
        float(binom.pmf(count, trials, true_fidelity))
        for count in range(trials + 1)
        if clopper_pearson_lower(count, trials, alpha) > threshold
    )
    margin_counts = {
        str(margin): math.ceil(
            (
                math.sqrt(math.log(1.0 / alpha))
                + math.sqrt(math.log(1.0 / 0.1))
            ) ** 2
            / (2.0 * margin**2)
        )
        for margin in (0.07, 0.05, 0.02)
    }
    spectrum_counts: list[dict[str, float | int]] = []
    for dimension, cut in ((8, 4), (32, 16), (1021, 510)):
        log_combinations = (
            math.lgamma(dimension + 1)
            - math.lgamma(cut + 1)
            - math.lgamma(dimension - cut + 1)
        )
        spectrum_counts.append(
            {
                "dimension": dimension,
                "cut": cut,
                "count_delta_0.05": math.ceil(
                    (log_combinations + math.log(1.0 / alpha))
                    / (2.0 * 0.05**2)
                ),
                "count_delta_0.01": math.ceil(
                    (log_combinations + math.log(1.0 / alpha))
                    / (2.0 * 0.01**2)
                ),
            }
        )
    write_csv("normGM-260814-v1-spectrum-sample-complexity.csv", spectrum_counts)
    return {
        "successes": successes,
        "trials": trials,
        "alpha": alpha,
        "estimate": estimate,
        "hoeffding_lower": hoeffding_lower,
        "clopper_pearson_lower": exact_lower,
        "hoeffding_E4_bound": reference_boundary(hoeffding_lower, threshold),
        "clopper_pearson_E4_bound": reference_boundary(exact_lower, threshold),
        "true_fidelity_for_power": true_fidelity,
        "hoeffding_power": hoeffding_power,
        "clopper_pearson_power": exact_power,
        "hoeffding_false_negative": 1.0 - hoeffding_power,
        "clopper_pearson_false_negative": 1.0 - exact_power,
        "hoeffding_sufficient_counts_by_margin": margin_counts,
        "spectrum_calibration_counts": spectrum_counts,
    }


def measurement_resource_example() -> dict[str, float | int | list[float]]:
    probabilities = [0.4, 0.3, 0.2, 0.1]
    rank = len(probabilities)
    shots_per_setting = 2000
    settings = rank * rank - rank + 1
    coefficient_square_sum = probabilities[0] ** 2 + 2.0 * sum(
        probabilities[first] * probabilities[second]
        for first in range(rank)
        for second in range(first + 1, rank)
    )
    failure_probability = 0.01
    radius = math.sqrt(
        2.0
        * math.log(1.0 / failure_probability)
        * coefficient_square_sum
        / shots_per_setting
    )
    return {
        "probabilities": probabilities,
        "settings": settings,
        "shots_per_setting": shots_per_setting,
        "total_shots": settings * shots_per_setting,
        "maximum_basis_switches": settings - 1,
        "failure_probability": failure_probability,
        "confidence": 1.0 - failure_probability,
        "coefficient_square_sum": coefficient_square_sum,
        "hoeffding_radius": radius,
        "diagonal_joint_outcome_bins": rank * rank,
        "off_diagonal_aggregated_joint_outcome_bins": 9,
    }


def main() -> None:
    pure_rows = run_pure_checks()
    noise_rows, dimension_rows = run_mixed_checks()
    weighted = run_weighted_comparison()
    two_qubit_rows = run_exact_two_qubit_benchmark()
    finite_sample = run_finite_sample_analysis()
    summary = {
        "software": {
            "python": sys.version.split()[0],
            "numpy": np.__version__,
            "scipy": scipy.__version__,
            "seed": SEED,
            "pure_optimizer": "SLSQP on general complex coefficient matrices",
            "mixed_optimizer": "rank-rho HJW ensembles with complex Givens search",
            "pure_starts": 10,
            "mixed_starts": 8,
            "mixed_proposals_per_start": 2500,
        },
        "pure_checks": pure_rows,
        "colored_noise": noise_rows,
        "dimension_scan": dimension_rows,
        "weighted_comparison": weighted,
        "exact_two_qubit_benchmark": two_qubit_rows,
        "finite_sample_analysis": finite_sample,
        "measurement_resource_example": measurement_resource_example(),
    }
    DATA_ROOT.mkdir(parents=True, exist_ok=True)
    with (DATA_ROOT / "normGM-260814-v1-numerics.json").open(
        "w", encoding="utf-8"
    ) as handle:
        json.dump(summary, handle, indent=2)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
