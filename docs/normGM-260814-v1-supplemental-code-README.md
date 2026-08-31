# Reproducibility files for `normGM-260814-v1`

## Scope

This archive reproduces the numerical quantities and plot data reported in the
manuscript. It does not contain setting-resolved experimental counts. The
published-data calculations are deterministic post-processing of quantities
reported in the cited articles, and the joint-confidence files are input
templates rather than experimental observations.

## Main calculation

Run from the archive directory:

```text
python normGM-260814-v1-numerics.py
```

The script requires Python 3.10 or later, NumPy, and SciPy. The calculations
were checked with Python 3.13.5, NumPy 2.3.2, and SciPy 1.16.1. The fixed random
seed is `20260806`.

The script writes:

- `normGM-260814-v1-numerics.json`: complete machine-readable summary;
- `normGM-260814-v1-pure-check.csv`: direct constrained pure-state searches;
- `normGM-260814-v1-two-qubit-exact.csv`: exact two-qubit mixed-state benchmark;
- `normGM-260814-v1-colored-noise.csv`: colored-noise lower/upper brackets;
- `normGM-260814-v1-dimension-scan.csv`: dimension dependence;
- `normGM-260814-v1-weighted-comparison.csv`: full-spectrum and largest-coefficient baselines;
- `normGM-260814-v1-spectrum-sample-complexity.csv`: distribution-free calibration counts.

The direct pure-state searches use ten SLSQP starts. The colored-noise upper
estimates use independent Hughston--Jozsa--Wootters ensemble searches and do not
call the manuscript's fixed-fidelity envelope or its KKT equation.

## Published-summary and confidence interfaces

- `normGM-260812-v1-self-testing.py` performs the deterministic screening
  calculation based on digitized published target spectra and reported
  tomography summaries.
- `normGM-260812-v2-experimental-interfaces.py` evaluates fixed input records
  for the experimental interfaces discussed in the manuscript.
- `normGM-260812-v3-joint-confidence.py` implements prospective held-out and
  joint-confidence calculations. Its JSON and CSV files are schemas with
  placeholders, not measured counts.

The corresponding README files document the required fields and assumptions.
The two `.dat` files reproduce the curves in Figs. 1 and 2.

## Interpretation limits

The pure-state optimization, the exact two-qubit family, and the analytic
isotropic calculation provide independent checks within their stated domains.
The higher-dimensional colored-noise ensemble search supplies an upper bound,
not the exact convex roof. No claim of strict mixed-state saturation, a
full-rank saturating family, an independent SDP lower bound, or an end-to-end
raw-count experimental demonstration is encoded in these files.

## Manuscript build

Compile `normGM-260814-v1.tex` with `normGM-260814-v1.bib` using a standard
REVTeX 4.2 installation. The two plot-data files must remain in the same
directory as the TeX source.
