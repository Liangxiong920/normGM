# Supplemental numerical material for normGM-260812-v1

This bundle reproduces the deterministic calculations reported in Sec. VI of
the manuscript. Both scripts use only the Python standard library.

## Files

- `normGM-260721-v1-numerics.py`: analytical benchmarks, witness comparisons,
  and data generation for the manuscript figures.
- `normGM-260812-v1-self-testing.py`: retrospective screening calculation for
  the data in Fig. 4 of W.-H. Zhang et al., *npj Quantum Information* 5, 4
  (2019), doi:10.1038/s41534-018-0120-0.
- `normGM-260721-v1-cross-isotropic.dat`: data used by the cross-isotropic
  comparison figure.
- `normGM-260721-v1-reference-scan.dat`: data used by the reference-spectrum
  scan figure.

## Reproduction

Run the scripts from this directory with Python 3.10 or later:

```text
python normGM-260721-v1-numerics.py
python normGM-260812-v1-self-testing.py
```

The first command regenerates the two `.dat` files and prints the benchmark
tables as JSON. The second prints the unrounded entries underlying Table II.

## Scope of the published-data calculation

The tomography fidelities and reported standard deviations are transcribed
from Fig. 4 of the cited experiment. The Schmidt-amplitude inputs are central
bar heights digitized from a 300-dpi rendering of that figure and are stored
explicitly in the script. The calculation does not use the experiment's raw
counts, reconstruct Bell correlations, or infer a high-dimensional
self-testing fidelity bound. In particular, `F_T` minus one reported standard
deviation is only a lower-shifted screening input; it is not a one-sided
confidence limit.

For a confidence-qualified device-independent application, the numerical
self-testing fidelity lower bound (or raw Bell data), the target coefficients,
and a joint statistical coverage statement are still required.
