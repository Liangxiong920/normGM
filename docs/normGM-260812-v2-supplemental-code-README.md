# Supplemental numerical code for normGM-260812-v2

The files in this bundle reproduce the numerical quantities reported in the
analytical and experimental-interface examples. They use only the Python
standard library and require Python 3.10 or later.

## Files

- `normGM-260721-v1-numerics.py`: analytical benchmarks, fixed-fidelity
  boundaries, multiscale profiles, and figure data.
- `normGM-260812-v1-self-testing.py`: deterministic screening calculation from
  the digitized target spectra and tomography fidelities reported by Zhang et
  al.
- `normGM-260812-v2-experimental-interfaces.py`: deterministic conversion of
  published nonuniform-projector, 32-dimensional spatial-mode, and
  1021-dimensional time-frequency fidelity inputs into Vidal-tail bounds.
- `normGM-260721-v1-cross-isotropic.dat` and
  `normGM-260721-v1-reference-scan.dat`: data plotted in the manuscript.
- `normGM-260812-v2-experimental-interfaces-output.json`: unrounded output from
  the experimental-interface script.

## Run

```text
python normGM-260721-v1-numerics.py
python normGM-260812-v1-self-testing.py
python normGM-260812-v2-experimental-interfaces.py
```

## Experimental input provenance

The experimental-interface script records the source DOI for every numerical
case. The nonuniform four-dimensional target and reported fidelity are from
Guo et al., Phys. Rev. A 97, 062309 (2018). The 32-dimensional fidelity is from
Hu et al., Phys. Rev. Lett. 125, 090503 (2020). The 1021-dimensional input is
from Chang et al., Sci. Adv. 12, eaee1333 (2026).

The values 0.988 and 0.932 are obtained by subtracting the displayed error from
the corresponding central value. They are used only as conservative screening
inputs and are not assigned a confidence interpretation. For the
1021-dimensional case, 0.650 is the lower endpoint of the source article's
displayed three-standard-deviation error bar and retains that article's
Poisson/Gaussian resampling assumptions.

## Scope

The scripts evaluate closed-form formulas using reported inputs. They do not
reanalyze raw detector counts, reconstruct Bell or local correlations, infer a
new fidelity confidence interval, propagate target-spectrum uncertainty, or
constitute an independent experimental validation of the analytical theorem.
