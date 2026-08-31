# Reproducibility material for `normGM`

This repository accompanies the manuscript *Multiscale Schmidt-Spectrum Bounds
for High-Dimensional Entanglement: Geometric Measures and Schmidt-Number
Witnesses*. It contains the scripts, machine-readable numerical outputs, and
statistical input templates referenced in the numerical and
experimental-interface section.

## Contents

- `code/`: numerical benchmarks and published-summary/statistical interfaces;
- `data/`: generated CSV, JSON, and plot-data files;
- `templates/`: empty schemas for setting-resolved counts and calibration data;
- `docs/`: detailed assumptions and run instructions.

The repository contains no setting-resolved experimental counts. The
published-summary calculations post-process quantities reported in the cited
articles. The raw-count templates contain placeholders and do not constitute
experimental observations. Integer bar heights digitized from a 300-dpi figure
rendering are included as deterministic screening inputs. They have no
calibrated pixel-to-amplitude uncertainty model and do not define a confidence
region; manuscript values derived from them are rounded to three decimal
places.

## Reproduce the numerical files

Python 3.10 or later, NumPy, and SciPy are required. From the repository root,
run

```text
cd code
python normGM-260814-v1-numerics.py
python normGM-260812-v1-self-testing.py --output ../data/normGM-260812-v1-self-testing-output.json
python normGM-260812-v2-experimental-interfaces.py --output ../data/normGM-260812-v2-experimental-interfaces-output.json
python normGM-260812-v3-joint-confidence.py --self-test
```

The main numerical script uses the fixed random seed `20260806`. The archived
outputs in `data/` permit direct comparison with a fresh run. See the files in
`docs/` for input assumptions and interpretation limits.

The final release checks record the Python, NumPy, and SciPy versions, command
outputs, and SHA-256 checksums in `REPRODUCIBILITY_RECORD.md` and
`SHA256SUMS.txt`.

## License

The repository uses the MIT License. `LICENSE_SCOPE.md` explains how the
license applies to authored materials, generated outputs, and numerical inputs
transcribed from third-party publications.

## Statistical scope

The joint-confidence program supports held-out target selection, a fixed but
uncertain implemented spectrum, and post-processing of a user-supplied joint
region. It does not prove that an experimental acquisition satisfies those
models. Basis alignment, phases, leakage, stopping rules, and target selection
must be certified by the experimenter.

## Citation and archival release

Please cite the associated manuscript and the archived `v1.0.0` release. The
software citation metadata are provided in `CITATION.cff`. The immutable
release is archived at <https://doi.org/10.5281/zenodo.22208160>; the concept
DOI for all versions is <https://doi.org/10.5281/zenodo.22208159>.
