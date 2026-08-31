# Reproducibility record

## Final local rerun

- Date: 2026-08-29
- Python: 3.12.13
- NumPy: 2.3.5
- SciPy: 1.18.1
- Fixed random seed: 20260806

Commands executed from `code/`:

```text
python normGM-260814-v1-numerics.py
python normGM-260812-v1-self-testing.py --output ../data/normGM-260812-v1-self-testing-output.json
python normGM-260812-v2-experimental-interfaces.py --output ../data/normGM-260812-v2-experimental-interfaces-output.json
python normGM-260812-v3-joint-confidence.py --self-test
```

All commands completed successfully. The joint-confidence program reported
`self-test passed`. Every generated data file remained byte-for-byte identical
to the archived release candidate, and the data entries in `SHA256SUMS.txt`
remained valid after the rerun.

The self-testing output stores full floating-point computational values, while
the manuscript reports the digitization-based probabilities and bounds to
three decimal places. Neither representation is a confidence region. No raw
experimental counts were used in these reruns.
