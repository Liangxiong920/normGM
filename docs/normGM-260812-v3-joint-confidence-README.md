# Joint one-sided confidence analysis

This supplement records what is needed to turn the manuscript's published-data screening calculations into confidence-qualified certificates. It contains no experimental count data.

## Files

- `normGM-260812-v3-joint-confidence.py` implements Propositions `held-out certificate for a selected reference` and `raw-count joint one-sided certificate` from Section III.
- `normGM-260812-v3-held-out-counts-template.json` is for a target selected from calibration data and frozen before independent fidelity-test rounds.
- `normGM-260812-v3-spectrum-counts-template.json` is for a fixed reference whose Schmidt probabilities are estimated from multinomial mode counts.
- `normGM-260812-v3-setting-counts-template.csv` is an archival layout for setting-resolved detector records. The Python script does not infer a fidelity bound from arbitrary local-setting rows; that conversion remains experiment-specific.

The `null` fields are deliberate. The program rejects a template until actual integer counts and provenance are supplied.

## Statistical distinction

### Selected and frozen target

Use `held_out_selected_reference` when a calibration record selects the normalized target spectrum. After selection, freeze the target and acquire fresh test rounds. Conditional on calibration, the target and every cumulative weight `A_r` are fixed, so one lower confidence bound for the held-out fidelity gives simultaneous Vidal-tail bounds for all requested thresholds. This certifies the declared selected target; it does not identify an unknown source spectrum.

### Fixed but uncertain implemented spectrum

Use `calibrated_unknown_spectrum` only when the reference is fixed independently of the observations. The script treats the mode counts as IID multinomial trials and obtains simultaneous upper bounds on the sums of the largest `r-1` probabilities. It combines these with a one-sided Hoeffding fidelity bound by a union bound. Dependence between the two confidence events is allowed, but each marginal sampling model must be valid.

Mode counts represent Schmidt probabilities only under a validated model for mode labels, basis alignment, phases, and leakage. The script cannot establish that model.

### Overlapping adaptive data

Do not use either pointwise mode if the same observations select a target and evaluate it. Supply an experiment-specific, selection-uniform joint region through `precomputed_joint_region`. The script will post-process such a region but cannot validate its coverage proof.

## Required raw-data record

Keep the unrounded integer counts, setting and outcome labels, trials or exposure times, acquisition batches and timestamps, stopping rule, randomization procedure, accidental/background policy, detector corrections, measurement operators, local-mode ordering, target-selection code, calibration/test split, confidence budget, and dimension/leakage evidence. Preserve a cryptographic hash of every raw file in the JSON provenance fields.

## Run

```text
python normGM-260812-v3-joint-confidence.py INPUT.json --output RESULT.json
python normGM-260812-v3-joint-confidence.py --self-test
```

The implementation uses conservative Hoeffding and union bounds. A sharper exact-binomial, likelihood, confidence-sequence, or martingale analysis may replace them only with a coverage argument matched to the acquisition protocol.
