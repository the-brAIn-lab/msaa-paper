# Data boundary

This repository intentionally starts from two precomputed HTFA inputs:

- `data/input/pieman_data.mat`
- `data/input/pieman_posterior_K700.mat`

The raw-fMRI-to-HTFA preprocessing step is outside the scope of this repository.
The two MAT files are large and are excluded from Git by default.

Place both files in `data/input/` before running any fitting workflow.

Run:

```bash
python scripts/check_inputs.py
```

to confirm that the expected conditions and posterior fields can be read.

The original analysis trims the rest condition by 100 time points and then
normalizes each condition across participants before MS-AA fitting. The cleaned
pipeline preserves that behavior.

## Standard reference atlas

Later network-localization/brain-map analyses use the standard Schaefer 1000-parcel, 7-network atlas. It can be fetched into `data/reference/` with:

```bash
make references
```

This atlas is treated as an external reference resource rather than a scientific starting dataset.
