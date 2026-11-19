# Data boundary

This repository intentionally starts from two precomputed HTFA inputs:

- `data/input/pieman_data.mat`
- `data/input/pieman_posterior_K700.mat`

The raw-fMRI-to-HTFA preprocessing step is outside the scope of this repository.
The two MAT files are large and are excluded from Git by default.

Place both files in `data/input/` before running any fitting workflow.

## Getting the inputs

The HTFA representation of the Pieman dataset is available as `pieman2_htfa.pkl`
from the Zenodo deposit by Owen & Manning, DOI
[10.5281/zenodo.24358585](https://doi.org/10.5281/zenodo.24358585), licensed
CC BY 4.0. The original fMRI data are from Simony et al. (2016).

Secondary mirror: <https://www.dropbox.com/scl/fi/vw5q3iiej9q288mng4s8w/pieman2_htfa.pkl>

The file holds the same per-scan 700-node HTFA weights and node centers/widths as the
two MAT files; we verified the values are identical. Its paragraph-scrambled
condition has 17 scans because a duplicated entry was removed.

Create the two inputs with:

```bash
uv run python scripts/pkl_to_mat.py <path-to-pieman2_htfa.pkl>
```

This writes `pieman_data.mat` and `pieman_posterior_K700.mat` to `data/input/`.
To match the original file layout, the converter keeps the paragraph-scrambled
condition at 18 columns with the duplicate in column 3. The four-condition
analysis drops that column, so it uses the same 17 scans.

Run:

```bash
uv run make check
```

to confirm that the expected conditions and posterior fields can be read.

The original analysis trims the rest condition by 100 time points and then
normalizes each condition across participants before MS-AA fitting. The cleaned
pipeline preserves that behavior.

## Standard reference atlas

Figure 5 and the node-to-network table use the standard Schaefer 1000-parcel, 7-network atlas. Fetch it into `data/reference/` with:

```bash
uv run make atlas
```

This atlas is treated as an external reference resource rather than a scientific starting dataset.
