# Cleanup decisions

## Preserved

- All 59 notebooks in the supplied code-only snapshot are retained under
  `legacy/original_notebooks/`.
- The original `helpers/` MS-AA implementation is retained unchanged.
- The local `hypertools/` normalization code is retained because the fitting
  notebook used it and exact preprocessing parity matters.
- Current paper/revision notebooks are copied into role-based folders.
- Original repository file/folder inventories are retained in `docs/`.

## Not bundled

The following were intentionally not present in the code-only upload and are not
reconstructed or fabricated here:

- the two required MAT starting inputs;
- multi-GB NPZ model sweeps;
- large raw/derived empirical data trees;
- generated paper figures and most intermediate CSVs.

The clean repository is designed to regenerate those outputs from the two MAT
starting inputs. Small historical outputs can optionally be copied from the old
repo for direct verification.

## Changed in curated notebook copies

Two mechanical path cleanups were applied only to curated copies:

1. `data/pieman/raw/pieman_data.mat` -> `data/input/pieman_data.mat`
2. `data/pieman/raw/pieman_posterior_K700.mat` -> `data/input/pieman_posterior_K700.mat`
3. hard-coded `/Users/lowen/Desktop/papers/archetypes/figures` output paths ->
   `outputs/notebook_figures`

The code-only source snapshot in `legacy/original_notebooks/` remains unpatched
for provenance.

## New code

`src/msaa_pipeline/` and `scripts/run_pipeline.py` are a cleaned reusable version
of the empirical fitting/decoding logic, written to support isolated smoke and
mini-paper validation before committing to a full sweep.

## Standard atlas dependency

Curated notebooks now read the Schaefer 1000-parcel/Yeo-7 files from `data/reference/`. `scripts/fetch_reference_atlas.py` fetches the standard atlas through Nilearn and generates the six-column label table format expected by the historical code.
