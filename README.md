# MS-AA archetypes paper analysis

Cleaned analysis repository for the archetypal-analysis manuscript.

## Reproducibility boundary

The repository starts from two precomputed HTFA products:

```text
data/input/pieman_data.mat
data/input/pieman_posterior_K700.mat
```

Generation of those HTFA products from raw neuroimaging data is intentionally
outside this repository. Everything from MS-AA fitting onward is retained here.

## What was cleaned

The original working repository mixed model fits, exploratory notebooks,
revision analyses, manuscript plotting code, large generated objects, and old
attempts at the top level. This version separates them into:

```text
analysis/                 scientific analyses, in pipeline order
paper/                    final manuscript and supplement plotting notebooks
src/msaa_pipeline/        reusable cleaned fit/decode/QC implementation
helpers/                  original MS-AA implementation used by the analyses
hypertools/               legacy normalization utilities required for parity
configs/                  full and small replication profiles
data/input/               the two required MAT inputs (not committed)
data/derived/             optional curated derived/source tables
outputs/                  isolated smoke / mini-paper outputs
reference/                optional original-repo references for comparisons
legacy/original_notebooks complete code-only snapshot of uploaded notebooks
docs/                     provenance, pipeline, and replication notes
```

No scientific notebook from the supplied code snapshot was deleted: the entire
set is retained under `legacy/original_notebooks/`. Curated copies are placed in
`analysis/` and `paper/` according to their role.

## Setup

```bash
conda env create -f environment.yml
conda activate msaa-paper
```

or:

```bash
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

Then place the two MAT files in `data/input/` and run:

```bash
make check
```

Later node-localization analyses also use the standard Schaefer-1000/Yeo-7 atlas. Fetch the exact reference layout expected by the historical notebooks with:

```bash
make references
```

This reference atlas is not an additional empirical study input.

## Fast end-to-end check

```bash
make smoke
```

This deliberately small run uses only 4 participants per condition, K=5, fewer
fit iterations, and fewer decoding repetitions. It verifies that the reorganized
code can traverse:

```text
MAT inputs -> normalization -> MS-AA -> full/single-archetype decoding -> participant-held-out top-m -> QC
```

It writes only to `outputs/smoke/` and cannot overwrite paper results.

## Small exact-replication check

```bash
make mini
```

`mini-paper` uses **all participants and the original fit/decoding settings**, but
only a few model orders. It runs K=5 across spatial/temporal and within/across
families, plus the representative spatial across-condition K=46 fit.

If the original repository is still available, compare directly:

```bash
make compare OLD_REPO=/path/to/original/MS_AA
```

The comparison checks matching decoding rows and first-subject reconstructed
fits. See `docs/QUICK_REPLICATION.md`.

## Full core sweep

The full K grids are recorded in `configs/profiles.json` from the actual fitted
model inventory in the original repository.

```bash
make full-fit     # fit only; expensive
make full-core    # fit + decoding; very expensive
```

The full runner writes canonical output directory names when invoked through
these Make targets so the historical downstream notebooks can find them.

## Later paper analyses

The held-out top-m cross-fit is now automated in `scripts/run_pipeline.py` and is the **only** canonical top-m analysis. After the core fits/decoding, the PCA benchmark,
K=46 interpretation, robustness tests, four-condition validation, simulations,
and manuscript figures are preserved as organized notebooks under `analysis/`
and `paper/`.

These notebooks can be executed from the repository root with:

```bash
make full-heldout
```

The original held-out notebook is retained for provenance, while the reusable runner mirrors its final algorithm. The helper forces the working directory to the repo root so historical
repo-relative paths continue to resolve after the notebooks were reorganized.

See `docs/PIPELINE.md` for the dependency chain and `docs/NOTEBOOK_MAP.md` for
where each analysis now lives.

## Important preservation rule

Do not discard the old large fit/output directories until the mini-paper
comparison succeeds. The clean repo intentionally does not bundle multi-GB NPZ
sweeps; instead it preserves the code/configuration needed to regenerate them
and provides isolated verification modes before migration.
