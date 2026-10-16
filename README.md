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

## Portable reuse of the validated original full sweep

The `mini-paper` rerun is the migration/reproducibility gate. Once representative
fits and decoding reproduce the original analysis, the recommended manuscript
workflow is to **reuse the already-computed expensive full K-grid outputs** rather
than spend days recomputing them solely for migration.

For a portable project folder (including student handoff), copy the legacy
derived outputs into the ignored local store:

```bash
make import-legacy OLD_REPO=/path/to/original/MS_AA
make validate-legacy
make activate-legacy
```

The imported bytes live under `reference_outputs/legacy_full/`, which is ignored
by Git. `configs/legacy_outputs_manifest.csv` is small and should be committed;
it records portable relative paths, sizes, and SHA-256 hashes. `make
activate-legacy` creates relative links at the historical output paths expected
by older notebooks.

A Git clone alone does **not** include ignored legacy data or the two large MAT
inputs. To give another person the complete working project, copy/archive the
whole project directory rather than relying only on Git, then have them run:

```bash
make handoff-check
make activate-legacy
```

See `reference_outputs/README.md` and
`docs/REPRODUCIBILITY_VALIDATION.md`.

## Optional full core sweep

A complete from-scratch rebuild remains available when desired. Deactivate the
legacy links first so regenerated outputs cannot overwrite or mingle with the
archived validated results:

```bash
make deactivate-legacy
make full-fit     # fit only; expensive
make full-core    # fit + decoding; very expensive
```

The full targets refuse to run while validated legacy links are active. The
archived bytes remain safely stored under `reference_outputs/legacy_full/`.

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

Do not commit `data/input/*.mat` or `reference_outputs/legacy_full/`. They are
intentionally local because they are large. The tracked input checksums, frozen
analysis configuration, and legacy-output SHA-256 manifest provide the audit
trail; the complete local project folder provides the convenient student
handoff.
