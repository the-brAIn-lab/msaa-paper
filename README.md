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
analysis/                 scientific analyses, in pipeline order (01-07 notebooks,
                          08_multiseed_stability scripts)
paper/                    final manuscript and supplement plotting notebooks
scripts/                  runners: run_pipeline.py (smoke/mini), run_parallel.py
                          (canonical full run), init freezing, node table, checks
src/msaa_pipeline/        reusable cleaned fit/decode/QC implementation
helpers/                  original MS-AA implementation used by the analyses
hypertools/               legacy normalization utilities required for parity
configs/                  run profiles, frozen FurthestSum picks
                          (furthestsum_init.csv.gz), node -> Yeo-7 table (node_yeo7.csv)
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
python3.12 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

or, with [uv](https://docs.astral.sh/uv/) (uses `uv.lock`, the exact environment of the paper runs):

```bash
uv sync
source .venv/bin/activate
```

All three install the same pinned versions (Python 3.12; numpy 2.5.4, scipy 1.18.1, pandas 3.1.0, ...).

Then place the two MAT files in `data/input/` and run:

```bash
make check
```

Node-localization analyses also use the standard Schaefer-1000/Yeo-7 atlas. Fetch the exact reference layout expected by the historical notebooks with:

```bash
make references
```

This reference atlas is not an additional empirical study input. The Yeo-7 label of
each of the 700 HTFA nodes is committed as `configs/node_yeo7.csv`; it is rebuilt from
the atlas with `python scripts/build_node_yeo7.py` (same rule as the interpretation
notebooks, and it matches the labels behind the paper figures for all 700 nodes).

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

If you have the original analysis outputs, they can be reused instead of recomputing
the full K grid. `mini-paper` checks that the code reproduces them on a few model orders.
Copy the legacy derived outputs into the ignored local store:

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
inputs. A copied project directory that includes them can be checked and linked with:

```bash
make handoff-check
make activate-legacy
```

See `reference_outputs/README.md` and
`docs/REPRODUCIBILITY_VALIDATION.md`.

## Canonical full run (the paper's numbers)

The fits, decoding and held-out top-m behind the revised paper were produced with
`scripts/run_parallel.py`, a parallel driver around the same `msaa_pipeline` functions.
Compared with the original analysis it:

- **fits every K** (profile `full`: spatial 2-700, temporal 2-300, for all four
  analysis/scope jobs; the original, sparser grid is kept as profile `paper_grid`);
- **freezes the FurthestSum starting picks** in `configs/furthestsum_init.csv.gz`
  (`scripts/freeze_init.py`, `scripts/freeze_everyk.py`) and does not re-sort archetypes
  after fitting, so archetype numbers are sequential and identical on every machine.
  Archetype IDs in the revised paper and figures are these 0-based IDs;
- **decodes top-m subsets with an exact low-rank method** (`decode_subjects(...,
  method="auto")`, at most 10 archetypes), about 12x faster and equal to the dense
  decoder up to rounding (`scripts/verify_lowrank.py`);
- scores correlations with a **tie tolerance** of 1e-8 (`--tie-tol`) so that exact ties,
  common in temporal fits, break the same way on every machine.

```bash
python scripts/run_parallel.py fit     --workers 8 [--jobs spatial_across] [--ks 5,46]
python scripts/run_parallel.py decode  --workers 8
python scripts/run_parallel.py heldout --workers 8   # participant-held-out top-m, 30 seeded K
python scripts/run_parallel.py merge
```

Each fit, decode and held-out task writes its own output file, so a run can be stopped and resumed. The full every-K
run is several thousand fits (days on a workstation). Held-out top-m at every K and the
10-seed stability analysis are in `analysis/08_multiseed_stability/` (see its README).

Across machines, spatial results reproduce exactly in numbering and to r >= 0.9999 in
the maps; decoding and held-out summaries agree within about 0.005. A few temporal
archetypes at high K can land in a slightly different local optimum (mostly a 1-3 TR
shift); nothing forces bit-identical BLAS results.

## Optional full core sweep (serial runner)

`make full-fit` / `make full-core` run `scripts/run_pipeline.py --profile full`
serially, on one core. For the full profile it uses the same settings as
`scripts/run_parallel.py`: frozen FurthestSum picks with no post-fit sort (the paper's
archetype numbering) and decode tie tolerance 1e-8. The `smoke` and `mini-paper` profiles
keep the legacy settings (FurthestSum recomputed, archetypes sorted after fitting, exact tie
comparison), because they are checked against the legacy outputs. `--numbering` and
`--tie-tol` override either default. `scripts/run_parallel.py` above is the faster way to
reproduce the paper.

Deactivate any legacy links first so regenerated outputs cannot overwrite or mingle with the
archived validated results:

```bash
make deactivate-legacy
make full-fit     # fit only; expensive
make full-core    # fit + decoding; very expensive
```

The full targets refuse to run while validated legacy links are active. The
archived bytes remain safely stored under `reference_outputs/legacy_full/`.

## Downstream analyses and figures

Participant-held-out top-m cross-fitting is the **only** canonical top-m analysis
(`docs/HELDOUT_TOPM.md`). After the core fits and decoding, the PCA benchmark, K=46
interpretation, robustness tests, four-condition validation, simulations and manuscript
figures are organized notebooks under `analysis/` and `paper/`. From the repository root:

```bash
make full-heldout   # held-out top-m cross-fit only, once fits exist
make paper-main     # execute the main-figure notebooks (needs the upstream result tables)
```

`scripts/execute_notebook.py` runs a notebook with the repo root as working directory, so
the notebooks' repo-relative paths resolve.

The multi-seed stability analysis (10 initialisations per K, every-K held-out top-m and
the stability figures) is plain Python; see `analysis/08_multiseed_stability/README.md`.

See `docs/PIPELINE.md` for the dependency chain and `docs/NOTEBOOK_MAP.md` for
where each analysis now lives.

## Important preservation rule

Do not commit `data/input/*.mat` or `reference_outputs/legacy_full/`. They are
intentionally local because they are large. The tracked input checksums, frozen
analysis configuration, and legacy-output SHA-256 manifest provide the audit
trail.
