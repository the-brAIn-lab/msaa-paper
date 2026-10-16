# Analysis pipeline

## Reproducibility boundary

The reproducible analysis starts from:

1. `pieman_data.mat`
2. `pieman_posterior_K700.mat`

No raw NIfTI preprocessing is required inside this repository.

## Main dependency chain

```text
pieman_data.mat + pieman_posterior_K700.mat
    |
    v
01_fit_msaa
    |  spatial / temporal
    |  across-condition / within-condition
    v
02_decode
    |  full reconstruction
    |  per-archetype decoding
    |  single-archetype diagnostics
    v
held-out top-m cross-fitting
    |  A ranks -> B evaluates; B ranks -> A evaluates
    |  only canonical sparse-selection analysis
    v
03_core_results
    |  matched-ratio / PCA summaries
    |  PCA comparison
    |  information concentration / clustering
    v
04_interpretation
    |  representative K=46 motif support
    |  temporal profiles / maps / Yeo-7 composition
    v
05_robustness
    |  edge trimming
    |  network localization / nulls
    v
06_four_condition
    |  paragraph-scrambled validation
    v
paper/main_figures + paper/supplement
```

Simulation analyses in `analysis/07_simulations/` are largely independent of the
empirical fit chain and provide method validation / robustness controls.

## What is automated in the cleaned runner?

`scripts/run_pipeline.py` automates the core, deterministic chain from the MAT
inputs through MS-AA fitting, full/single-archetype decoding, and the canonical held-out top-m cross-fit. This is the part needed for the fast
`smoke` and `mini-paper` checks.

The later revision analyses remain as curated notebooks because they contain
paper-specific analysis choices (K=46 interpretation,
edge-trim sensitivity, network nulls, and simulation variants). They are now
organized by stage rather than mixed with exploratory attempts.

This is deliberate: the cleanup preserves the scientific code rather than
silently rewriting every paper analysis into a new implementation at once.

## Canonical full model grids

The `full` profile records the K values found in the fitted-model inventory from
the original repository. It therefore documents the actual model-order sweeps
that produced the existing outputs, rather than relying on whichever temporary
K list happened to be left in a notebook cell.

Run the full fit sweep with:

```bash
make full-fit
```

and the full core fit + decoding sweep with:

```bash
make full-core
```

Both are computationally expensive.


See `docs/HELDOUT_TOPM.md` for the canonical sparse-analysis definition. Historical same-data top-m analyses are provenance-only.

## Portable validated legacy full-sweep mode

After `mini-paper` reproduces representative fits/decoding, the cleaned repo can
reuse the original expensive full K-grid outputs as explicitly validated derived
artifacts. They are copied **inside the project folder** but outside Git:

```text
reference_outputs/legacy_full/
```

Set up once with:

```bash
make import-legacy OLD_REPO=/path/to/original/MS_AA
make validate-legacy-hash   # optional full cryptographic verification
make activate-legacy
```

Activation creates relative symlinks at the historical top-level output paths so
curated legacy notebooks do not need machine-specific paths. The links point
only within the same project folder, so they can be recreated after moving or
copying the repository.

The optional `full` profile remains available for a future from-scratch rebuild,
but it is not required merely to establish migration reproducibility. Run `make
deactivate-legacy` before any full write-producing target.

See `docs/REPRODUCIBILITY_VALIDATION.md` for provenance and handoff details.
