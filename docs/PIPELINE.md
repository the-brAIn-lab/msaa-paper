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
