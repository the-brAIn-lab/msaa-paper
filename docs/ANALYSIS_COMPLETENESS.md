# Analysis completeness

This cleaned repository preserves the full supplied **code** lineage from the two
HTFA MAT files onward. It separates what is newly automated from what is retained
as the original scientific notebook logic.

| Stage | Included? | Primary location | Notes |
|---|---|---|---|
| HTFA MAT input loading | Yes | `src/msaa_pipeline/data.py` | Starts from the two user-designated MAT files. |
| Spatial MS-AA fitting | Yes | `src/msaa_pipeline/fit.py`, `analysis/01_fit/` | New runner mirrors original settings; original notebook retained. |
| Temporal MS-AA fitting | Yes | same | Uses the same transpose/orientation convention as the original flexible fit notebook. |
| Across-condition fitting | Yes | same | Automated in smoke/mini/full profiles. |
| Within-condition fitting | Yes | same | Automated in mini/full profiles. |
| Full reconstruction decoding | Yes | `src/msaa_pipeline/decode.py` | Resumable CSV checkpoints. |
| Per-archetype decoding/ranking | Yes | same | Condition-specific ranking preserved. |
| Top-m decoding | Yes | same | Automated and resumable. |
| Across/within postprocessing | Yes | `analysis/02_decode/` | Original analysis notebooks retained. |
| Held-out top-m cross-fitting | Yes | `analysis/03_core_results/revision_04_heldout_topm_msaa_crossfit.ipynb` | Paper-specific notebook preserved. |
| PCA / matched-ratio benchmark | Yes | `analysis/03_core_results/revision_05_matched_pca_aa_equivalence_integrated_benchmark.ipynb` | Preserved. |
| Information-concentration / geometry analyses | Yes | `analysis/03_core_results/` | Preserved. |
| K=46 interpretation/support | Yes | `analysis/04_interpretation/` | Preserved. |
| Edge-trim sensitivity | Yes | `analysis/05_robustness/` | Preserved. |
| Network localization/nulls | Yes | `analysis/05_robustness/` | Preserved; standard Schaefer/Yeo atlas can be fetched with `make references`. |
| Four-condition paragraph validation | Yes | `analysis/06_four_condition/` | Fit, decode, and figure notebooks preserved. |
| Simulation benchmarks | Yes | `analysis/07_simulations/` | All current simulation notebooks from the supplied code snapshot retained. |
| Main paper Figures 2–6 | Yes | `paper/main_figures/` | Curated copies preserved. |
| Supplementary analyses/figures | Yes | `paper/supplement/` | Current supplement notebooks preserved. |
| Every uploaded notebook version | Yes | `legacy/original_notebooks/` | Full code-only snapshot preserved for provenance. |

## Important limitation

The uploaded code-only ZIP intentionally excluded the multi-GB fitted NPZs and
most intermediate CSVs. They are therefore not embedded in this cleaned ZIP.
The repository contains the code/configuration to regenerate the core fits and
decoding from the two MAT inputs, while the later paper-specific notebooks
regenerate their own derived products from those core outputs.

The `mini-paper` profile is the recommended migration gate: do not delete the old
large output directories until selected reconstructed fits and decoding outputs
match the original repository.
