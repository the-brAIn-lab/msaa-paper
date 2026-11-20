# Reproducing the paper

Figure numbers in this repository (file names, docs, build steps) follow the figure-bundle numbering: Fig 2 = main Fig 3, Fig 3 = main Fig 4, Fig 4A/B = main Fig 5A/B, Fig 4C = main Fig 6A, Fig 4D = main Fig 6B, Fig 5 = main Fig 7, Fig 6/Fig 8 = main Fig 8. Supplementary numbers are unchanged.

Everything starts from the two HTFA inputs in `data/input/` (see `data/README.md`). The chain has
two parts: **results** (fits, decoding, held-out top-m and the follow-up analyses; expensive) and
**figures** (`scripts/build_figures.py`; minutes, reads only result tables).

Run every command from the repo root with the project environment (`uv sync`, then `uv run ...`).
Results land in the repo root (gitignored), the "run root". To draw figures from results kept
somewhere else, pass that folder as `--root` (or set `MSAA_ROOT`); `build_figures.py` links the
code folders into it.

All decoding uses the 1e-8 tie rule (`MSAA_TIE_TOL`, default 1e-8). 10-seed results average the two
directions within each outer split, then the 20 splits within a seed, then the 10 fit seeds.

## 1. Results

| run root entry                                                                                                                       | made by                                                                                                                                                        | used by                          |
| ------------------------------------------------------------------------------------------------------------------------------------ | -------------------------------------------------------------------------------------------------------------------------------------------------------------- | -------------------------------- |
| `msaa_flexible_outputs_npz/`, `msaa_condrank_decoding_outputs_<analysis>_across/`, `msaa_heldout_topm_outputs_spatial_across/` | `analysis/08_multiseed_stability/run_local.py --steps seed0` (seed-0 fits with the frozen picks in `configs/`, decoding, held-out top-m at the 30 paper K) | Fig 2, 5, everything downstream |
| `seeds/`                                                                                                                           | `run_local.py --steps picks,fits,decode,merge` (seeds 1-9)                                                                                                   | S04-S06 tables                   |
| `heldout-everyk-<analysis>-<seed0\|all>-results/`                                                                                   | `run_local.py --steps heldoutK` (held-out top-m at every K, 10 seeds)                                                                                        | Fig 4A/B, 4C, S01, S09, S12      |
| `heldout-highm-results/`                                                                                                           | `run_local.py --steps highm` (m = 7-100 at the 15 matched spatial/temporal K pairs)                                                                          | Fig 4C, S01                      |
| `heldout-highm-landscape-results/`                                                                                                 | `run_local.py --steps landscape` (m = 7-50, spatial K = 10-700)                                                                                              | S12                              |
| `outputs/multiseed_stability/matched_fraction/aggregate/`                                                                          | `analysis/08_multiseed_stability/README.md`, section 3b                                                                                                      | Fig 4D, S03                      |
| `outputs/multiseed_stability/k50_cross_seed_stability/`                                                                            | `plots/loadx.py`, `plots/fig5_per_seed.py 50 28,39,16,37,49`, `k50_stability/k50_cross_seed_stability.py`                                                | Fig 6A, S04                      |
| `outputs/multiseed_stability/topm_stability_neighbourhood/`                                                                        | `plots/topm_stability_neighbourhood.py 45,...,55`                                                                                                            | S05, S06, tables                 |
| `pca_aa_equivalence_outputs/`                                                                                                      | `analysis/03_core_results/revision_05_matched_pca_aa_equivalence_integrated_benchmark.ipynb`                                                                 | Fig 2                            |
| `outputs/temporal_compact_dictionary_K14_local/`                                                                                   | `analysis/03_core_results/revision_10a_*.ipynb`, then `revision_10c_*.ipynb`                                                                               | S02                              |
| `revision_05_archetype_interpretation_and_support/`                                                                                | `analysis/04_interpretation/revision_06_*.ipynb`                                                                                                             | Fig 6B/C                         |
| `outputs/interpretation/condition_influence_and_seed_stability_K50/`                                                               | `analysis/04_interpretation/revision_09a_*.ipynb` (leave-one-condition-out refits), then `revision_09_*.ipynb`                                             | Fig 6D, S08                      |
| `msaa_edge_trim_sensitivity_K50_top5/`                                                                                             | `analysis/05_robustness/revision_07_*.ipynb`                                                                                                                 | S10                              |
| `four_condition_ratio_decoding_outputs/`                                                                                           | `analysis/06_four_condition/rerun_four_condition.py all --out .`                                                                                             | S11                              |
| `four_condition_comparison/`                                                                                                       | `analysis/06_four_condition/compare_four_condition.py` (against the August 2026 tables)                                                                      | S11 table                        |
| `simulation_pca_msaa_extended_outputs/`                                                                                            | `analysis/07_simulations/revision_simulation.ipynb`, then `s13_sensitivity_sweep.py`                                                                       | Fig 3, S13                       |
| `data/reference/`                                                                                                                  | `make atlas`                                                                                                                                                 | Fig 5, S07                       |

Run a notebook with `uv run python scripts/execute_notebook.py <notebook>` (it runs from the repo
root). `uv run make results` runs every `run_local.py` step; that is weeks of CPU time on one machine
(the paper used about 100 cluster cores for several days), so `run_local.py --ks/--seeds/--jobs`
restricts it for testing. Every step is resumable, and the high-m steps check themselves against
the every-K run (`VALIDATION PASS` in their `summary.txt`).

**Held-out top-m** is always participant-held-out: archetypes are ranked on one half of the
participants and the top m are scored on the other half (20 outer splits x 2 directions).

## 2. Figures

```bash
uv run python scripts/build_figures.py --root <run root> --out <bundle folder>
uv run python scripts/build_figures.py --only fig5,s07          # some figures
uv run python scripts/build_figures.py --no-draw                 # just reassemble the bundle
```

Each step draws one figure (or a few) from the result tables and writes it under
`<run root>/outputs/`; the bundle gets PNG + PDF of every figure, the source tables in `tables/`,
and `paper/FIGURE_BUNDLE_README.txt` as its README. Figure captions are in the paper.
Every figure is drawn at the journal's print width (about 5 in), so the type sizes are the printed
sizes.

| figure          | drawn by                                                                                                 |
| --------------- | -------------------------------------------------------------------------------------------------------- |
| Fig 2           | `paper/main_figures/paper_figure02_*.ipynb` (source table), `figure02_plot.py --sharey`              |
| Fig 3           | `paper/main_figures/paper_figure03_*.ipynb`                                                            |
| Fig 4A/B/D, S09 | `analysis/08_multiseed_stability/plots/fig4_heldout_panels.py`                                         |
| Fig 4C, S01     | `plots/fig4c_matched_ratio_heatmaps.py`                                                                |
| Fig 5           | `paper/main_figures/paper_figure05_*.ipynb`                                                            |
| Fig 6           | `paper/main_figures/paper_figure06_*.ipynb`                                                            |
| Fig 8 (main, A-F) | `paper/main_figures/figure08_support_robustness.py` (build step `fig8`)                              |
| S02             | `plots/k14_supplement.py`                                                                              |
| S03             | `matched_fraction/k14_local_contrast.py`                                                               |
| S04             | `k50_stability/m0m4_summary.py`                                                                        |
| S05a/b, S06     | `plots/topm_stability_neighbourhood.py --plot-only`; tables from `k50_stability/make_supp_tables.py` |
| S07, S08        | `paper/supplement/supp_atlas_condition_plot.py` (after Fig 5)                                          |
| S10             | `paper/supplement/supp_05_edge_trim_diagnostics.ipynb`                                                 |
| S11             | `paper/supplement/supp_four_condition_plot.py`                                                         |
| S12             | `plots/spatial_topm_landscape.py`                                                                      |
| S13             | `paper/supplement/supp_13_simulation_sensitivity.py`                                                   |

Figure 1 is a schematic and is not generated here.
