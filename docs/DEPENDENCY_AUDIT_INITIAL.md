# MS_AA paper-figure dependency audit

## Main figures

### Figure 2 — empirical PCA / MS-AA benchmark
Inputs:
- `pca_aa_equivalence_outputs/raw_htfa_summary.csv`
- `pca_aa_equivalence_outputs/pca_summary.csv`
- `msaa_condrank_decoding_outputs_spatial_across/full_summary.csv`
- `msaa_condrank_decoding_outputs_temporal_across/full_summary.csv`

### Figure 3 — PCA vs MS-AA simulation benchmark
Input checkpoint:
- `simulation_pca_msaa_extended_outputs/pca_msaa_extended_results_paper.csv`

Code dependency for the representative simulation replicate:
- `helpers/MultiSubject_AA.py` and helper modules

### Figure 4 — held-out sparse MS-AA
Primary input:
- `msaa_heldout_topm_outputs_spatial_across/heldout_topm_crossfit_runs.csv`

Keep this raw cross-fit table rather than only the old summary, because the current figure recomputes the corrected uncertainty after averaging A→B and B→A within outer split.

### Figure 5 — selected K=46 archetypes
The current notebook directly reads:
- `msaa_flexible_outputs_npz/spatialAA_across_acrossCond_K46.npz`
- `data/pieman/raw/pieman_posterior_K700.mat` if geometry is absent from the fit
- Schaefer/Yeo atlas files
- held-out cross-fit results

However, the notebook already exports everything needed to draw the final figure:
- `paper_figure_outputs/Figure5_spatial_maps_allnodes.csv`
- `paper_figure_outputs/Figure5_temporal_profiles.csv`
- `paper_figure_outputs/Figure5_Yeo7_allnode_mass.csv`
- `paper_figure_outputs/Figure5_crossfit_selection_frequency.csv`

Recommended cleanup: refactor Figure 5 into a presentation-only notebook that reads those four source-data tables. This removes the 4.6-GB fit sweep and 223-MB posterior from the figure-reproduction dependency chain.

### Figure 6 — support and robustness
Inputs:
- `revision_05_archetype_interpretation_and_support/selected_archetype_support_by_independent_participant.csv`
- `revision_05_archetype_interpretation_and_support/selected_archetype_map_similarity_by_condition.csv`
- `revision_05_archetype_interpretation_and_support/selected_archetype_condition_summary_independent.csv`
- `msaa_edge_trim_sensitivity_K46_top5/edge_trim_accuracy_summary.csv`
- `msaa_edge_trim_sensitivity_K46_top5/edge_trim_selection_frequency_comparison.csv`

## Supplementary figures

- Supp 1: same inputs as Figure 2.
- Supp 2: `heldout_topm_crossfit_runs.csv`.
- Supp 3: four `full_summary.csv` files for spatial/temporal × across/within.
- Supp 4: spatial-across `topm_summary.csv` + held-out cross-fit table.
- Supp 5: edge-trim summary/Jaccard/frequency/single-component CSVs.
- Supp 6: four-condition raw/PCA/AA retained-information summaries.
- Supp 7: two network-localization CSVs listed below.
- Supp 8: clustering cache `cluster_summary_df.csv` + spatial-across per-archetype decoding CSV.
- Supps 9–11: simulation paper checkpoint.

## Reproducibility issue found

`013_network_localization_analysis_strict_null_matched_ratio.ipynb` sets:

`FIG_ROOT = "/Users/lowen/Desktop/papers/archetypes/figures"`

and exports:
- `strict_random_archetype_null_across.csv`
- `matched_ratio_spatial_temporal_summary_across.csv`

there, outside the repository. The final Supplement 7 PNG/PDF is present in `paper_figure_outputs`, but the source CSVs are absent from the repository file inventory. These should be copied into a repo-local derived-data directory and the notebook should be changed to use a relative path.

## Large directories: recommendation for now

### Preserve/archive, but do not include in the figure-reproduction package
- `msaa_flexible_outputs_npz/` (~4.6 GB): fitted model sweep. Not needed for most final figures; K=46 is currently used by Figure 5 but can be replaced by exported source data at the figure layer.
- `msaa_four_condition_ratio_within_outputs_npz/` (~1.5 GB): upstream fitted objects for the four-condition sensitivity analysis.
- `data/` (~1.1 GB): raw/source data and atlases. Keep separately as analysis inputs; not necessary for a compact paper-figure package after Figure 5 is refactored.

### Likely cleanup/archive candidates after verification
- `notebook_attempts/`
- `.ipynb_checkpoints/`
- `__pycache__/`
- old `*_old` output directories
- redundant historical plotting directories and exported PNG/PDF copies

Do not delete these until the cleaned figure notebooks have been run successfully against the compact source-data bundle.
