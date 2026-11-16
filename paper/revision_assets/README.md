# MS-AA paper revision figures and source tables

Source figures and machine-readable tables for the revised manuscript:

**Sparse Distributed Archetypal Motifs Capture Shared Network Structure in Naturalistic fMRI**

Updated: November 2026

## Contents

- `main-pdf/` — vector PDFs used to assemble the main manuscript figures
- `main-png/` — PNG versions of the same figure panels
- `supplement-pdf/` — supplementary figure PDFs
- `supplement-png/` — supplementary figure PNGs
- `tables/` — machine-readable source data underlying the figures and robustness analyses

All held-out analyses use the `1e-8` tie rule. For 10-seed results, the two cross-fit directions are averaged within each outer split, the 20 outer splits are then averaged within each fit seed, and results are finally summarized across the 10 fit seeds.

## Main-manuscript figure mapping

The filenames retain the analysis-pipeline numbering used during revision. Their final manuscript locations are:

| Source filename prefix | Final manuscript figure |
|---|---|
| `Figure2_PCA_vs_MSAA_benchmark` | Figure 3 — empirical PCA / MS-AA / raw-HTFA benchmark |
| `Figure3_simulation_motif_vs_span` | Figure 4 — simulation motif recovery vs. span preservation |
| `Figure4A_spatial_fixed_topm` | Figure 5A — spatial held-out top-m |
| `Figure4B_temporal_fixed_topm` | Figure 5B — temporal held-out top-m |
| `Figure4C_spatial_minus_temporal_matched_ratio` | Figure 5C — matched spatial-minus-temporal comparison |
| `Figure4D_temporal_matched_retained_fraction` | Figure 5D — temporal retained-fraction analysis |
| `Figure5_selected_spatial_motifs_M0-M4` | Figure 6 — representative spatial motifs M0–M4 |
| `Figure6_motif_support_and_robustness` | Figure 7 — motif support and robustness |

## Supplementary figures

- `S01` — matched spatial-versus-temporal comparison, all conditions
- `S02` — temporal K=14 accumulation, generator-anchor masking, and higher-K correspondence
- `S03` — temporal K=12/14/16 seed contrasts
- `S04` — M0–M4 cross-seed stability
- `S05a`, `S05b` — source panels for the combined Supplementary Figure S5 (top-5 and top-10 stability)
- `S06` — M0–M4 selected/recovered across fit seeds
- `S07` — atlas-mapping audit
- `S08` — condition contributions and leave-one-condition-out stability
- `S09` — full-model held-out accuracy vs. model order
- `S10` — edge-trim diagnostics
- `S11` — four-condition paragraph-scrambled sensitivity analysis
- `S12` — complete participant-held-out spatial K x m landscape
- `S13` — simulation sensitivity analyses

## Key reproducibility tables

- `M0-M4_raw_ID_crosswalk.csv` — mapping from manuscript labels M0–M4 to raw K=46, seed-0 component IDs and legacy labels
- `K46_M0-M4_supplementary_table.csv` — recovery and selection summaries for M0–M4
- `M0-M4_cross_seed_stability_summary.csv` — cross-seed spatial, temporal, and Yeo-7 correspondence
- `M0-M4_selection_by_K.csv` — representative-motif selection across neighboring model orders
- `Figure6A_cross_seed_map_r_source.csv` — source data for cross-seed map recovery
- `Figure6B_C_support_source.csv` — source data for motif support/deployment
- `Figure6D_S08_LOCO_map_r_source.csv` — leave-one-condition-out source data
- `Figure5_temporal_profiles.csv` — temporal profiles for representative motifs
- `Figure5_Yeo7_allnode_mass.csv` — Yeo-7 composition for representative motifs
- `spatial_vs_temporal_matched_ratio_summary.csv` / `spatial_vs_temporal_matched_ratio_per_seed.csv` — matched-ratio analyses
- `temporal_matched_fraction_seed_level.csv` / `temporal_K14_local_contrast_delta14.csv` — temporal retained-fraction analyses
- `temporal_K14_anchor_masked_vs_unmasked.csv` — generator-anchor masking
- `S07_atlas_mapping_audit_summary.csv` — atlas-mapping audit summary
- `S09_full_model_accuracy_source.csv` — full-model accuracy curves
- `S10_edge_trim_source.csv` — edge-trim sensitivity
- `S11_*` — four-condition sensitivity source data and legacy comparison
- `S12_spatial_topm_landscape.csv` — full spatial K x m landscape
- `S13_simulation_sensitivity_source.csv` — summarized S13 sensitivity results
- `S13_simulation_sensitivity_runs.csv` — replicate-level S13 simulation output

## Notes

Raw archetype indices are arbitrary across independently initialized fits. Manuscript labels M0–M4 refer to aligned representative motifs; raw component IDs should be used only for reproducibility.

The figure PDFs/PNGs are presentation outputs. The CSV files in `tables/` provide the machine-readable values underlying the plotted analyses.
