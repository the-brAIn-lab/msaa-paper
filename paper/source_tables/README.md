# Source tables

These are the exact tables behind the published figures, copied unchanged from the figure bundle of
17 Nov 2026 (`tables/`). Figure and supplementary-figure numbers below follow the manuscript. File
names keep the older bundle numbering, which is one lower than the manuscript for the main figures
(for example `Figure2_*` is main Figure 3).

Archetype indices (the `archetype`, `k46_seed0_id` and similar columns) refer to the spatial MS-AA
fit at K = 46, fit seed 0, and are 0-based. They are arbitrary up to permutation: a different seed
or software version can number the same archetypes differently. Motif labels M0-M4 are the stable
names; `M0-M4_raw_ID_crosswalk.csv` maps them to the raw indices.

| File | Underlies |
| --- | --- |
| `Figure2_source.csv` | Main Fig. 3: empirical benchmark, PCA reconstruction vs. spatial and temporal MS-AA, with raw HTFA reference |
| `Figure4A_source.csv` | Main Fig. 5A: spatial MS-AA held-out top-m accuracy vs. K (mean, SD over 10 fit seeds) |
| `Figure4B_source.csv` | Main Fig. 5B: temporal MS-AA held-out top-m accuracy vs. K |
| `Figure4D_source.csv` | Main Fig. 5D: temporal MS-AA retained fraction at 25/50/70% of K |
| `temporal_matched_fraction_seed_level.csv` | Main Fig. 5D and Supp. Fig. S3: per-seed values behind the retained-fraction curves |
| `temporal_K14_local_contrast_delta14.csv` | Supp. Fig. S3: K = 12/14/16 local contrast (Delta14), per seed with bootstrap CIs; text result on K = 14 |
| `temporal_K14_anchor_masked_vs_unmasked.csv` | Supp. Fig. S2A-B: temporal K = 14 accumulation over m and generator-anchor masking |
| `spatial_vs_temporal_matched_ratio_summary.csv` | Main Fig. 5C and Supp. Fig. S1: spatial minus temporal accuracy at matched ratios, seed means and seed counts |
| `spatial_vs_temporal_matched_ratio_per_seed.csv` | Main Fig. 5C and Supp. Fig. S1: per-seed comparisons (4,410 rows) behind the summary |
| `Figure5_temporal_profiles.csv` | Main Fig. 6: temporal profiles of M0-M4 (mean and SEM over scans) |
| `Figure5_Yeo7_allnode_mass.csv` | Main Fig. 6: Yeo-7 composition of M0-M4, all 700 nodes; text results on network shares |
| `Figure5_crossfit_selection_frequency.csv` | Main Fig. 6: cross-fitted top-5 selection frequency by condition (10-seed mean) |
| `Figure6A_cross_seed_map_r_source.csv` | Main Fig. 7A: seed-0 vs. seeds 1-9 spatial-map r for M0-M4 |
| `Figure6B_C_support_source.csv` | Main Fig. 7B-C: per-scan spatial support; its `component_fro` column is the rank-1 contribution in Supp. Fig. S8A |
| `Figure6D_S08_LOCO_map_r_source.csv` | Main Fig. 7D and Supp. Fig. S8B: leave-one-condition-out spatial-map r |
| `M0-M4_raw_ID_crosswalk.csv` | Motif crosswalk: M0-M4 to raw K = 46 seed-0 archetype indices, earlier labels, dominant Yeo-7 networks |
| `M0-M4_cross_seed_stability_summary.csv` | Motif stability: Supp. Fig. S4 summary (map, temporal-profile and Yeo-7 r across seeds; recovered counts), text results on recovery |
| `K46_M0-M4_supplementary_table.csv` | Motif stability: Supp. Figs. S4 and S6, recovery and selected-and-recovered counts per motif and condition at K = 46 |
| `M0-M4_selection_by_K.csv` | Motif stability: M0-M4 recovery and selection at K = 45-55 by condition and m (context for Supp. Figs. S5-S6; no single panel plots it directly) |
| `K45-55_selection_stability.csv` | Supp. Fig. S5a-b: recovery-aware top-5/top-10 selection stability, K = 45-55, with the full set of diagnostics |
| `S07_atlas_mapping_audit_summary.csv` | Supp. Fig. S7: atlas mapping audit (node counts, distances, agreement rates in the caption) |
| `S09_full_model_accuracy_source.csv` | Supp. Fig. S9: full-model held-out accuracy vs. K; denominators for the retained-fraction measures |
| `S10_edge_trim_source.csv` | Supp. Fig. S10A-B: edge-trim accuracy and fraction preserved; text result on run boundaries |
| `S11_four_condition_msaa.csv` | Supp. Fig. S11: four-condition sensitivity, spatial and temporal MS-AA |
| `S11_four_condition_pca.csv` | Supp. Fig. S11: four-condition sensitivity, PCA |
| `S11_four_condition_raw_htfa.csv` | Supp. Fig. S11: four-condition sensitivity, raw HTFA reference |
| `S11_new_vs_legacy_comparison.csv` | Supp. Fig. S11: check of the current-pipeline rerun against the earlier version (caption/Supplementary Results) |
| `S12_spatial_topm_landscape.csv` | Supp. Fig. S12: full spatial K x m held-out landscape (accuracy and retained fraction, seed mean and SD) |
| `S12_landscape_run_validation.txt` | Supp. Fig. S12: validation record of the high-m landscape run |
| `S13_simulation_sensitivity_source.csv` | Supp. Fig. S13: plotted simulation-sensitivity values |
| `S13_simulation_sensitivity_runs.csv` | Supp. Fig. S13: every fit (680), including the combined-mixing check in the text |

Not covered here: main Figs. 1-2 are schematics, main Fig. 4 (simulation benchmark) has no table in
this set, and Supp. Figs. S2C and S10C-D are drawn from results that have no separate table.
