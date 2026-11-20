# Source tables

These are the exact tables behind the published figures, copied unchanged from the figure bundle of
17 Nov 2026 (`tables/`), with every motif-dependent table (main Figs. 7-8, Supp. Figs. S4-S8 and
S10) taken from the K = 50 figure bundle of 19 Nov 2026. Figure and supplementary-figure numbers below follow the
manuscript. File names keep the older bundle numbering (for example `Figure2_*` is main Figure 3,
`Figure4*` is main Figures 5-6, `Figure5_*` is main Figure 7 and `Figure6*` is main Figure 8).

Archetype indices (the `archetype`, `seed0_id` and similar columns) refer to the spatial MS-AA
fit at K = 50, fit seed 0, and are 0-based. They are arbitrary up to permutation: a different seed
or software version can number the same archetypes differently. Motif labels M0-M4 are the stable
names; `M0-M4_raw_ID_crosswalk.csv` maps them to the raw indices.

| File | Underlies |
| --- | --- |
| `Figure2_source.csv` | Main Fig. 3: empirical benchmark, PCA reconstruction vs. spatial and temporal MS-AA, with raw HTFA reference |
| `Figure4A_source.csv` | Main Fig. 5A: spatial MS-AA held-out top-m accuracy vs. K (mean, SD over 10 fit seeds) |
| `Figure4B_source.csv` | Main Fig. 5B: temporal MS-AA held-out top-m accuracy vs. K |
| `Figure4D_source.csv` | Main Fig. 6B: temporal MS-AA retained fraction at 25/50/70% of K |
| `temporal_matched_fraction_seed_level.csv` | Main Fig. 6B and Supp. Fig. S3: per-seed values behind the retained-fraction curves |
| `temporal_K14_local_contrast_delta14.csv` | Supp. Fig. S3: K = 12/14/16 local contrast (Delta14), per seed with bootstrap CIs; text result on K = 14 |
| `temporal_K14_anchor_masked_vs_unmasked.csv` | Supp. Fig. S2A-B: temporal K = 14 accumulation over m and generator-anchor masking |
| `spatial_vs_temporal_matched_ratio_summary.csv` | Main Fig. 6A and Supp. Fig. S1: spatial minus temporal accuracy at matched ratios, seed means and seed counts |
| `spatial_vs_temporal_matched_ratio_per_seed.csv` | Main Fig. 6A and Supp. Fig. S1: per-seed comparisons (4,410 rows) behind the summary |
| `Figure5_temporal_profiles.csv` | Main Fig. 7: temporal profiles of M0-M4 (mean and SEM over scans) |
| `Figure5_Yeo7_allnode_mass.csv` | Main Fig. 7: Yeo-7 composition of M0-M4, all 700 nodes; text results on network shares |
| `Figure5_crossfit_selection_frequency.csv` | Main Fig. 7: cross-fitted top-5 selection frequency by condition (10-seed mean) |
| `K50_all_archetype_selection_ranking.csv` | Methods (selection of M0-M4): held-out top-5 selection frequency of every K = 50 archetype by condition, per fit seed and 10-seed mean (seed-0 ids in `archetype0`); the five motifs, the next-ranked archetype (0.46) and the excluded, non-recovered archetype (raw id 31, seventh) are read from it |
| `Figure5_spatial_maps_allnodes.csv` | Main Fig. 7: participant-averaged spatial maps of M0-M4 on all 700 nodes, with node coordinates and Yeo-7 labels; source of the mean pairwise cosine (0.17, maximum 0.44) in the Supplementary Results |
| `Figure6A_cross_seed_map_r_source.csv` | Main Fig. 8A: seed-0 vs. seeds 1-9 spatial-map r for M0-M4 |
| `Figure6B_C_support_source.csv` | Main Fig. 8B-C: per-scan spatial support; its `component_fro` column is the rank-1 contribution in Supp. Fig. S8A |
| `Figure6D_S08_LOCO_map_r_source.csv` | Main Fig. 8D and Supp. Fig. S8B: leave-one-condition-out spatial-map r |
| `M0-M4_raw_ID_crosswalk.csv` | Motif crosswalk: M0-M4 to raw K = 50 seed-0 archetype indices, dominant Yeo-7 networks |
| `M0-M4_cross_seed_stability_summary.csv` | Motif stability: Supp. Fig. S4 summary (map, temporal-profile and Yeo-7 r across seeds; recovered counts), text results on recovery |
| `K50_M0-M4_supplementary_table.csv` | Motif stability: Supp. Figs. S4 and S6, recovery and selected-and-recovered counts per motif and condition at K = 50 |
| `Figure8E_K50_retention_source.csv` | Main Fig. 8E: held-out top-3/5/10 retention at K = 50 (10-seed mean and SD) |
| `Figure8F_K50_edge_trim_source.csv` | Main Fig. 8F: edge-trim fraction of above-chance accuracy preserved at K = 50 |
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
| `S13_simulation_sensitivity_runs.csv` | Supp. Fig. S13: every fit (680), including the combined-mixing check in the text. The `diag_*` columns are run diagnostics that the paper does not use; `diag_mean_dwell` and `diag_median_dwell` in particular do not track the dwell setting and should not be interpreted |

Not covered here: main Figs. 1-2 are schematics, main Fig. 4 (simulation benchmark) has no table in
this set, and Supp. Figs. S2C and S10C-D are drawn from results that have no separate table.
