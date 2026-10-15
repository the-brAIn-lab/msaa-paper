# Archetypes manuscript — supplementary figure notebooks

These notebooks are figure-only layers. They are intended to sit on top of the
cleaned analysis notebooks and exported CSV/checkpoint outputs.

## Included notebooks

1. `supp_01_PCA_validation_condensed.ipynb`
   - PCA projected vs reconstructed decoding, all conditions overlaid.

2. `supp_02_joint_intact_word_convergence.ipynb`
   - joint intact/word held-out convergence moved out of the main Figure 4.

3. `supp_03_full_across_within_ratio_landscapes.ipynb`
   - full MS-AA decoding across fit scope and orientation.

4. `supp_04_same_data_vs_heldout_topm.ipynb`
   - same-data sparse ranking/evaluation vs participant-held-out cross-fitting.

5. `supp_05_edge_trim_diagnostics.ipynb`
   - detailed K=46 boundary-trim diagnostics.

6. `supp_06_paragraph_four_condition_validation.ipynb`
   - paragraph-scrambled four-condition sensitivity, absolute accuracy and
     fractional retention.

7. `supp_07_network_localization_nulls.ipynb`
   - strict random-archetype null and matched-ratio localization/network
     comparisons.

8. `supp_08_clustering_geometry.ipynb`
   - cached clustering score, purity/balance, and decoding-vs-clustering
     summary.

9. `supp_09_simulation_overlap_dynamics.ipynb`
   - finalized overlap / dynamical-mixing simulation sensitivity.

10. `supp_10_simulation_robustness.ipynb`
    - finalized noise / heterogeneity / fitted-K simulation robustness.

11. `supp_11_simulation_temporal.ipynb`
    - finalized temporal-recovery calibration and metastable paired difference.

## Design principle

The notebooks do not refit models or silently change the analyses. They read
the analysis outputs and create publication-ready figures. This keeps the
reproducibility layer separate from figure styling.
