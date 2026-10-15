# Analysis notebooks

The numbered folders reflect scientific dependency rather than the chronology
of exploratory notebook creation.

- `01_fit/`: empirical spatial/temporal MS-AA model fitting.
- `02_decode/`: reconstruction and decoding analyses.
- `03_core_results/`: matched-ratio analyses, held-out sparse cross-fitting,
  PCA comparison, information concentration, and clustering-related results.
- `04_interpretation/`: representative K=46 archetype interpretation/support.
- `05_robustness/`: edge trimming and network-localization/null analyses.
- `06_four_condition/`: paragraph-scrambled four-condition sensitivity branch.
- `07_simulations/`: simulation benchmarks and alternative PCA/ICA controls.

The reusable `scripts/run_pipeline.py` implementation covers fitting, full/single-archetype decoding, and the canonical participant-held-out top-m cross-fit used for quick verification. These notebooks preserve the original paper-specific
analysis logic for the later stages.

Historical same-data top-m analyses are not canonical and are retained only under `legacy/original_notebooks/`.
