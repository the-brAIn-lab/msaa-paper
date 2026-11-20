# Analyses

Each folder holds what one part of the paper needs; `docs/REPRODUCE.md` says which output each
makes and which figure uses it.

- `03_core_results/`: PCA vs MS-AA benchmark tables (Fig 2), temporal K = 14 accumulation and
  generator-anchor masking (S02).
- `04_interpretation/`: K = 50 motif support (Fig 6B/C), leave-one-condition-out refits and
  condition deployment (Fig 6D, S08).
- `05_robustness/`: edge-trim sensitivity of the K = 50 held-out top-5 (S10).
- `06_four_condition/`: four-condition (paragraph-scrambled) sensitivity (S11).
- `07_simulations/`: simulation benchmark (Fig 3) and its sensitivity sweep (S13).
- `08_multiseed_stability/`: the seed-0 and 10-seed runs (fits, decoding, held-out top-m at every
  K and at high m), the matched-fraction analysis and the cross-seed stability tables, plus the
  scripts that draw Fig 4, S01-S06, S09 and S12 (see its README).

Notebooks run from the repo root: `uv run python scripts/execute_notebook.py <notebook>`.
