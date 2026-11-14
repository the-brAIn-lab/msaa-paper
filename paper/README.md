# Paper figures

`main_figures/` and `supplement/` hold the notebooks and scripts that draw Figures 2, 3, 5, 6 and
S07, S08, S10, S11, S13 from saved result tables (the rest are drawn by
`analysis/08_multiseed_stability/`). None of them fit or decode.

`scripts/build_figures.py` runs all of them and assembles the figure bundle; `docs/REPRODUCE.md`
maps each figure to its script and inputs.

- `FIGURE_CAPTIONS.tex`: LaTeX captions for every figure (`\caption{...}` + `\label{...}`), with the uncertainty definitions and conventions.
- `FIGURE_BUNDLE_README.txt`: the README that goes into the bundle.
