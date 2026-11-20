# Paper figures

`main_figures/` and `supplement/` hold the notebooks and scripts that draw Figures 2, 3, 5, 6, 8 and
S07, S08, S10, S11, S13 from saved result tables (the rest are drawn by
`analysis/08_multiseed_stability/`). None of them fit or decode. Figure 8 (panels A-F) is drawn by
`main_figures/figure08_support_robustness.py`, which adds the retention (E) and edge-trim (F)
panels to the Figure 6 panels A-D.

These numbers are the figure-bundle numbering. In the manuscript: Fig 2 = main Fig 3,
Fig 3 = main Fig 4, Fig 4A/B = main Fig 5A/B, Fig 4C = main Fig 6A, Fig 4D = main Fig 6B,
Fig 5 = main Fig 7, Fig 8 = main Fig 8 (Fig 6 is its panels A-D). Main Figs. 1-2 are schematics.
Supplementary numbers are the same in both.

`scripts/build_figures.py` runs all of them and assembles the figure bundle; `docs/REPRODUCE.md`
maps each figure to its script and inputs.

- `final_figures/`: every figure exactly as published, under manuscript numbers (EPS + PDF).
- `source_tables/`: the data behind each figure (see its README).
- `FIGURE_BUNDLE_README.txt`: the README that goes into the bundle.
