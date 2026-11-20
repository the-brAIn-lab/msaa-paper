paper figures (figure-bundle numbering; manuscript numbers in brackets)

each figure is a PNG and a vector PDF, drawn at 5" width
the source data behind each figure are in tables/
all held-out numbers use the 1e-8 tie rule.
10-seed results average the two directions within each outer split, then the 20 splits within a seed, then the 10 fit seeds.

Main figures (main-png/, main-pdf/)

    Figure2     [main Fig 3]   PCA vs spatial and temporal MS-AA, with raw HTFA and chance
    Figure3     [main Fig 4]   simulation: motif recovery vs span preservation
    Figure4A    [main Fig 5A]  spatial held-out top-m accuracy
    Figure4B    [main Fig 5B]  temporal held-out top-m accuracy
    Figure4C    [main Fig 6A]  spatial minus temporal at matched component ratios (intact)
    Figure4D    [main Fig 6B]  temporal retained accuracy at matched fractions of K
    Figure5     [main Fig 7]   selected spatial motifs M0-M4 at K = 50
    Figure6     [main Fig 8 A-D] support and robustness of M0-M4
    Figure8     [main Fig 8]   support and robustness of M0-M4, panels A-F
                               (adds E, K = 50 retention, and F, edge trim)

Main Figs 1-2 are schematics and are not in the bundle.

Supplement (supplement-png/, supplement-pdf/)

    S01     spatial minus temporal at matched ratios, all conditions
    S02     temporal K = 14: accumulation, generator-anchor masking, higher-K splitting
    S03     temporal K = 12/14/16 per-seed contrasts
    S04     M0-M4 across further fit seeds
    S05a    K = 45-55 selection stability, top-5
    S05b    K = 45-55 selection stability, top-10
    S06     M0-M4 selected and recovered in each seed's top-5 / top-10
    S07     atlas mapping audit
    S08     condition contributions
    S09     full-model held-out accuracy vs K
    S10     edge-trim diagnostics
    S11     four-condition sensitivity (adds paragraph-scrambled)
    S12     full participant-held-out spatial K x m landscape
    S13     simulation sensitivity of motif recovery

tables/

    M0-M4_raw_ID_crosswalk.csv is the only place raw archetype IDs appear
    the other files are the plotted data, named by figure.
