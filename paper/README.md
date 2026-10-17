# Paper figure notebooks

`main_figures/` contains the current Figure 2–6 notebooks. Figure 1 is conceptual
and has no analysis dependency in the supplied code snapshot.

`supplement/` contains the cleaned/curated supplementary plotting and analysis
notebooks from the code snapshot.

Most paper notebooks intentionally consume saved analysis tables rather than
refitting models. This keeps plotting changes separate from scientific fitting
and statistical analysis.

Run a nested notebook from the repository root with:

```bash
python scripts/execute_notebook.py paper/main_figures/<notebook>.ipynb
```

The required upstream tables are documented in `docs/PIPELINE.md` and the
initial dependency audit.

The manuscript sparse-information figures must use `msaa_heldout_topm_outputs_spatial_across/`; same-data top-m outputs are not part of the cleaned paper pipeline.

## Generated figure outputs

Fresh main-text figures are written under `outputs/paper_figures/`. Fresh
supplementary figures and source tables are written under
`outputs/paper_supplement/`. These generated outputs are intentionally ignored by
Git and are kept separate from the read-only validated legacy store.
