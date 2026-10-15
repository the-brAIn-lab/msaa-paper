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
