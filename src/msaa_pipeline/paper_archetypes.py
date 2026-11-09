"""The five K=46 spatial archetypes interpreted in the paper (M0-M4), in either numbering.

M0-M4 are the five archetypes of the representative K=46 spatial across-condition fit that
the paper interprets, chosen with participant-held-out top-5 selection. Archetype numbers
depend on how the fit was made:

  canonical  fits with frozen FurthestSum picks and no post-fit sort: scripts/run_parallel.py,
             and scripts/run_pipeline.py --profile full (make full-fit). This is the numbering
             used in the paper text and the multi-seed analysis.
  legacy     the original analysis outputs (reference_outputs/legacy_full/)

A fit without recorded FurthestSum picks is treated as legacy. Legacy-style refits
(run_pipeline.py --numbering legacy, or the mini-paper profile) sort archetypes after fitting
but are not guaranteed to reproduce the legacy order; check them with
scripts/compare_k46_archetypes.py.

The same five archetypes in both numberings (matched maps r >= 0.9999):

  label  canonical  legacy
  M0     28         4
  M1     35         5
  M2     16         7
  M3     39         8
  M4     45         9

Usage in a notebook (run from the repo root):

    import sys; sys.path.insert(0, "src")
    from msaa_pipeline.paper_archetypes import paper_archetypes
    SELECTED_ARCHETYPES, LABELS = paper_archetypes()          # reads the K=46 fit in msaa_flexible_outputs_npz/
    LABEL = dict(zip(SELECTED_ARCHETYPES, LABELS))
"""
from __future__ import annotations

from pathlib import Path

import numpy as np

PAPER_K = 46
LABELS = ["M0", "M1", "M2", "M3", "M4"]
IDS = {"canonical": [28, 35, 16, 39, 45], "legacy": [4, 5, 7, 8, 9]}
DEFAULT_FIT = Path("msaa_flexible_outputs_npz") / f"spatialAA_across_acrossCond_K{PAPER_K}.npz"


def numbering(fit_file: str | Path = DEFAULT_FIT) -> str:
    """'canonical' if the fit records its frozen FurthestSum picks (init_index), else 'legacy'."""
    p = Path(fit_file)
    if not p.exists():
        raise FileNotFoundError(f"{p} not found: run the K={PAPER_K} spatial across fit first "
                                "(scripts/run_parallel.py fit --jobs spatial_across --ks 46)")
    with np.load(p, allow_pickle=True) as f:
        return "canonical" if "init_index" in f.files else "legacy"


def paper_archetypes(fit_file: str | Path = DEFAULT_FIT) -> tuple[list[int], list[str]]:
    """(archetype ids in the fit's own numbering, labels M0-M4)."""
    return list(IDS[numbering(fit_file)]), list(LABELS)


def label_map(fit_file: str | Path = DEFAULT_FIT) -> dict[int, str]:
    ids, labels = paper_archetypes(fit_file)
    return dict(zip(ids, labels))
