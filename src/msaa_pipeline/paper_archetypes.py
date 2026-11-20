"""The five K=50 spatial archetypes interpreted in the paper (M0-M4), in either numbering.

M0-M4 are the five archetypes of the representative K=50 spatial across-condition fit that
the paper interprets. They were chosen by participant-held-out top-5 selection frequency
averaged over the 10 fit seeds (each aligned to seed 0 by Hungarian matching), in the best
condition, after excluding archetypes not reproducibly recovered across seeds (canonical 31,
recovered in 5/9 seeds); see analysis/08_multiseed_stability/DEFINITIONS.md Sec. 4. The ids
below are seed-0 ids. Ranking (10-seed mean frequency, best condition): 28 .995 word, 16 .885
intact, 39 .698 word, 49 .622 intact, 37 .508 intact; next 36 .460, 31 .393 (excluded).
Archetype numbers depend on how the fit was made:

  canonical  fits with frozen FurthestSum picks and no post-fit sort: scripts/run_parallel.py,
             and scripts/run_pipeline.py --profile full (make full-fit). This is the numbering
             used in the paper text and the multi-seed analysis.
  legacy     the original analysis outputs (reference_outputs/legacy_full/)

A fit without recorded FurthestSum picks is treated as legacy. There are no legacy-numbered
K=50 paper ids, so a legacy K=50 fit is refused.

  label  canonical K=50
  M0     28
  M1     39
  M2     16
  M3     37
  M4     49

Usage in a notebook (run from the repo root):

    import sys; sys.path.insert(0, "src")
    from msaa_pipeline.paper_archetypes import paper_archetypes
    SELECTED_ARCHETYPES, LABELS = paper_archetypes()          # reads the K=50 fit in msaa_flexible_outputs_npz/
    LABEL = dict(zip(SELECTED_ARCHETYPES, LABELS))
"""
from __future__ import annotations

from pathlib import Path

import numpy as np

PAPER_K = 50
LABELS = ["M0", "M1", "M2", "M3", "M4"]
IDS = {"canonical": [28, 39, 16, 37, 49]}
DEFAULT_FIT = Path("msaa_flexible_outputs_npz") / f"spatialAA_across_acrossCond_K{PAPER_K}.npz"


def numbering(fit_file: str | Path = DEFAULT_FIT) -> str:
    """'canonical' if the fit records its frozen FurthestSum picks (init_index), else 'legacy'."""
    p = Path(fit_file)
    if not p.exists():
        raise FileNotFoundError(f"{p} not found: run the K={PAPER_K} spatial across fit first "
                                "(scripts/run_parallel.py fit --jobs spatial_across --ks 50)")
    with np.load(p, allow_pickle=True) as f:
        return "canonical" if "init_index" in f.files else "legacy"


def paper_archetypes(fit_file: str | Path = DEFAULT_FIT) -> tuple[list[int], list[str]]:
    """(archetype ids in the fit's own numbering, labels M0-M4)."""
    n = numbering(fit_file)
    if n not in IDS:
        raise ValueError(f"{fit_file}: no {n}-numbered M0-M4 ids at K={PAPER_K}; use a canonical (frozen-pick) fit")
    return list(IDS[n]), list(LABELS)


def label_map(fit_file: str | Path = DEFAULT_FIT) -> dict[int, str]:
    ids, labels = paper_archetypes(fit_file)
    return dict(zip(ids, labels))
