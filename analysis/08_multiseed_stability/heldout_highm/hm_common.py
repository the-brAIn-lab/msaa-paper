"""Shared settings for the high-m held-out run (Fig 4C / S01 and the S12 K x m landscape, 10 seeds).

Scores extra top-m values at the 15 exact spatial/temporal matched-ratio K pairs of
the matched-ratio heatmaps, reusing the per-split archetype rankings the every-K held-out run
already saved (no refitting, no re-ranking). Splits, evaluation participants and decoder
seeds are the every-K run's, so m=10 and the full model must reproduce it exactly."""
import os
from pathlib import Path

CONDS = ["intact", "word", "rest"]
NREP = 20
NFOLDS, NREPS, BASE_SEED = 2, 20, 42
OLD_SUB = "msaa_heldout_topm_outputs_spatial_across"     # 30-K held-out results (spatial), reused by every-K

# (spatial K, temporal K) with exactly equal K/700 and K/300 (fig4c_matched_ratio_heatmaps.py EXACT_PAIRS)
EXACT_PAIRS = [(7, 3), (14, 6), (21, 9), (35, 15), (42, 18), (63, 27), (70, 30), (105, 45), (126, 54),
               (140, 60), (175, 75), (189, 81), (210, 90), (252, 108), (700, 300)]
NEW_M = [7, 15, 20, 25, 30, 40, 50, 75, 100]    # heatmap m values not in the every-K run (1,2,3,5,10)
CHECK_M = 10                                     # rescored too, as the exact-match check


# HM_GRID selects what to score:
#   pairs      (default) the 15 matched spatial/temporal K pairs, m up to 100 (Fig 4C / S01)
#   landscape  spatial only, K = 10, 20, ..., 700, m = 7, 15, 20, 25, 30, 40, 50 (S12 K x m landscape)
GRID = os.environ.get("HM_GRID", "pairs")
LANDSCAPE_KS = list(range(10, 701, 10))
LANDSCAPE_M = [7, 15, 20, 25, 30, 40, 50]


def m_values(at, K):
    """m to score at this K (plus m=10 as the exact-match check)."""
    if GRID == "landscape":
        return [m for m in sorted(LANDSCAPE_M + [CHECK_M]) if m <= K]
    lim = {p[0] if at == "spatial" else p[1]: min(p) for p in EXACT_PAIRS}[K]   # m must fit both K of the pair
    return [m for m in sorted(NEW_M + [CHECK_M]) if m <= lim]


def pair_ks(at):
    if GRID == "landscape":
        return LANDSCAPE_KS if at == "spatial" else []
    return [p[0] if at == "spatial" else p[1] for p in EXACT_PAIRS]


def seed_root(run, s):
    return Path(run) if s == 0 else Path(run) / "seeds" / "full" / f"seed{s:02d}"


def everyk_split_dir(run, at, s, K, c, r):
    return Path(run) / "heldout-everyk" / "parts" / at / f"seed{s:02d}" / f"K{K}_{c}_r{r}"


def base_dir(run, mode):
    name = "heldout-highm" + ("-landscape" if GRID == "landscape" else "")
    return Path(run) / (name if mode == "full" else name + "-smoke")


def out_dir(base, at, s, K, c, r):
    return Path(base) / "parts" / at / f"seed{s:02d}" / f"K{K}_{c}_r{r}"


def tie_tol_ok():
    return os.environ.get("MSAA_TIE_TOL", "") not in ("", "0", "0.0")
