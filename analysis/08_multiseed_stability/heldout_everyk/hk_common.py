"""Shared settings for the every-K held-out top-m kit (see README.md)."""
import os
from pathlib import Path

import numpy as np
import pandas as pd

CONDS = ["intact", "word", "rest"]
NREP = 20                       # outer participant splits per (condition, K), as in the paper
TOP_M = [1, 2, 3, 5, 10]        # profiles.json "full" top_m
NFOLDS, NREPS, BASE_SEED = 2, 20, 42
OLD_SUB = "msaa_heldout_topm_outputs_spatial_across"   # existing 30-K held-out results per seed
NEW_SEED_SALT = 20261102        # outer seeds for K not in configs/heldout_outer_seeds.csv

# phase -> (analysis type, seeds); spatial entirely first, then temporal (ben 2026-11-02)
PHASES = {"s0": ("spatial", [0]), "sA": ("spatial", list(range(1, 10))),
          "t0": ("temporal", [0]), "tA": ("temporal", list(range(1, 10)))}
PHASE_NAMES = {"s0": "SPATIAL SEED 0", "sA": "SPATIAL SEEDS 1-9", "t0": "TEMPORAL SEED 0", "tA": "TEMPORAL SEEDS 1-9"}
KS_FULL = {"spatial": list(range(2, 701)), "temporal": list(range(2, 301))}
KS_SMOKE = {"spatial": [3, 46, 47], "temporal": [3, 47]}     # 46 tests reuse of the existing 30-K results
SMOKE_SEEDS = {"s0": [0], "sA": [1], "t0": [0], "tA": [1]}


def phase(mode, ph):
    at, seeds = PHASES[ph]
    if mode == "smoke":
        seeds = SMOKE_SEEDS[ph]
    ks = (KS_SMOKE if mode == "smoke" else KS_FULL)[at]
    # HK_SEEDS / HK_KS (comma lists, set by analysis/08_multiseed_stability/run_local.py) restrict a full run
    if os.environ.get("HK_TYPES") and at not in os.environ["HK_TYPES"].split(","):
        seeds = []
    if os.environ.get("HK_SEEDS"):
        keep = {int(x) for x in os.environ["HK_SEEDS"].split(",")}
        seeds = [s for s in seeds if s in keep]
    if os.environ.get("HK_KS"):
        keep = {int(x) for x in os.environ["HK_KS"].split(",")}
        ks = [k for k in ks if k in keep]
    return at, seeds, ks


def base_dir(run, mode):
    return Path(run) / ("heldout-everyk" if mode == "full" else "heldout-everyk-smoke")


def seed_root(run, s):
    return Path(run) if s == 0 else Path(run) / "seeds" / "full" / f"seed{s:02d}"


def split_dir(base, at, s, K, c, r):
    return Path(base) / "parts" / at / f"seed{s:02d}" / f"K{K}_{c}_r{r}"


def fit_file(run, at, s, K):
    return seed_root(run, s) / "msaa_flexible_outputs_npz" / f"{at}AA_across_acrossCond_K{K}.npz"


def outer_seed_table(run):
    """Frozen paper seeds for the 30 K, plus deterministic new seeds for every other K (2-700)."""
    old = pd.read_csv(Path(run) / "configs" / "heldout_outer_seeds.csv")[["condition", "K", "outer_rep", "outer_seed"]]
    have = set(old.K.astype(int))
    rows = [dict(condition=c, K=K, outer_rep=r,
                 outer_seed=int(np.random.default_rng([NEW_SEED_SALT, ci, K, r]).integers(0, 2**31 - 1)))
            for K in range(2, 701) if K not in have for ci, c in enumerate(CONDS) for r in range(NREP)]
    t = pd.concat([old.assign(source="paper_frozen"), pd.DataFrame(rows).assign(source=f"new_{NEW_SEED_SALT}")])
    t = t.astype({"K": int, "outer_rep": int, "outer_seed": "int64"}).sort_values(["K", "condition", "outer_rep"])
    assert not t.duplicated(["condition", "K", "outer_rep"]).any()
    return t.reset_index(drop=True)


def load_outer_seeds(base):
    t = pd.read_csv(Path(base) / "configs" / "heldout_outer_seeds_everyk.csv")
    return {(str(x.condition), int(x.K), int(x.outer_rep)): int(x.outer_seed) for x in t.itertuples()}


def old_complete_ks(run, s):
    """K whose 60 outer splits (both directions) are already in the seed's existing spatial held-out results."""
    d = seed_root(run, s) / OLD_SUB
    p, q = d / "heldout_topm_crossfit_runs.csv", d / "heldout_per_archetype_ranking_runs.csv"
    if not (p.exists() and q.exists()):
        return set()
    ok = None
    for f in (p, q):
        x = pd.read_csv(f, usecols=["condition", "K", "outer_rep", "direction"])
        n = x.groupby(["K", "condition", "outer_rep"]).direction.nunique().eq(2).groupby("K").sum()
        ks = set(n[n == 3 * NREP].index.astype(int))
        ok = ks if ok is None else ok & ks
    return ok
