#!/usr/bin/env python3
"""K<K>_top<m>.csv from selection_by_seed.csv: held-out top-m selection frequency of every seed-0 archetype
in each condition and fit seed, plus the 10-seed mean and the number of seeds >= 0.5 (the table hk_merge.py
writes for the representative K; the Figure 5 notebook and the M0-M4 selection rule read it at the representative K).

Usage: python hk_topm_table.py K [--top-m 5] [--sel selection_by_seed.csv] [--out DIR]
Defaults (paths.py): heldout-everyk-spatial-all-results/selection_by_seed.csv, written next to it.
"""
import argparse
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))   # analysis/08_multiseed_stability, for paths.py
import paths  # noqa: E402

ap = argparse.ArgumentParser(); ap.add_argument("K", type=int); ap.add_argument("--top-m", type=int, default=5)
ap.add_argument("--sel", default=None); ap.add_argument("--out", default=None); a = ap.parse_args()
sel_f = Path(a.sel or paths.heldout_everyk("spatial") / "selection_by_seed.csv")
out = Path(a.out or sel_f.parent); out.mkdir(parents=True, exist_ok=True)
sel = pd.read_csv(sel_f)
k = sel[(sel.K == a.K) & (sel.top_m == a.top_m)]
assert k.seed.nunique() == paths.N_SEEDS, sorted(k.seed.unique())
t = k.pivot_table(index=["archetype0", "condition"], columns="seed", values="selection_frequency", fill_value=0)
t["mean"] = t.mean(axis=1); t["seeds_ge_0.5"] = (t.drop(columns="mean") >= 0.5).sum(axis=1)
t = t.sort_values("mean", ascending=False); t.to_csv(out / f"K{a.K}_top{a.top_m}.csv")
print(t.head(12).round(3).to_string())
