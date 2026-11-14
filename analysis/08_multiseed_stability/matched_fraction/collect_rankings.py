#!/usr/bin/env python3
"""Collect the temporal held-out per-archetype rankings temporal_matched_fraction.py needs into ONE csv.

    python collect_rankings.py RUN_DIR OUT_CSV KS CONDS SEEDS [NREP]

RUN_DIR is the run root (the repo root after run_local.py). KS, CONDS, SEEDS are ':'-separated,
e.g. 5:8:10  intact:word  0:1:2.
Source: the every-K held-out kit's per-split files
    RUN_DIR/heldout-everyk/parts/temporal/seedNN/K<K>_<cond>_r<rep>/ranks.csv
(one row per archetype x direction; these are exactly the rankings the held-out top-m used,
including the outer_seed of each participant split). Adds a fit_seed column so each row is
mapped to the right seed's fit. Exits 2 if anything is missing.
"""
import sys
from pathlib import Path

import pandas as pd

def main():
    run, out = Path(sys.argv[1]), Path(sys.argv[2])
    ks = [int(x) for x in sys.argv[3].split(":")]
    conds = sys.argv[4].split(":")
    seeds = [int(x) for x in sys.argv[5].split(":")]
    NREP = int(sys.argv[6]) if len(sys.argv) > 6 else 20   # outer splits 0..NREP-1 (smoke uses 2)
    parts = run / "heldout-everyk" / "parts" / "temporal"
    frames, missing, bad = [], [], []
    for s in seeds:
        for K in ks:
            for c in conds:
                for r in range(NREP):
                    f = parts / f"seed{s:02d}" / f"K{K}_{c}_r{r}" / "ranks.csv"
                    if not f.exists():
                        missing.append(str(f))
                        continue
                    x = pd.read_csv(f)
                    n = x.groupby("direction").archetype.nunique()
                    if len(n) != 2 or (n != K).any():
                        bad.append(f"{f}: {n.to_dict()}")
                    frames.append(x.assign(fit_seed=s))
    print(f"read {len(frames)} ranking files ({len(seeds)} seeds x {len(ks)} K x {len(conds)} conditions x {NREP} splits"
          f" = {len(seeds) * len(ks) * len(conds) * NREP} expected)")
    if missing or bad:
        print(f"MISSING {len(missing)}, INCOMPLETE {len(bad)}; first few:")
        for line in (missing + bad)[:10]:
            print("  ", line)
        sys.exit(2)
    df = pd.concat(frames, ignore_index=True)
    out.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out, index=False)
    print(f"wrote {out} ({len(df)} rows)")


if __name__ == "__main__":
    main()
