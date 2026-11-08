#!/usr/bin/env python3
"""Combine each seed directory's picks parts into its configs/furthestsum_init.csv.gz (the file
the seqnum code's load_frozen_init reads), and check every stage fit has K picks.
  python picks_merge.py STAGE
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import pandas as pd  # noqa: E402
from seeds_common import fits_of, seed_dir, stage_seeds  # noqa: E402

KEY = ["analysis_type", "fit_scope", "condition", "K"]


def main():
    stage = sys.argv[1]
    want = fits_of(stage)
    bad = 0
    for s in stage_seeds(stage):
        d = seed_dir(stage, s)
        parts = sorted((d / "picks_parts").glob("*.csv"))
        if not parts:
            print(f"seed {s}: no picks parts"); bad += 1; continue
        df = pd.concat([pd.read_csv(p) for p in parts], ignore_index=True)
        df = df.sort_values(KEY + ["archetype"]).reset_index(drop=True)
        n = df.groupby(KEY).size()
        missing = [f for f in want if f[0] == s and
                   n.get((f[1], f[2], f[4] or "acrossCond", f[3]), 0) != f[3]]
        df.to_csv(d / "configs" / "furthestsum_init.csv.gz", index=False)
        print(f"seed {s}: {n.size} fits with picks, missing {len(missing)}")
        bad += bool(missing)
    print("picks OK" if not bad else f"picks INCOMPLETE in {bad} seed(s)")
    sys.exit(1 if bad else 0)


if __name__ == "__main__":
    main()
