#!/usr/bin/env python3
"""Check the matched-fraction numbers against the existing held-out top-m where they overlap.

    python validate.py RUN_DIR RESULTS_DIR

Whenever round(fraction*K) is one of the held-out top-m values (1,2,3,5,10), e.g. K=14 at 70% -> m=10,
K=20 at 50% -> m=10, K=10 at 50% -> m=5, the same split, ranking, evaluation participants and decoder
seed were already scored by the every-K held-out run:
    RUN_DIR/heldout-everyk/merged/temporal/seedNN/heldout_topm_crossfit_runs.csv.gz
Those cells should agree up to rounding (temporal_matched_fraction.py decodes dense, the pipeline
low-rank for m<=10).
Writes RESULTS_DIR/validation.csv + validation.txt and prints the summary. Exit 1 if a cell differs
by more than 1e-3 or a selected set differs.
"""
import ast
import glob
import sys
from pathlib import Path

import pandas as pd

KEY = ["condition", "K", "outer_rep", "direction", "top_m"]


def main():
    run, res = Path(sys.argv[1]), Path(sys.argv[2])
    files = sorted(glob.glob(str(res / "tasks" / "*.csv")))
    if not files:
        print("no task outputs yet"); sys.exit(1)
    new = pd.concat([pd.read_csv(f) for f in files], ignore_index=True)
    rows = []
    for s, sub in new.groupby("fit_seed"):
        p = run / "heldout-everyk" / "merged" / "temporal" / f"seed{int(s):02d}" / "heldout_topm_crossfit_runs.csv.gz"
        if not p.exists():
            print(f"seed {s}: no existing held-out runs at {p}"); continue
        old = pd.read_csv(p)
        old = old[old.analysis_type == "temporal"] if "analysis_type" in old else old
        m = sub[sub.top_m < sub.K].merge(old[KEY + ["mean", "full_mean_same_eval", "selected_archetypes"]],
                                         on=KEY, suffixes=("", "_old"))
        m["fit_seed"] = s
        rows.append(m)
    if not rows:
        print("nothing to compare"); sys.exit(1)
    m = pd.concat(rows, ignore_index=True)
    m["d_mean"] = (m["mean"] - m["mean_old"]).abs()
    m["d_full"] = (m["full_mean_same_eval"] - m["full_mean_same_eval_old"]).abs()
    m["same_set"] = [sorted(ast.literal_eval(str(a))) == sorted(ast.literal_eval(str(b)))
                     for a, b in zip(m.selected_archetypes, m.selected_archetypes_old)]
    m.to_csv(res / "validation.csv", index=False)
    ok = m.d_mean.max() <= 1e-3 and m.d_full.max() <= 1e-3 and m.same_set.all()
    by = m.groupby(["condition", "K", "top_m"]).agg(cells=("d_mean", "size"), max_d_mean=("d_mean", "max"),
                                                    max_d_full=("d_full", "max"), same_set=("same_set", "mean"))
    txt = (f"VALIDATION {'PASS' if ok else 'FAIL'}: {len(m)} overlapping cells, max |diff| top-m {m.d_mean.max():.2e}, "
           f"full model {m.d_full.max():.2e}, selected sets identical {m.same_set.mean():.0%}\n" + by.to_string() + "\n")
    (res / "validation.txt").write_text(txt)
    print(txt)
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
