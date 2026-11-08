"""Bootstrap 95% CI for the K at which held-out top-m accuracy peaks (every K), per condition x top-m.

Resampling: outer splits (20, with replacement; the two directions of a split stay together) and, when the file has
several fit seeds, the fit seeds too (with replacement). Each bootstrap curve = mean over the resampled seeds and
splits, smoothed with a centred rolling mean over W consecutive K; its argmax is one draw of the peak K.
CI = 2.5-97.5 percentiles of those draws. Point estimate = argmax of the smoothed full-data curve.
Caveat: all splits/seeds reuse the same 36 participants, so this is split + initialisation uncertainty,
not uncertainty over new participants.
Usage: python -I peakk_ci.py HELDOUT_DIR OUT.csv [W=21] [B=2000]
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

D = Path(sys.argv[1]); OUT = Path(sys.argv[2])
W = int(sys.argv[3]) if len(sys.argv) > 3 else 21
B = int(sys.argv[4]) if len(sys.argv) > 4 else 2000
rng = np.random.default_rng(20261104)
R = pd.read_csv(D / "heldout_topm_runs_all_seeds.csv.gz",
                usecols=["top_m", "mean", "condition", "K", "outer_rep", "fit_seed", "analysis_type"])
AT = R.analysis_type.iloc[0]
R = R[R.top_m.isin([1, 2, 3, 5, 10])]
R = R.groupby(["condition", "top_m", "fit_seed", "K", "outer_rep"])["mean"].mean()  # average the 2 directions


def smooth_argmax(curves, Ks):
    # curves: (n, nK); centred rolling mean with shrinking window at the edges
    c = np.cumsum(np.pad(curves, ((0, 0), (1, 0))), axis=1)
    h = W // 2; idx = np.arange(curves.shape[1])
    lo = np.clip(idx - h, 0, None); hi = np.clip(idx + h + 1, None, curves.shape[1])
    sm = (c[:, hi] - c[:, lo]) / (hi - lo)
    return Ks[sm.argmax(axis=1)]


rows = []
for (cond, m), s in R.groupby(level=["condition", "top_m"]):
    A = s.droplevel(["condition", "top_m"]).unstack("K")          # index (seed, rep), columns K
    A = A.loc[:, A.columns > m]
    Ks = A.columns.to_numpy()
    seeds = A.index.get_level_values(0).unique().to_numpy(); reps = A.index.get_level_values(1).unique().to_numpy()
    X = np.stack([A.loc[sd].reindex(reps).to_numpy() for sd in seeds])   # (nseed, nrep, nK)
    point = smooth_argmax(np.nanmean(X, axis=(0, 1))[None], Ks)[0]
    draws = np.empty(B, dtype=int)
    for b in range(B):
        si = rng.integers(0, len(seeds), len(seeds)) if len(seeds) > 1 else np.zeros(1, int)
        ri = rng.integers(0, len(reps), len(reps))
        draws[b] = smooth_argmax(np.nanmean(X[si][:, ri], axis=(0, 1))[None], Ks)[0]
    lo, hi = np.percentile(draws, [2.5, 97.5])
    vals, cnt = np.unique(draws, return_counts=True)
    rows.append(dict(analysis_type=AT, condition=cond, top_m=m, n_seeds=len(seeds), W=W, B=B, peak_K=int(point),
                     ci_lo=int(np.floor(lo)), ci_hi=int(np.ceil(hi)), draw_median=int(np.median(draws)),
                     share_in_ci=float(((draws >= lo) & (draws <= hi)).mean()),
                     modal_K=int(vals[cnt.argmax()]), modal_share=float(cnt.max() / B)))
    print(rows[-1], flush=True)
pd.DataFrame(rows).to_csv(OUT, index=False)
print("wrote", OUT)
