"""Leave-one-seed-out validation of "good-K" bands on held-out top-m accuracy (every K).

For each held-out seed s: average the other 9 seeds' held-out accuracy per K, smooth with a centred rolling
mean over W consecutive K, and define the good-K band = all K whose smoothed value is >= (1 - DELTA) x peak
(may be several segments). Then test the band on seed s (same smoothing):
  peak_in_band   : is seed s's own smoothed peak K inside the 9-seed band?
  jaccard        : overlap of seed s's own band (same rule) with the 9-seed band
  band_minus_out : seed s mean accuracy inside the band minus outside it
  window_pct     : percentile of seed s's band mean among all same-width sliding windows of K (1 = best window)
Usage: python goodk_loso.py ANALYSIS [--dir HELDOUT_ALL_DIR] [--out OUTDIR] [--W 9] [--deltas 0.05,0.10]
Default HELDOUT_ALL_DIR heldout-everyk-<ANALYSIS>-all-results/ (needs seedNN_heldout_topm_summary.csv for every seed),
OUTDIR outputs/multiseed_stability/goodk_loso/<ANALYSIS>_W<W>/.
"""
import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))   # analysis/08_multiseed_stability, for paths.py
import paths  # noqa: E402
ap = argparse.ArgumentParser(); ap.add_argument("analysis", choices=["spatial", "temporal"])
ap.add_argument("--dir"); ap.add_argument("--out"); ap.add_argument("--W", type=int, default=9)
ap.add_argument("--deltas", default="0.05,0.10")
a = ap.parse_args()
W = a.W
DELTAS = [float(x) for x in a.deltas.split(",")]
D = Path(a.dir or paths.heldout_everyk(a.analysis, "all"))
OUT = Path(a.out or paths.OUT / "goodk_loso" / f"{a.analysis}_W{W}"); OUT.mkdir(parents=True, exist_ok=True)
S = pd.concat([pd.read_csv(f).assign(seed=int(f.name[4:6])) for f in sorted(D.glob("seed??_heldout_topm_summary.csv"))])
seeds = sorted(S.seed.unique())


def smooth(v):
    return v.rolling(W, center=True, min_periods=1).mean()


def band(v, delta):
    return set(v.index[v >= (1 - delta) * v.max()])


def segs(ks):
    ks = sorted(ks); out = []
    for k in ks:
        if out and k == out[-1][1] + 1: out[-1][1] = k
        else: out.append([k, k])
    return " ".join(f"{a}-{b}" if a != b else f"{a}" for a, b in out)


rows, bands_all = [], []
for (cond, m), g in S[S.condition.isin(["intact", "word"]) & S.top_m.isin([5, 10])].groupby(["condition", "top_m"]):
    P = g[g.K > m].pivot_table(index="K", columns="seed", values="mean").sort_index()
    for delta in DELTAS:
        full = band(smooth(P.mean(axis=1)), delta)
        bands_all.append(dict(condition=cond, top_m=m, delta=delta, fold="all 10", band=segs(full), n_K=len(full)))
        for s in seeds:
            train = smooth(P.drop(columns=s).mean(axis=1)); test = smooth(P[s])
            b = band(train, delta); bt = band(test, delta)
            inb = P.index.isin(sorted(b))
            n = len(b)
            bm = P[s][inb].mean()
            win = P[s].rolling(n, min_periods=n).mean().dropna()   # all contiguous same-width windows, seed s raw
            rows.append(dict(condition=cond, top_m=m, delta=delta, held_out_seed=s, band=segs(b), n_K=n,
                             test_peak_K=int(test.idxmax()), peak_in_band=int(test.idxmax()) in b,
                             jaccard=len(b & bt) / len(b | bt), band_minus_out=bm - P[s][~inb].mean(),
                             window_pct=float((win <= bm).mean()),
                             test_band_mean=bm))
            bands_all.append(dict(condition=cond, top_m=m, delta=delta, fold=f"without seed {s}", band=segs(b), n_K=n))
R = pd.DataFrame(rows); R.to_csv(OUT / "goodk_loso_folds.csv", index=False)
pd.DataFrame(bands_all).to_csv(OUT / "goodk_bands.csv", index=False)
summ = R.groupby(["condition", "top_m", "delta"]).agg(
    peak_in_band=("peak_in_band", "mean"), jaccard_mean=("jaccard", "mean"), jaccard_min=("jaccard", "min"),
    band_minus_out=("band_minus_out", "mean"), band_minus_out_min=("band_minus_out", "min"),
    window_pct_mean=("window_pct", "mean"), window_pct_min=("window_pct", "min"),
    n_K_min=("n_K", "min"), n_K_max=("n_K", "max")).reset_index()
summ.to_csv(OUT / "goodk_loso_summary.csv", index=False)
pd.set_option("display.width", 200)
print(f"W={W}"); print(summ.round(4).to_string(index=False)); print()
print(pd.DataFrame(bands_all).query("fold == 'all 10'").to_string(index=False))
