"""All seeds on one figure: Yeo-7 composition bars of the selected spatial archetypes, one panel per archetype, one bar per seed.
Matching as glassbrains_per_seed.py (Hungarian to seed 0 on Pearson r of Smean, so panels keep seed-0 IDs). Composition as paper
Fig 5: |group map| summed per Yeo-7 network / total, over all nodes (no display threshold).
Usage: python yeo7_bars_all_seeds.py K IDS(comma) [--nodes NODES.csv] [--out OUTDIR]
Reads every seed's compact K fit from seeds/ (see paths.py); default OUTDIR outputs/multiseed_stability/yeo7_bars_all_seeds.
"""
import argparse, re, sys
from pathlib import Path
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np, pandas as pd
from scipy.optimize import linear_sum_assignment
plt.rcParams["font.family"] = "DejaVu Sans"

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))   # analysis/08_multiseed_stability, for paths.py
import paths  # noqa: E402
ap = argparse.ArgumentParser(); ap.add_argument("K", type=int); ap.add_argument("ids")
ap.add_argument("--nodes", default=str(paths.NODES)); ap.add_argument("--out", default=str(paths.OUT / "yeo7_bars_all_seeds"))
a = ap.parse_args(); K = a.K; SEL = [int(x) for x in a.ids.split(",")]
NODES, OUT = Path(a.nodes), Path(a.out); OUT.mkdir(parents=True, exist_ok=True)
NETS = ["Visual", "Somatomotor", "Dorsal attention", "Ventral attention", "Limbic", "Frontoparietal", "Default mode"]
COL = {"Visual": "#D7DF23", "Somatomotor": "#39B54A", "Dorsal attention": "#00A79D", "Ventral attention": "#27AAE1",
       "Limbic": "#1C75BC", "Frontoparietal": "#92278F", "Default mode": "#EE2A7B"}
net = pd.read_csv(NODES).drop_duplicates("node").set_index("node").sort_index()["network"].to_numpy()

S = {}
for s, p in paths.seed_fits("spatial", K).items():
    S[s] = np.nan_to_num(np.load(p)["Smean"].astype(float))
seeds = sorted(S); print("seeds:", seeds); assert seeds == list(range(10)) and S[0].shape == (K, len(net))

def zr(M): return (M - M.mean(1, keepdims=True)) / M.std(1, keepdims=True)
def mass(v): a = np.abs(v); return np.array([a[net == n].sum() for n in NETS]) / a.sum()

rows = []
for s in seeds:
    R = zr(S[0]) @ zr(S[s]).T / S[0].shape[1]; i, j = linear_sum_assignment(-R); m = dict(zip(i, j))
    for k in SEL:
        b = m[k]; ms = mass(S[s][b])
        rows.append(dict(seed=s, archetype0=k, matched=int(b), map_r=R[k, b], top_net=NETS[ms.argmax()], **dict(zip(NETS, ms))))
P = pd.DataFrame(rows); P.to_csv(OUT / f"K{K}_allseeds_yeo7_bars.csv", index=False)

fig, axs = plt.subplots(1, len(SEL), figsize=(15, 4.6), sharey=True)
for ax, k in zip(axs, SEL):
    d = P[P.archetype0 == k].sort_values("seed"); left = np.zeros(len(d))
    for n in NETS:
        ax.barh(d.seed, d[n], left=left, color=COL[n], height=0.8, edgecolor="white", lw=0.5); left += d[n].to_numpy()
    for y, r in zip(d.seed, d.map_r):
        ax.text(1.02, y, f"{r:.2f}", va="center", fontsize=9, color="red" if r < 0.8 else "0.3")
    ax.text(1.02, -0.85, "map r", va="center", fontsize=9, color="0.3")
    ax.set_title(f"A{k}", fontsize=14, fontweight="bold", pad=14); ax.set_xlim(0, 1); ax.set_xticks([0, 0.5, 1])
    ax.set_xlabel("Fraction of |spatial weight|", fontsize=10); ax.tick_params(labelsize=10, length=0)
    ax.spines[["top", "right", "left"]].set_visible(False)
axs[0].set_yticks(seeds); axs[0].set_yticklabels([f"seed {s}" for s in seeds]); axs[0].invert_yaxis()
fig.legend([plt.Rectangle((0, 0), 1, 1, color=COL[n]) for n in NETS], NETS, loc="lower center", ncol=7, frameon=False, fontsize=11)
fig.suptitle(f"K={K} spatial archetypes: Yeo-7 composition across 10 seeds (seed-0 IDs)", fontsize=15, fontweight="bold")
fig.tight_layout(rect=(0, 0.07, 1, 0.98), w_pad=2.6)
for e in ("png", "pdf"): fig.savefig(OUT / f"K{K}_allseeds_yeo7_bars.{e}", dpi=300 if e == "png" else None, bbox_inches="tight")
pd.set_option("display.width", 200)
print(P.groupby("archetype0")[NETS].agg(["min", "max"]).round(3).loc[SEL].to_string())
print(P.pivot(index="archetype0", columns="seed", values="top_net").loc[SEL].to_string())
