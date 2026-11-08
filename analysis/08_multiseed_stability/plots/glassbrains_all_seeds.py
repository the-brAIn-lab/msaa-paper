"""All seeds on one figure: Yeo-7 glass brains of the selected spatial archetypes, rows = seeds 0-9, columns = archetypes.
Matching and drawing as glassbrains_per_seed.py: Hungarian-matched to seed 0 on Pearson r of Smean, so columns keep seed-0
IDs; nodes >= 0.30 x max of the map, marker size by weight, Yeo-7 colour, views l/z/r.
Usage: python glassbrains_all_seeds.py K IDS(comma) [--nodes NODES.csv] [--out OUTDIR]
Reads every seed's compact K fit from seeds/ (see paths.py); default OUTDIR outputs/multiseed_stability/glassbrains_all_seeds.
"""
import argparse, io, re, sys
from pathlib import Path
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np, pandas as pd
from nilearn import plotting as niplot
from PIL import Image, ImageChops
from scipy.optimize import linear_sum_assignment
plt.rcParams["font.family"] = "DejaVu Sans"

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))   # analysis/08_multiseed_stability, for paths.py
import paths  # noqa: E402
ap = argparse.ArgumentParser(); ap.add_argument("K", type=int); ap.add_argument("ids")
ap.add_argument("--nodes", default=str(paths.NODES)); ap.add_argument("--out", default=str(paths.OUT / "glassbrains_all_seeds"))
a = ap.parse_args(); K = a.K; SEL = [int(x) for x in a.ids.split(",")]
NODES, OUT = Path(a.nodes), Path(a.out); OUT.mkdir(parents=True, exist_ok=True)
THR = 0.30
NETS = ["Visual", "Somatomotor", "Dorsal attention", "Ventral attention", "Limbic", "Frontoparietal", "Default mode"]
COL = {"Visual": "#D7DF23", "Somatomotor": "#39B54A", "Dorsal attention": "#00A79D", "Ventral attention": "#27AAE1",
       "Limbic": "#1C75BC", "Frontoparietal": "#92278F", "Default mode": "#EE2A7B"}
nodes = pd.read_csv(NODES).drop_duplicates("node").set_index("node").sort_index()
xyz = nodes[["x", "y", "z"]].to_numpy(float); net = nodes["network"].to_numpy()

S = {}
for s, p in paths.seed_fits("spatial", K).items():
    S[s] = np.nan_to_num(np.load(p)["Smean"].astype(float))
seeds = sorted(S); print("seeds:", seeds); assert seeds == list(range(10)) and S[0].shape == (K, len(nodes))

def zr(M): return (M - M.mean(1, keepdims=True)) / M.std(1, keepdims=True)
def mass(v): a = np.abs(v); return np.array([a[net == n].sum() for n in NETS]) / a.sum()

def brain(coords, nets, sizes, view):
    f = plt.figure(figsize=(4.0, 3.4))
    d = niplot.plot_glass_brain(None, display_mode=view, figure=f, black_bg=False, colorbar=False, plot_abs=False, annotate=True)
    for n in NETS:
        mk = nets == n
        if mk.any(): d.add_markers(marker_coords=coords[mk], marker_color=COL[n], marker_size=sizes[mk].tolist())
    buf = io.BytesIO(); f.savefig(buf, format="png", dpi=220, facecolor="white", bbox_inches="tight", pad_inches=0.02)
    plt.close(f); buf.seek(0); img = Image.open(buf).convert("RGB")
    bb = ImageChops.difference(img, Image.new("RGB", img.size, "white")).convert("L").point(lambda p: 255 if p > 8 else 0).getbbox()
    return np.asarray(img.crop(bb) if bb else img)

fig = plt.figure(figsize=(26, 17.5))
outer = fig.add_gridspec(len(seeds), len(SEL), left=0.04, right=0.995, top=0.94, bottom=0.05, hspace=0.26, wspace=0.05)
for ri, s in enumerate(seeds):
    R = zr(S[0]) @ zr(S[s]).T / S[0].shape[1]; i, j = linear_sum_assignment(-R); m = dict(zip(i, j))
    for ci, k in enumerate(SEL):
        b = m[k]; g = S[s][b]; keep = g >= THR * g.max(); v = g[keep]
        sz = np.full(len(v), 61.5) if np.isclose(v.max(), v.min()) else 18 + (v - v.min()) / (v.max() - v.min()) * 87
        inner = outer[ri, ci].subgridspec(1, 3, width_ratios=[1.15, 0.9, 1.15], wspace=0.02)
        for vi, view in enumerate("lzr"):
            ax = fig.add_subplot(inner[0, vi]); ax.imshow(brain(xyz[keep], net[keep], sz, view), interpolation="antialiased"); ax.set_axis_off()
            if vi == 1:
                ax.set_title(f"#{b} · r={R[k, b]:.2f} · {NETS[mass(g).argmax()]}", fontsize=10.5,
                             color="red" if R[k, b] < 0.8 else "0.25", pad=3)
            if vi == 0 and ci == 0:
                ax.text(-0.14, 0.5, f"seed {s}", transform=ax.transAxes, fontsize=14, fontweight="bold", ha="right", va="center")
        if ri == 0:
            pos = outer[0, ci].get_position(fig)
            fig.text((pos.x0 + pos.x1) / 2, pos.y1 + 0.022, f"A{k}", fontsize=18, fontweight="bold", ha="center")
    print("seed", s, "done", flush=True)
fig.legend([plt.Rectangle((0, 0), 1, 1, color=COL[n]) for n in NETS], NETS, loc="lower center", ncol=7, frameon=False, fontsize=14)
fig.suptitle(f"K={K} spatial archetypes across 10 seeds, Yeo-7 glass brains (seed-0 IDs; panel text: seed's own # · map r with seed 0 · top network)",
             fontsize=17, fontweight="bold", y=0.99)
for e in ("png", "pdf"): fig.savefig(OUT / f"K{K}_allseeds_glassbrains.{e}", dpi=300 if e == "png" else None, facecolor="white")
