"""Per-seed version of the (new) Fig 5 layout + cross-seed similarity of every archetype.

For each seed s = 0..9 of the 10-seed spatial-across fits at one K:
  * seed s's archetypes are matched one-to-one to seed 0 (Hungarian, maximising the summed Pearson r of the
    group-mean maps Smean = mean_i S_i over the 108 scans);
  * the 5 selected archetypes (seed-0 IDs) are drawn from seed s's matched archetype: temporal motif =
    mean +- SEM over the 108 scans of X_i C_s[:, m] (exactly the fit's sXC, rebuilt from the data and seed s's C),
    glass brains from Smean_s[m] (same 0.30 x max display threshold), Yeo-7 |weight| composition from Smean_s[m],
    heatmap = seed s's own held-out top-5 selection frequency.
Similarity (all K archetypes, seed s vs seed 0, matched pairs):
  map_r      Pearson r of Smean rows over the 700 nodes (the matching criterion)
  map_r_next best r of that seed-0 archetype with any OTHER seed-s archetype (margin = map_r - map_r_next)
  tc_r       per-scan Pearson r of X_i C_0[:,k] and X_i C_s[:,m], averaged over the 108 scans via Fisher z
             (zero lag: in spatial AA both time courses are X_i times a node-weight vector, so no lag can arise)
  tc_r_min   lowest per-scan r   (scans where the time course is constant -- archetype sits on nodes that are
             flat in that scan -- are skipped; n_scans_tc says how many were used)
  motif_r    r of the Fig-5 motif (mean over scans) between seeds
  net_r      r of the 7-network |weight| profiles; same_net = same dominant network
Usage: python fig5_per_seed.py K ARCHETYPES(comma) [--x X.npz] [--fit FULLFIT_seed0.npz] [--nodes NODES.csv]
                               [--sel selection_by_seed.csv] [--out OUTDIR]
Defaults (paths.py): X.npz from plots/loadx.py, the seed-0 fit in msaa_flexible_outputs_npz/, every seed's compact
fit from seeds/, heldout-everyk-spatial-all-results/selection_by_seed.csv, OUTDIR outputs/multiseed_stability/fig5_per_seed/K<K>.
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
for o in ("x", "fit", "nodes", "sel", "out"): ap.add_argument(f"--{o}", default=None)
a = ap.parse_args(); K = a.K; SELECTED = [int(x) for x in a.ids.split(",")]
XF = Path(a.x or paths.X_NPZ); FIT = Path(a.fit or paths.full_fit0("spatial", K)); NODES = Path(a.nodes or paths.NODES)
SELF = Path(a.sel or paths.heldout_everyk("spatial") / "selection_by_seed.csv")
OUT = Path(a.out or paths.OUT / "fig5_per_seed" / f"K{K}"); OUT.mkdir(parents=True, exist_ok=True)
CONDITIONS = ["intact", "word", "rest"]; DISPLAY_THR_FRAC = 0.30; DPI = 200
NETS = ["Visual", "Somatomotor", "Dorsal attention", "Ventral attention", "Limbic", "Frontoparietal", "Default mode"]
COL = {"Visual": "#D7DF23", "Somatomotor": "#39B54A", "Dorsal attention": "#00A79D", "Ventral attention": "#27AAE1",
       "Limbic": "#1C75BC", "Frontoparietal": "#92278F", "Default mode": "#EE2A7B"}

fit = np.load(FIT, allow_pickle=True); assert int(fit["K"]) == K
centers = np.asarray(fit["centers"], float)
nodes = pd.read_csv(NODES).drop_duplicates("node").set_index("node").sort_index()
assert np.allclose(nodes[["x", "y", "z"]].to_numpy(), centers); net = nodes["network"].to_numpy()
Xd = np.load(XF); lab = Xd["lab"]; X = [Xd[f"arr_{i}"].astype(float) for i in range(len(lab))]
assert list(lab) == list(fit["condition_labels_str"])

S, C = {}, {}
for s, p in paths.seed_fits("spatial", K).items():
    d = np.load(p)
    S[s] = np.nan_to_num(d["Smean"].astype(float)); C[s] = d["C"].astype(float)
seeds = sorted(S); assert seeds == list(range(10)), seeds
# check: seed-0 slim C reproduces the full fit's sXC
r0 = fit["results_subj"].tolist()
print("seed0 sXC max abs diff (slim C vs full fit):", max(np.abs(X[i] @ C[0] - np.asarray(r0[i]["sXC"], float)).max() for i in (0, 50, 107)))

def z(M, ax):
    sd = M.std(ax, keepdims=True)
    with np.errstate(invalid='ignore', divide='ignore'): return (M - M.mean(ax, keepdims=True)) / np.where(sd > 0, sd, np.nan)
def mass(v): a = np.abs(v); return np.array([a[net == n].sum() for n in NETS]) / a.sum()
TC = {s: np.stack([x @ C[s] for x in X]) for s in seeds}          # scans x T x K
Z0 = z(TC[0], 1)

match, rows = {}, []
for s in seeds:
    R = z(S[0], 1) @ z(S[s], 1).T / S[0].shape[1]
    i, j = linear_sum_assignment(-R); m = dict(zip(i, j)); match[s] = m
    Zs = z(TC[s], 1)
    for k in range(K):
        b = m[k]; others = np.delete(R[k], b)
        per = (Z0[:, :, k] * Zs[:, :, b]).mean(1)                   # per-scan r
        mot0, mots = TC[0][:, :, k].mean(0), TC[s][:, :, b].mean(0)
        m0, ms = mass(S[0][k]), mass(S[s][b])
        rows.append(dict(K=K, seed=s, archetype0=k, matched=int(b), map_r=R[k, b], map_r_next=others.max(),
                         margin=R[k, b] - others.max(),
                         tc_r=np.tanh(np.nanmean(np.arctanh(np.clip(per, -0.999999, 0.999999)))), tc_r_min=np.nanmin(per),
                         n_scans_tc=int(np.isfinite(per).sum()),
                         motif_r=np.corrcoef(mot0, mots)[0, 1], net_r=np.corrcoef(m0, ms)[0, 1],
                         net0=NETS[m0.argmax()], net=NETS[ms.argmax()], same_net=m0.argmax() == ms.argmax(),
                         selected=k in SELECTED))
A = pd.DataFrame(rows); A.to_csv(OUT / f"K{K}_crossseed_similarity_all.csv", index=False)

sel = pd.read_csv(SELF); sel = sel[(sel.K == K) & (sel.top_m == 5)]
# our matching must equal the one used to put the held-out selection into seed-0 IDs
chk = sel[["seed", "archetype", "archetype0"]].drop_duplicates()
assert all(match[r.seed][r.archetype0] == r.archetype for r in chk.itertuples()), "matching differs from selection file"
F = sel.groupby(["seed", "condition", "archetype0"]).selection_frequency.sum()

def sizes_of(v, lo=18, hi=105):
    if np.isclose(v.max(), v.min()): return np.full(len(v), (lo + hi) / 2)
    return lo + (v - v.min()) / (v.max() - v.min()) * (hi - lo)

def brain(coords, nets, sizes, view):
    f = plt.figure(figsize=(4.0, 3.4))
    d = niplot.plot_glass_brain(None, display_mode=view, figure=f, black_bg=False, colorbar=False, plot_abs=False, annotate=True)
    for n in NETS:
        mk = nets == n
        if mk.any(): d.add_markers(marker_coords=coords[mk], marker_color=COL[n], marker_size=sizes[mk].tolist())
    buf = io.BytesIO(); f.savefig(buf, format="png", dpi=160, facecolor="white", bbox_inches="tight", pad_inches=0.02)
    plt.close(f); buf.seek(0); img = Image.open(buf).convert("RGB")
    bb = ImageChops.difference(img, Image.new("RGB", img.size, "white")).convert("L").point(lambda p: 255 if p > 8 else 0).getbbox()
    return np.asarray(img.crop(bb) if bb else img)

for s in seeds:
    m = match[s]; a_s = A[A.seed == s].set_index("archetype0")
    fig = plt.figure(figsize=(17, 22))
    outer = fig.add_gridspec(2, 1, height_ratios=[6.4, 1.8], hspace=0.22)
    top = outer[0].subgridspec(len(SELECTED), 4, width_ratios=[1.55, 1.15, 0.90, 1.15], wspace=0.06, hspace=0.34)
    comp = []
    for ri, k in enumerate(SELECTED):
        b = m[k]; tc = TC[s][:, :, b]; mu = tc.mean(0); se = tc.std(0, ddof=1) / np.sqrt(len(tc)); t = np.arange(len(mu))
        ax = fig.add_subplot(top[ri, 0])
        if s: ax.plot(t, TC[0][:, :, k].mean(0), color="0.6", lw=1.0, label="seed 0")
        ax.plot(t, mu, lw=1.8, color="C0", label=f"seed {s}"); ax.fill_between(t, mu - se, mu + se, alpha=0.2, color="C0")
        ax.axhline(0, lw=0.7, ls=":"); ax.set_xlim(0, len(mu) - 1); ax.tick_params(labelsize=12)
        ax.set_ylabel(f"A{k}", rotation=0, labelpad=28, fontsize=14, fontweight="bold", va="center")
        r = a_s.loc[k]; bad = r.map_r < 0.8 or r.net_r < 0.9 or r.tc_r < 0.8
        ax.text(0.01, 0.97, f"seed-{s} #{b}   map r={r.map_r:.2f} (next {r.map_r_next:.2f})   time r={r.tc_r:.2f}" + ("" if r.same_net else f"   top net {r.net0[:7]}->{r.net[:7]}"),
                transform=ax.transAxes, va="top", fontsize=10, color="red" if bad else "0.25",
                fontweight="bold" if bad else "normal")
        if ri == len(SELECTED) - 1: ax.set_xlabel("Timepoint", fontsize=14)
        else: ax.set_xticklabels([])
        if ri == 0:
            ax.set_title("Temporal motif", fontweight="bold", fontsize=15, pad=10)
            if s: ax.legend(frameon=False, fontsize=9, loc="lower right")
        g = S[s][b]; keep = g >= DISPLAY_THR_FRAC * g.max(); sz = sizes_of(g[keep])
        bax = []
        for ci, v in enumerate("lzr", start=1):
            axb = fig.add_subplot(top[ri, ci]); bax.append(axb)
            axb.imshow(brain(centers[keep], net[keep], sz, v), interpolation="antialiased"); axb.set_axis_off()
        if ri == 0: bax[1].set_title("Spatial deployment", fontweight="bold", fontsize=15, pad=10)
        comp.append(mass(g))
    bottom = outer[1].subgridspec(1, 2, width_ratios=[1.70, 1.00], wspace=0.18)
    ax = fig.add_subplot(bottom[0, 0]); y = np.arange(len(SELECTED)); left = np.zeros(len(SELECTED)); comp = np.array(comp)
    for j, n in enumerate(NETS):
        ax.barh(y, comp[:, j], left=left, color=COL[n], height=0.72); left += comp[:, j]
    ax.set_yticks(y); ax.set_yticklabels([f"A{k}" for k in SELECTED], fontsize=15); ax.invert_yaxis(); ax.set_xlim(0, 1)
    ax.set_xlabel("Fraction of |spatial weight|", fontsize=15); ax.tick_params(axis="x", labelsize=15)
    ax.set_title(f"Yeo-7 composition (seed {s})", fontweight="bold", fontsize=17)
    ax = fig.add_subplot(bottom[0, 1])
    H = np.array([[F.get((s, c, k), 0.0) for c in CONDITIONS] for k in SELECTED])
    im = ax.imshow(H, aspect="auto", vmin=0, vmax=1)
    ax.set_xticks(range(3)); ax.set_xticklabels(["Intact", "Word", "Rest"], fontsize=15)
    ax.set_yticks(y); ax.set_yticklabels([f"A{k}" for k in SELECTED], fontsize=15)
    ax.set_title(f"Held-out top-5 selection (seed {s})", fontweight="bold", fontsize=17)
    for r_ in range(H.shape[0]):
        for c_ in range(3): ax.text(c_, r_, f"{H[r_, c_]:.2f}", ha="center", va="center", fontsize=13, color="white" if H[r_, c_] >= 0.55 else "black")
    cb = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04); cb.set_label("Selection frequency", fontsize=15)
    fig.legend([plt.Rectangle((0, 0), 1, 1, color=COL[n]) for n in NETS], NETS, loc="lower center", ncol=4,
               frameon=False, fontsize=14, bbox_to_anchor=(0.5, 0.018))
    fig.suptitle(f"K={K} spatial archetypes, seed {s} (matched to seed-0 IDs)", fontsize=19, fontweight="bold", y=0.995)
    fig.subplots_adjust(bottom=0.115, top=0.965)
    fig.savefig(OUT / f"K{K}_seed{s}_networks.png", dpi=DPI, bbox_inches="tight"); fig.savefig(OUT / f"K{K}_seed{s}_networks.pdf", bbox_inches="tight")
    plt.close(fig); print("seed", s, "done", flush=True)

pd.set_option("display.width", 220)
A1 = A[A.seed > 0]
summ = A1.groupby("archetype0").agg(map_r_med=("map_r", "median"), map_r_min=("map_r", "min"), margin_min=("margin", "min"),
                                    tc_r_med=("tc_r", "median"), tc_r_min=("tc_r", "min"), motif_r_min=("motif_r", "min"),
                                    net_r_min=("net_r", "min"), same_net=("same_net", "sum"), net0=("net0", "first"))
summ.to_csv(OUT / f"K{K}_crossseed_similarity_by_archetype.csv")
print(summ.loc[SELECTED].round(3).to_string())
print("ALL archetypes, seeds 1-9: map r median %.3f, tc r median %.3f, Spearman(map_r, tc_r) %.2f, frac map_r>=0.9 %.2f, tc_r>=0.9 %.2f"
      % (A1.map_r.median(), A1.tc_r.median(), A1[["map_r", "tc_r"]].corr("spearman").iloc[0, 1], (A1.map_r >= .9).mean(), (A1.tc_r >= .9).mean()))
