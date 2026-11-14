"""Top-5 vs top-10 held-out selection stability, spatial K neighbourhood.

For each spatial K and m in {5, 10}, using each seed's OWN participant-held-out ranking (40 splits = 20 outer reps x 2
directions; each split ranks on 18 participants, is scored on the other 18):
  seed top-m set   = the m archetypes that seed's 40 splits select most often at top_m=m (ties: mean rank position),
                     mapped to that K's seed-0 IDs by Hungarian matching on Pearson r of Smean maps (700 nodes).
  jaccard_seedsets_*          45 seed pairs of seed top-m sets
  overlap_mean                mean |A & B| over those pairs (motifs in common)
  jaccard_same_split_*        per-split sets, identical participant split, seed pairs (10 C 2 x 40)
  jaccard_across_splits_*     per-split sets, different splits, same seed
  jaccard_random              E[J] of two random m-subsets of K (exact, hypergeometric); overlap_random = m^2/K
  *_cc                        chance-corrected Jaccard (J - J_random) / (1 - J_random)
  core_10of10 / n_core        seed-0 IDs in all 10 seed top-m sets
  yeo7_r_*                    pairwise r of the aggregate (mean) Yeo-7 |weight| share of each seed's top-m set
  selected_not_recovered      count of (seed 1-9, pick) occurrences in seed top-m sets whose Hungarian match to the
                              seed-0 archetype has map r < 0.9 (max 9*m); frac_ = / (9*m)
  RECOVERY-AWARE (main): a pick whose Hungarian match to seed 0 has r < 0.9 gets a unique seed-specific label (cannot count as
  overlap); *_forced columns use the forced Hungarian label (upper bound). Yeo-7 uses the actual selected maps.
  M0-M4 (K=46 seed-0 #28,#35,#16,#39,#45): forced label in each (K, seed) = Hungarian of all 46 K=46 seed-0 maps vs that
  fit's maps; recovered = r >= 0.9 (M0M4_selection_by_K.csv: recovered / selected&recovered / selected-forced counts).

Usage: python topm_stability_neighbourhood.py K[,K...] [--nodes NODES.csv] [--heldout DIR] [--out OUTDIR]
Defaults (paths.py): every seed's compact spatial fits from seeds/ (K=46 always needed for M0-M4),
heldout-everyk-spatial-all-results/, OUTDIR outputs/multiseed_stability/topm_stability_neighbourhood.
"""
import argparse
import os
import re
import sys
from itertools import combinations
from math import comb
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import matplotlib.colorbar as _mcb
_mcb.Colorbar.n_rasterize = 10**9  # keep colorbars vector in PDFs
_add_solids = _mcb.Colorbar._add_solids


def _add_solids_no_seams(self, *a, **k):
    # stroke each colour band in its own colour: no hairline gaps in PDF viewers (matplotlib colorbar docs)
    _add_solids(self, *a, **k)
    if getattr(self, "solids", None) is not None:
        self.solids.set_edgecolor("face")
        self.solids.set_linewidth(0.3)


_mcb.Colorbar._add_solids = _add_solids_no_seams


def vheat(ax, arr, origin="upper", **kw):
    """imshow look-alike drawn as vector cells (pcolormesh), so the PDF holds no raster image."""
    a = np.ma.masked_invalid(np.asarray(arr, dtype=float))
    n, m = a.shape
    im = ax.pcolormesh(np.arange(m + 1) - 0.5, np.arange(n + 1) - 0.5, a, shading="flat",
                       edgecolors="face", linewidth=0.3, **kw)  # face-coloured edges: no seams
    ax.set_xlim(-0.5, m - 0.5)
    ax.set_ylim((-0.5, n - 0.5) if origin == "lower" else (n - 0.5, -0.5))
    return im


try:
    from scipy.optimize import linear_sum_assignment
except ImportError:  # device VM has no scipy
    def linear_sum_assignment(cost):
        cost = np.asarray(cost, float); n = cost.shape[0]
        u = np.zeros(n + 1); v = np.zeros(n + 1); p = np.zeros(n + 1, int); way = np.zeros(n + 1, int)
        for i in range(1, n + 1):
            p[0] = i; j0 = 0; minv = np.full(n + 1, np.inf); used = np.zeros(n + 1, bool)
            while True:
                used[j0] = True; i0 = p[j0]
                cur = cost[i0 - 1] - u[i0] - v[1:]
                upd = ~used[1:] & (cur < minv[1:]); minv[1:][upd] = cur[upd]; way[1:][upd] = j0
                free = np.where(~used[1:])[0] + 1; j1 = free[np.argmin(minv[free])]; delta = minv[j1]
                u[p[used]] += delta; v[used] -= delta; minv[~used] -= delta; j0 = j1
                if p[j0] == 0: break
            while True:
                j1 = way[j0]; p[j0] = p[j1]; j0 = j1
                if j0 == 0: break
        rows = p[1:] - 1; order = np.argsort(rows); return rows[order], np.arange(n)[order]

plt.rcParams["font.family"] = "DejaVu Sans"
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))   # analysis/08_multiseed_stability, for paths.py
import paths  # noqa: E402
ap = argparse.ArgumentParser(); ap.add_argument("ks"); ap.add_argument("--nodes", default=str(paths.NODES))
ap.add_argument("--heldout", default=str(paths.heldout_everyk("spatial")))
ap.add_argument("--out", default=str(paths.OUT / "topm_stability_neighbourhood"))
ap.add_argument("--plot-only", action="store_true", help="redraw the figures from the CSVs already in OUTDIR")
a = ap.parse_args(); OUT = Path(a.out); OUT.mkdir(parents=True, exist_ok=True)
KS = [int(k) for k in a.ks.split(",")]
SEEDS = list(range(10)); CONDS = ["intact", "word", "rest"]; MS = [5, 10]; REP_R = 0.9
NETS = ["Visual", "Somatomotor", "Dorsal attention", "Ventral attention", "Limbic", "Frontoparietal", "Default mode"]
NSHORT = ["Vis", "SomMot", "DorsAttn", "VentAttn", "Limbic", "FP", "DMN"]
MOT = {"M0": 28, "M1": 35, "M2": 16, "M3": 39, "M4": 45}       # K=46 seed-0 ids
LEG = {"M0": "A4", "M1": "A5", "M2": "A7", "M3": "A8", "M4": "A9"}
nodes = pd.read_csv(a.nodes)
nodes = nodes.drop_duplicates("node").set_index("node").sort_index()["network"].to_numpy()
HE = Path(a.heldout)
# --plot-only only redraws: it must not truncate the summary written by the full run
log = open(os.devnull if a.plot_only else OUT / "summary.txt", "w")
def say(*a):
    s = " ".join(str(x) for x in a); print(s); log.write(s + "\n"); log.flush()


def rowcorr(A, B):
    A = A - A.mean(1, keepdims=True); B = B - B.mean(1, keepdims=True)
    return (A / np.linalg.norm(A, axis=1)[:, None]) @ (B / np.linalg.norm(B, axis=1)[:, None]).T


def mass(v):
    a = np.abs(np.nan_to_num(v)); return np.array([a[nodes == n].sum() for n in NETS]) / a.sum()


def jac(a, b): a, b = set(a), set(b); return len(a & b) / len(a | b)


def jrand(m, K):  # exact E[J] for two independent uniform m-subsets of K
    return sum(comb(m, i) * comb(K - m, m - i) / comb(K, m) * i / (2 * m - i) for i in range(m + 1))


def maps(K):
    M = {s: np.load(p)["Smean"].astype(float) for s, p in paths.seed_fits("spatial", K).items()}
    assert sorted(M) == SEEDS, (K, sorted(M)); return M


if a.plot_only:
    S = pd.read_csv(OUT / "topm_stability_neighbourhood.csv"); MT = pd.read_csv(OUT / "M0M4_selection_by_K.csv")
    S = S[S.K.isin(KS)]
else:
    # held-out per-split runs for these K, top_m 5 and 10 (one pass over the big file)
    cache = OUT / f".runs_K{'_'.join(map(str, KS))}.csv.gz"
    if cache.exists():
        RUNS = pd.read_csv(cache)
    else:
        parts = [ch[ch.K.isin(KS) & ch.top_m.isin(MS)] for ch in pd.read_csv(
            HE / "heldout_topm_runs_all_seeds.csv.gz", chunksize=500000,
            usecols=["top_m", "selected_archetypes", "condition", "K", "outer_rep", "direction", "fit_seed"])]
        RUNS = pd.concat(parts); RUNS.to_csv(cache, index=False)
    FREQ = pd.concat([pd.read_csv(HE / f"seed{s:02d}_heldout_archetype_selection_stability.csv").assign(seed=s) for s in SEEDS])
    FREQ = FREQ[FREQ.K.isin(KS) & FREQ.top_m.isin(MS)]

    M46 = maps(46); MOTIDX = list(MOT.values())
    rows, mrows, sets_rows, motrows = [], [], [], []
    for K in KS:
        M = maps(K); match, R, inv, rec, mfor = {}, {}, {}, {}, {}
        for s in SEEDS:
            R[s] = rowcorr(M[0], M[s]); i, j = linear_sum_assignment(-R[s]); match[s] = dict(zip(i.tolist(), j.tolist()))
            inv[s] = {b: a for a, b in match[s].items()}
            rec[s] = {a: (s == 0) or R[s][a, b] >= REP_R for a, b in match[s].items()}
            # forced label of M0-M4 in (K, seed): Hungarian of all 46 K=46 seed-0 maps vs this fit's K maps (padded square)
            rm = rowcorr(M46[0], M[s]); n = max(rm.shape); P = np.full((n, n), -1.0); P[:rm.shape[0], :rm.shape[1]] = rm
            i, j = linear_sum_assignment(-P); hm = dict(zip(i.tolist(), j.tolist()))
            mfor[s] = [(hm[a], float(rm[a, hm[a]]) if hm[a] < rm.shape[1] else np.nan) for a in MOTIDX]
        def lab(s, own_id):  # recovery-aware label: seed-0 id if recovered (r >= 0.9), else unique seed-specific label
            a = inv[s][own_id]; return a if rec[s][a] else f"u{s}_{own_id}"
        for m in MS:
            jr = jrand(m, K)
            for c in CONDS:
                f = FREQ[(FREQ.K == K) & (FREQ.top_m == m) & (FREQ.condition == c)].sort_values(
                    ["seed", "selection_frequency", "mean_rank_position"], ascending=[True, False, True])
                own = {s: f[f.seed == s].head(m).archetype.astype(int).tolist() for s in SEEDS}
                topF = {s: [inv[s][a] for a in own[s]] for s in SEEDS}          # forced Hungarian labels
                topR = {s: [lab(s, a) for a in own[s]] for s in SEEDS}          # recovery-aware labels
                rc = RUNS[(RUNS.K == K) & (RUNS.top_m == m) & (RUNS.condition == c)]
                spl = {s: g.set_index(["outer_rep", "direction"]).selected_archetypes.map(lambda t: [int(x) for x in re.findall(r"\d+", t)])
                       for s, g in rc.groupby("fit_seed")}
                assert sorted(spl) == SEEDS and all(len(g) == 40 for g in spl.values()), (K, m, c)
                keys = spl[0].index
                row = dict(K=K, top_m=m, condition=c, seed0_set=sorted(topF[0]))
                for tag, top, L in [("", topR, lab), ("_forced", topF, lambda s, b: inv[s][b])]:
                    pj = np.array([jac(top[a], top[b]) for a, b in combinations(SEEDS, 2)])
                    ov = np.array([len(set(top[a]) & set(top[b])) for a, b in combinations(SEEDS, 2)])
                    same = np.array([jac([L(a, x) for x in spl[a][k]], [L(b, x) for x in spl[b][k]]) for k in keys for a, b in combinations(SEEDS, 2)])
                    # hierarchical: seed pairs within (split, direction) -> directions within outer split -> outer splits
                    Sm = same.reshape(len(keys), -1).mean(1)
                    same_h = pd.Series(Sm, index=keys).groupby(level="outer_rep").mean().mean()
                    assert abs(same_h - same.mean()) < 1e-12, (K, m, c, tag)
                    cnt = pd.Series([a for s in SEEDS for a in top[s]]).value_counts(); core = sorted(cnt[cnt == 10].index.tolist())
                    row.update({f"jaccard_seedsets_mean{tag}": pj.mean(), f"jaccard_seedsets_median{tag}": np.median(pj),
                                f"jaccard_seedsets_min{tag}": pj.min(), f"jaccard_seedsets_mean_cc{tag}": (pj.mean() - jr) / (1 - jr),
                                f"overlap_mean{tag}": ov.mean(), f"frac_pairs_identical{tag}": np.mean(pj == 1),
                                f"jaccard_same_split_across_seeds{tag}": same.mean(), f"jaccard_same_split_across_seeds_hier{tag}": same_h,
                                f"jaccard_same_split_across_seeds_cc{tag}": (same.mean() - jr) / (1 - jr),
                                f"n_core_10of10{tag}": len(core), f"core_10of10{tag}": core, f"n_selected_ever{tag}": len(cnt)})
                # within-seed labels are consistent under either labelling, so this is shared
                within = np.array([jac(g.iloc[i], g.iloc[j]) for g in spl.values() for i, j in combinations(range(40), 2)])
                within_h = np.mean([np.mean([jac(g.iloc[i], g.iloc[j]) for i, j in combinations(range(40), 2)]) for g in spl.values()])
                assert abs(within_h - within.mean()) < 1e-12 and all(g.index.equals(keys) for g in spl.values())
                unrec = sum(isinstance(x, str) for s in SEEDS for x in topR[s])
                comp = np.array([np.mean([mass(M[s][b]) for b in own[s]], axis=0) for s in SEEDS])   # actual selected maps
                py = np.corrcoef(comp)[np.triu_indices(10, 1)]
                row.update(jaccard_across_splits_within_seed=within.mean(), jaccard_across_splits_within_seed_hier=within_h, jaccard_across_splits_within_seed_cc=(within.mean() - jr) / (1 - jr),
                           jaccard_random=jr, overlap_random=m * m / K, selected_not_recovered=unrec, frac_selected_not_recovered=unrec / (9 * m),
                           yeo7_r_mean=py.mean(), yeo7_r_min=py.min())
                rows.append(row)
                for s in SEEDS:
                    sets_rows.append(dict(K=K, top_m=m, condition=c, seed=s, set_own_ids=sorted(own[s]), set_forced_seed0_ids=sorted(topF[s]),
                                          set_recovery_aware=sorted(map(str, topR[s])), n_not_recovered=sum(isinstance(x, str) for x in topR[s])))
                # M0-M4 per-motif selection, split into genuine recoveries and forced (r < 0.9) labels
                for mi, mlab in enumerate(MOT):
                    ok = [s for s in SEEDS if mfor[s][mi][1] >= REP_R]
                    sel = [s for s in SEEDS if mfor[s][mi][0] in own[s]]
                    fq = f.set_index(["seed", "archetype"]).selection_frequency
                    sf = [fq.get((s, mfor[s][mi][0]), 0.0) for s in ok]
                    motrows.append(dict(K=K, condition=c, top_m=m, motif=mlab, legacy=LEG[mlab], k46_seed0_id=MOT[mlab],
                                        n_recovered=len(ok), n_selected_and_recovered=len(set(sel) & set(ok)),
                                        n_selected_forced_r_lt_090=len(set(sel) - set(ok)),
                                        frac_selected_given_recovered=len(set(sel) & set(ok)) / len(ok) if ok else np.nan,
                                        mean_split_selection_freq_given_recovered=np.mean(sf) if ok else np.nan,
                                        forced_r_by_seed=[round(mfor[s][mi][1], 3) for s in SEEDS]))
        for i, mlab in enumerate(MOT):
            rs = np.array([mfor[s][i][1] for s in SEEDS])
            mrows.append(dict(K=K, motif=mlab, legacy=LEG[mlab], k46_seed0_id=MOT[mlab], forced_match_seed0=mfor[0][i][0],
                              r_seed0=rs[0], r_median=np.median(rs), r_min=rs.min(), seeds_recovered=int((rs >= REP_R).sum()),
                              r_by_seed=[round(x, 3) for x in rs]))

    S = pd.DataFrame(rows); MR = pd.DataFrame(mrows); MT = pd.DataFrame(motrows)
    S["frac_M0M4_recovered"] = S.K.map(MR.groupby("K").seeds_recovered.sum() / 50)
    S.to_csv(OUT / "topm_stability_neighbourhood.csv", index=False); MR.to_csv(OUT / "M0M4_recovery_by_K.csv", index=False)
    MT.to_csv(OUT / "M0M4_selection_by_K.csv", index=False); pd.DataFrame(sets_rows).to_csv(OUT / "topm_sets_by_seed.csv", index=False)
    pd.set_option("display.width", 250); pd.set_option("display.max_columns", 40)
    cols = ["K", "top_m", "condition", "jaccard_seedsets_mean", "jaccard_seedsets_mean_forced", "jaccard_seedsets_mean_cc", "overlap_mean",
            "jaccard_same_split_across_seeds", "jaccard_same_split_across_seeds_forced", "jaccard_across_splits_within_seed", "jaccard_random",
            "n_core_10of10", "n_core_10of10_forced", "selected_not_recovered", "yeo7_r_mean", "yeo7_r_min", "frac_M0M4_recovered"]
    say(S[cols].round(3).to_string(index=False))
    say(""); say(MT.drop(columns="forced_r_by_seed").round(3).to_string(index=False))
    say(""); say(MR.round(3).to_string(index=False))


# ---------------------------------------------------------------- S05: one figure per m, the five prioritised metrics
# Every metric is on 0..1, so one colour scale serves all rows; the full set of diagnostics stays in the CSV.
# Drawn at the journal's print width (Network Neuroscience reproduces figures ~5 in across), so the point sizes are
# the printed sizes: one panel per condition, stacked, metrics as rows and K as columns.
FS_TITLE, FS_LABEL, FS_TICK, FS_ANNOT = 7.5, 7, 6.5, 6
plt.rcParams.update({"axes.linewidth": 0.6, "xtick.major.width": 0.6, "ytick.major.width": 0.6,
                     "xtick.major.size": 2, "ytick.major.size": 2, "xtick.major.pad": 1.5, "ytick.major.pad": 1.5})
CLAB = {"intact": "Intact narrative", "word": "Word-scrambled", "rest": "Rest"}
mets = [("jaccard_seedsets_mean", "Recovery-aware Jaccard"), ("jaccard_seedsets_mean_cc", "Chance-corrected Jaccard"),
        ("jaccard_same_split_across_seeds", "Same split, across seeds"),
        ("jaccard_across_splits_within_seed", "Across splits, within seed"), ("frac_M0M4_recovered", "M0\u2013M4 recovered")]
for m in MS:
    fig, axr = plt.subplots(3, 1, figsize=(5.0, 0.75 + 0.135 * 3 * (len(mets) + 1.6)), sharex=True, layout="constrained")
    fig.get_layout_engine().set(w_pad=0.02, h_pad=0.02, hspace=0.03)
    for ax, c in zip(axr, CONDS):
        d = S[(S.top_m == m) & (S.condition == c)].set_index("K").loc[KS]
        X = d[[k for k, _ in mets]].to_numpy(float).T          # metrics x K
        im = vheat(ax, X, vmin=0, vmax=1, cmap="viridis")
        for i in range(X.shape[0]):
            for j in range(X.shape[1]):
                ax.text(j, i, f"{X[i, j]:.2f}".replace("0.", ".", 1), ha="center", va="center", fontsize=FS_ANNOT,
                        color="white" if X[i, j] < 0.55 else "black")
        ax.set_yticks(range(len(mets))); ax.set_yticklabels([t for _, t in mets], fontsize=FS_TICK)
        ax.set_xticks(range(len(KS))); ax.set_xticklabels(KS, fontsize=FS_TICK)
        ax.tick_params(length=0)
        ax.set_title(CLAB[c], fontsize=FS_TITLE, fontweight="bold", pad=2)
    axr[-1].set_xlabel("Spatial model order (K)", fontsize=FS_LABEL, labelpad=1)
    cb = fig.colorbar(im, ax=axr, shrink=0.5, aspect=20, pad=0.02, ticks=[0, 0.5, 1])
    cb.set_label("Similarity / fraction", fontsize=FS_LABEL, labelpad=2)
    cb.ax.tick_params(labelsize=FS_TICK, length=2, width=0.6); cb.outline.set_linewidth(0.6)
    fig.suptitle(f"Stability of held-out top-{m} selection across fit seeds", fontsize=9, fontweight="bold")
    for e in ("png", "pdf"):
        fig.savefig(OUT / f"topm_stability_neighbourhood_top{m}.{e}", dpi=600, bbox_inches="tight", pad_inches=0.02)
    plt.close(fig)

# S06: M0-M4 per-motif selection at K=46: solid = selected & recovered (r >= .9), hatched = forced label selected (r < .9)
if 46 in KS:
    from matplotlib.patches import Patch
    from matplotlib.lines import Line2D
    COL = {5: "#4C9BE8", 10: "#E4572E"}
    plt.rcParams["hatch.linewidth"] = 0.6
    fig, axs = plt.subplots(3, 1, figsize=(5.0, 3.9), sharex=True, layout="constrained")
    fig.get_layout_engine().set(w_pad=0.02, h_pad=0.02, hspace=0.04)
    for ci, c in enumerate(CONDS):
        ax = axs[ci]; d = MT[(MT.K == 46) & (MT.condition == c)]; x = np.arange(5)
        nrec = d[d.top_m == 5].set_index("motif").loc[list(MOT)].n_recovered.to_numpy()
        for o, m in [(-0.2, 5), (0.2, 10)]:
            col = COL[m]
            dm = d[d.top_m == m].set_index("motif").loc[list(MOT)]
            a_, b = dm.n_selected_and_recovered.to_numpy(), dm.n_selected_forced_r_lt_090.to_numpy()
            ax.bar(x + o, a_, 0.4, color=col)
            ax.bar(x + o, b, 0.4, bottom=a_, color="white", edgecolor=col, hatch="////", lw=0.6)
            for xi, ai, bi, nr in zip(x + o, a_, b, nrec):   # label sits above the bar and the dashed line
                ax.text(xi, max(ai + bi, nr) + 0.25, f"{ai}" + (f"+{bi}" if bi else ""), ha="center", va="bottom", fontsize=FS_ANNOT)
        for xi, n in zip(x, nrec): ax.plot([xi - 0.44, xi + 0.44], [n, n], color="k", lw=0.8, ls="--")
        ax.set_xticks(x); ax.set_xticklabels(list(MOT), fontsize=FS_TICK)
        ax.set_yticks([0, 5, 10]); ax.tick_params(labelsize=FS_TICK)
        ax.set_ylim(0, 12.2); ax.set_title(CLAB[c], fontsize=FS_TITLE, fontweight="bold", pad=2, loc="left")
        ax.spines[["top", "right"]].set_visible(False)
    for ax in axs: ax.set_ylabel("Fit seeds (of 10)", fontsize=FS_LABEL, labelpad=1)
    h = [Patch(color=COL[5], label="Top-5"), Patch(color=COL[10], label="Top-10"),
         Patch(facecolor="white", edgecolor="0.4", hatch="////", lw=0.6, label="Selected, not recovered (r < 0.90)"),
         Line2D([], [], color="k", ls="--", lw=0.8, label="Seeds with motif recovered")]
    fig.suptitle("K = 46: M0\u2013M4 in each seed's own top-5 and top-10", fontsize=9, fontweight="bold")
    fig.legend(handles=h, frameon=False, fontsize=FS_TICK, loc="outside lower center", ncol=2, handlelength=1.6,
               handleheight=0.9, columnspacing=1.2)
    for e in ("png", "pdf"): fig.savefig(OUT / f"M0M4_top5_vs_top10_K46.{e}", dpi=600, bbox_inches="tight", pad_inches=0.02)
    plt.close(fig)
print("figures:", ", ".join(sorted(p.name for p in OUT.glob("*.png"))))
