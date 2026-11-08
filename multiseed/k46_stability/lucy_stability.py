"""Lucy's 2026-11-07 asks, from existing 10-seed outputs (no refits, no redecode).

1. K=46 alignment for the five paper archetypes (legacy A4,A5,A6,A12,A8 = seed-0/seqnum 28,35,37,31,39).
2. Stability of each seed's OWN participant-held-out top-5 selection (K=46 and K=51 spatial, K=14 temporal):
   selection frequency per seed-0 archetype, pairwise Jaccard of the 10 top-5 sets (by condition), Jaccard on the
   identical participant splits across seeds vs. across splits within a seed, and r of aggregate Yeo-7 composition.
3. Temporal K=14: is the A5+A9 subspace / reconstruction preserved in seeds 4 and 9 where they merge?

Matching: Hungarian on Pearson r, each seed to seed 0 (spatial: Smean maps over 700 nodes; temporal: C columns,
zero lag) -- the same matching used for selection_by_seed.csv and the per-seed network figures.
"Seed top-5 set" = the 5 archetypes that seed's 40 held-out splits (20 reps x 2 directions) select most often at
top_m=5 (ties broken by mean rank position). Each split ranks on 18 participants and is scored on the other 18.

Usage: python lucy_stability.py SEEDS_RESULTS_DIR AA_DIR OUTDIR
"""
import re
import sys
from itertools import combinations
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

try:
    from scipy.optimize import linear_sum_assignment
except ImportError:  # device VM has no scipy: O(n^3) Hungarian, min-cost, square
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
SR, AA, OUT = map(Path, sys.argv[1:4]); OUT.mkdir(parents=True, exist_ok=True)
SEEDS = list(range(10)); CONDS = ["intact", "word", "rest"]
NETS = ["Visual", "Somatomotor", "Dorsal attention", "Ventral attention", "Limbic", "Frontoparietal", "Default mode"]
NSHORT = ["Vis", "SomMot", "DorsAttn", "VentAttn", "Limbic", "FP", "DMN"]
NCOL = ["#D7DF23", "#39B54A", "#00A79D", "#27AAE1", "#1C75BC", "#92278F", "#EE2A7B"]
PAPER = {28: "A4", 35: "A5", 37: "A6", 31: "A12", 39: "A8"}   # seqnum id -> legacy paper label (legacy_full fit, map r = 1.0)
LEGACY = {16: "A7", 45: "A9", **PAPER}
import os
if os.environ.get("PAPER_SET"):   # e.g. PAPER_SET=28,35,37,16,39 (ben 2026-11-07: A12 -> #16 = legacy A7)
    PAPER = {int(a): LEGACY.get(int(a), f"#{a}") for a in os.environ["PAPER_SET"].split(",")}
REP_R = 0.9   # matched r >= 0.9 = "same archetype recovered" (matched r is bimodal)
nodes = pd.read_csv(AA / "msaa-paper/reference_outputs/legacy_full/paper_figure_outputs/Figure5_spatial_maps_allnodes.csv")
nodes = nodes.drop_duplicates("node").set_index("node").sort_index()["network"].to_numpy()
HE = SR / "heldout-everyk"
log = open(OUT / "summary.txt", "w")
def say(*a):
    s = " ".join(str(x) for x in a); print(s); log.write(s + "\n"); log.flush()


def seed_of(p):
    m = re.search(r"seed(\d\d)", str(p)); return int(m.group(1))


def rowcorr(A, B):  # rows of A vs rows of B
    A = A - A.mean(1, keepdims=True); B = B - B.mean(1, keepdims=True)
    na, nb = np.linalg.norm(A, axis=1), np.linalg.norm(B, axis=1)
    return (A / np.where(na > 0, na, 1)[:, None]) @ (B / np.where(nb > 0, nb, 1)[:, None]).T


def mass(v):
    a = np.abs(np.nan_to_num(v)); return np.array([a[nodes == n].sum() for n in NETS]) / a.sum()


def load(mode, K):
    if mode == "spatial":
        root = {46: SR / "k46-networks/seeds", 51: SR / "k51-networks/maps/seeds"}[K]
        M = {seed_of(p): np.load(p)["Smean"].astype(float) for p in root.rglob(f"spatialAA_across_acrossCond_K{K}.npz")}
        F = M                                      # match on maps
    else:
        M, F = {}, {}
        for p in (SR / "k-temporal-maps/seeds").rglob(f"temporalAA_across_acrossCond_K{K}.npz"):
            d = np.load(p); s = seed_of(p); F[s] = d["C"].astype(float).T; M[s] = d["Smean"].astype(float)
        f = np.load(AA / f"msaa-rerun-seqnum-2026-10-27/msaa_flexible_outputs_npz/temporalAA_across_acrossCond_K{K}.npz", allow_pickle=True)
        F[0] = np.asarray(f["C"].tolist(), float).T
        M[0] = np.stack([np.asarray(r["S"], float) for r in f["results_subj"].tolist()]).mean(0)
    assert sorted(F) == SEEDS, (mode, K, sorted(F))
    match, R = {}, {}
    for s in SEEDS:
        R[s] = rowcorr(F[0], F[s]); i, j = linear_sum_assignment(-R[s]); match[s] = dict(zip(i.tolist(), j.tolist()))
    return M, F, match, R


# ---------------------------------------------------------------- 1. K=46 paper archetypes
say("=" * 100); say("1. K=46 alignment of the paper archetypes (legacy label = seqnum/seed-0 id)"); say("=" * 100)
M46, _, match46, R46 = load("spatial", 46)
sim = pd.read_csv(SR / "per-seed-networks/K46/K46_crossseed_similarity_all.csv")
for s in SEEDS:  # same matching as the per-seed figure
    sub = sim[sim.seed == s]; assert all(match46[s][a] == m for a, m in zip(sub.archetype0, sub.matched)), s
rows = []
for a, lab in PAPER.items():
    m0 = mass(M46[0][a])
    for s in SEEDS:
        b = match46[s][a]; ms = mass(M46[s][b]); r = R46[s][a].copy(); rb = r[b]; r[b] = -np.inf
        t = sim[(sim.seed == s) & (sim.archetype0 == a)].iloc[0]
        rows.append(dict(paper_label=lab, archetype_seed0=a, seed=s, matched_archetype=b, map_r=rb,
                         best_other_map_r=r.max(), recovered=rb >= REP_R, timecourse_r=t.tc_r,
                         yeo7_profile_r=np.corrcoef(m0, ms)[0, 1], dominant_net_seed0=NETS[m0.argmax()],
                         dominant_net=NETS[ms.argmax()], same_dominant_net=m0.argmax() == ms.argmax(),
                         **{f"yeo_{n}": v for n, v in zip(NSHORT, ms)}))
A = pd.DataFrame(rows); A.to_csv(OUT / "K46_paper_archetypes_alignment.csv", index=False)
sim.to_csv(OUT / "K46_all_archetypes_alignment.csv", index=False)
o = A[A.seed > 0]
say(o.groupby(["paper_label", "archetype_seed0"], sort=False).agg(
    map_r_median=("map_r", "median"), map_r_min=("map_r", "min"), recovered=("recovered", "sum"),
    tc_r_min=("timecourse_r", "min"), yeo_r_min=("yeo7_profile_r", "min"), same_net=("same_dominant_net", "sum"),
    net=("dominant_net_seed0", "first")).round(3).to_string())
say("(recovered / same_net = out of 9 seeds; recovered = matched map r >= 0.9)")
say("Failures (map r < 0.9):"); say(o[~o.recovered][["paper_label", "archetype_seed0", "seed", "map_r", "best_other_map_r", "dominant_net"]].round(3).to_string(index=False))
allK = sim[sim.seed > 0]
say(f"All 46 archetypes, seeds 1-9: {np.mean(allK.map_r >= REP_R):.0%} of matches have map r >= 0.9; median {allK.map_r.median():.3f}")

# ---------------------------------------------------------------- 2. held-out top-5 stability
say(""); say("=" * 100); say("2. Stability of each seed's own held-out top-5 selection"); say("=" * 100)


def jac(a, b): a, b = set(a), set(b); return len(a & b) / len(a | b)


def topm_runs(mode, K):
    cache = OUT / f".runs_{mode}_K{K}.csv"
    if cache.exists(): return pd.read_csv(cache)
    parts = []
    for ch in pd.read_csv(HE / f"heldout-everyk-{mode}-all-results/heldout_topm_runs_all_seeds.csv.gz",
                          usecols=["top_m", "selected_archetypes", "condition", "K", "outer_rep", "direction", "fit_seed"], chunksize=500000):
        parts.append(ch[(ch.K == K) & (ch.top_m == 5)])
    d = pd.concat(parts); d.to_csv(cache, index=False); return d


ALL_SUM, STAB = [], {}
for mode, K in [("spatial", 46), ("spatial", 51), ("temporal", 14)]:
    M, F, match, R = (M46, None, match46, R46) if K == 46 else load(mode, K)
    inv = {s: {b: a for a, b in match[s].items()} for s in SEEDS}
    rec = {s: {a: R[s][a, b] >= REP_R for a, b in match[s].items()} for s in SEEDS}
    fr = []
    for s in SEEDS:
        d = pd.read_csv(HE / f"heldout-everyk-{mode}-all-results/seed{s:02d}_heldout_archetype_selection_stability.csv")
        d = d[(d.K == K) & (d.top_m == 5)].copy(); d["seed"] = s; d["archetype0"] = d.archetype.map(inv[s])
        fr.append(d)
    fr = pd.concat(fr)
    if mode == "spatial":  # same mapping as hsummary's selection_by_seed.csv
        chk = pd.read_csv(HE / "heldout-everyk-spatial-all-results/K46_top5.csv") if K == 46 else None
        if chk is not None:
            p = fr.pivot_table(index=["archetype0", "condition"], columns="seed", values="selection_frequency", fill_value=0)
            for _, r in chk.iterrows():
                assert np.allclose([p.loc[(r.archetype0, r.condition), s] for s in SEEDS], [r[str(s)] for s in SEEDS]), r.archetype0
    runs = topm_runs(mode, K)
    runs["sel0"] = [tuple(sorted(inv[s][int(x)] for x in re.findall(r"\d+", t))) for s, t in zip(runs.fit_seed, runs.selected_archetypes)]
    rng = np.random.default_rng(0); null = np.mean([jac(rng.choice(K, 5, replace=False), rng.choice(K, 5, replace=False)) for _ in range(20000)])
    for c in CONDS:
        f = fr[fr.condition == c].sort_values(["seed", "selection_frequency", "mean_rank_position"], ascending=[True, False, True])
        top = {s: f[f.seed == s].head(5).archetype0.tolist() for s in SEEDS}
        J = np.ones((10, 10))
        for a, b in combinations(SEEDS, 2): J[a, b] = J[b, a] = jac(top[a], top[b])
        pj = J[np.triu_indices(10, 1)]
        rc = runs[runs.condition == c]; byk = {s: g.set_index(["outer_rep", "direction"]).sel0 for s, g in rc.groupby("fit_seed")}
        keys = byk[0].index
        same_split = [jac(byk[a][k], byk[b][k]) for k in keys for a, b in combinations(SEEDS, 2)]
        within = [jac(g.iloc[i], g.iloc[j]) for g in byk.values() for i, j in combinations(range(len(g)), 2)]
        cnt = pd.Series([a for s in SEEDS for a in top[s]]).value_counts()
        meanfreq = f.groupby("archetype0").selection_frequency.mean()
        unrec = sum(not rec[s][a] for s in SEEDS[1:] for a in top[s])
        row = dict(analysis=mode, K=K, condition=c, consensus_set_seed0=sorted(top[0]),
                   jaccard_seedsets_mean=pj.mean(), jaccard_seedsets_min=pj.min(), jaccard_seedsets_median=np.median(pj),
                   frac_pairs_identical=np.mean(pj == 1), jaccard_same_split_across_seeds=np.mean(same_split),
                   jaccard_across_splits_within_seed=np.mean(within), jaccard_random=null,
                   n_selected_ever=len(cnt), core_10of10=sorted(cnt[cnt == 10].index.tolist()),
                   selected_not_recovered=unrec)
        if mode == "spatial":
            comp = np.array([np.mean([mass(M[s][match[s][a]]) for a in top[s]], axis=0) for s in SEEDS])
            Y = np.corrcoef(comp); py = Y[np.triu_indices(10, 1)]
            row.update(yeo7_r_mean=py.mean(), yeo7_r_min=py.min(), **{f"yeo_{n}_mean": v for n, v in zip(NSHORT, comp.mean(0))},
                       **{f"yeo_{n}_sd": v for n, v in zip(NSHORT, comp.std(0, ddof=1))})
            pd.DataFrame(comp, columns=NSHORT, index=pd.Index(SEEDS, name="seed")).to_csv(OUT / f"{mode}_K{K}_{c}_top5_yeo7_by_seed.csv")
            pd.DataFrame(Y, index=SEEDS, columns=SEEDS).to_csv(OUT / f"{mode}_K{K}_{c}_top5_yeo7_r_pairwise.csv")
        else:
            Y = None
        ALL_SUM.append(row)
        sf = pd.DataFrame(dict(archetype0=cnt.index, n_seeds_in_top5=cnt.values))
        sf["mean_split_selection_freq"] = sf.archetype0.map(meanfreq)
        sf["seeds_recovered_r09"] = [sum(rec[s][a] for s in SEEDS) for a in sf.archetype0]
        if K == 46: sf["paper_label"] = sf.archetype0.map(PAPER).fillna("")
        sf.insert(0, "condition", c); sf.insert(0, "K", K); sf.insert(0, "analysis", mode)
        sf.sort_values(["n_seeds_in_top5", "mean_split_selection_freq"], ascending=False).to_csv(OUT / f"{mode}_K{K}_{c}_top5_selection_frequency.csv", index=False)
        pd.DataFrame(J, index=SEEDS, columns=SEEDS).to_csv(OUT / f"{mode}_K{K}_{c}_top5_jaccard_pairwise.csv")
        pd.DataFrame([dict(seed=s, top5_seed0_ids=sorted(top[s]), top5_own_ids=sorted(match[s][a] for a in top[s]),
                           all_recovered=all(rec[s][a] for a in top[s])) for s in SEEDS]).to_csv(OUT / f"{mode}_K{K}_{c}_top5_sets_by_seed.csv", index=False)
        STAB[(mode, K, c)] = dict(top=top, J=J, Y=Y, freq=f, cnt=cnt, same=np.array(same_split), within=np.array(within), null=null,
                                  comp=comp if mode == "spatial" else None, rec=rec)
S = pd.DataFrame(ALL_SUM); S.to_csv(OUT / "top5_stability_summary.csv", index=False)
cols = ["analysis", "K", "condition", "consensus_set_seed0", "core_10of10", "n_selected_ever", "jaccard_seedsets_mean", "jaccard_seedsets_min",
        "frac_pairs_identical", "jaccard_same_split_across_seeds", "jaccard_across_splits_within_seed", "jaccard_random", "yeo7_r_mean", "yeo7_r_min", "selected_not_recovered"]
pd.set_option("display.width", 250); pd.set_option("display.max_columns", 30)
say(S[cols].round(3).to_string(index=False))
say("jaccard_seedsets = the 10 seeds' consensus top-5 sets, 45 pairs; same_split = per-split sets on the identical participant split, seed pairs;")
say("across_splits_within_seed = per-split sets, different splits, same seed; random = two random 5-sets of K. selected_not_recovered = picks (seeds 1-9) whose seed-0 match has r < 0.9.")

# ---------------------------------------------------------------- 3. temporal merge check
say(""); say("=" * 100); say("3. Temporal K=14: A5+A9 subspace in seeds 4 and 9"); say("=" * 100)
MT, FT, matT, RT = load("temporal", 14)
PAIR = [5, 9]
C0 = FT[0][PAIR].T                                           # T x 2
Q, _ = np.linalg.qr(C0)
fitd = {}
for s in SEEDS:
    p = (SR / f"k-temporal-maps/seeds/full/seed{s:02d}/slim/temporalAA_across_acrossCond_K14.npz")
    fitd[s] = float(np.load(p)["varexpl"]) if p.exists() else float(np.load(AA / "msaa-rerun-seqnum-2026-10-27/msaa_flexible_outputs_npz/temporalAA_across_acrossCond_K14.npz", allow_pickle=True)["varexpl"])
op0 = C0 @ MT[0][PAIR]                                        # T x T: A5+A9 part of the reconstruction operator X -> X C S
trows = []
for s in SEEDS:
    Cs = FT[s].T                                              # T x K
    r2 = 1 - np.linalg.norm(Cs - Q @ (Q.T @ Cs), axis=0) ** 2 / np.linalg.norm(Cs, axis=0) ** 2   # uncentred: X C is linear in C
    inside = [int(k) for k in np.where(r2 >= 0.5)[0]]
    ops = Cs[:, inside] @ MT[s][inside]
    rop = np.corrcoef(op0.ravel(), ops.ravel())[0, 1]
    # principal angles between span{C0_5, C0_9} and the seed's archetypes in that region
    Qs, _ = np.linalg.qr(Cs[:, inside]); sv = np.clip(np.linalg.svd(Q.T @ Qs, compute_uv=False), 0, 1)
    expr = np.corrcoef(MT[0][PAIR].sum(0), MT[s][inside].sum(0))[0, 1]
    trows.append(dict(seed=s, match_A5=matT[s][5], r_A5=RT[s][5, matT[s][5]], match_A9=matT[s][9], r_A9=RT[s][9, matT[s][9]],
                      seed_archetypes_in_A5A9_span=inside, R2_in_span=[round(float(r2[k]), 3) for k in inside],
                      cos_principal_angles=[round(float(x), 3) for x in sv], pair_reconstruction_r=rop,
                      summed_expression_r=expr, varexpl=fitd[s]))
T = pd.DataFrame(trows); T.to_csv(OUT / "temporal_K14_A5A9_subspace.csv", index=False)
say(T.round(3).to_string(index=False))
say("R2_in_span: share of each seed archetype's generator column (C) inside span{seed-0 C_A5, C_A9}; archetypes with R2 >= 0.5 listed.")
say("pair_reconstruction_r: r between seed-0 C[:,{5,9}] S[{5,9},:] and the same T x T operator built from the listed seed archetypes")
say("(the part of the group-mean reconstruction X C S those archetypes carry). summed_expression_r: r of summed group-mean expression.")

# ---------------------------------------------------------------- figures
NC = plt.cm.viridis


def heat(ax, Mx, rl, cl, title, vmin=0, vmax=1, fmt="{:.2f}", cmap="viridis", bad=None):
    im = ax.imshow(Mx, vmin=vmin, vmax=vmax, cmap=cmap, aspect="auto")
    for i in range(Mx.shape[0]):
        for j in range(Mx.shape[1]):
            v = Mx[i, j]; ax.text(j, i, fmt.format(v), ha="center", va="center", fontsize=8.5,
                                  color="white" if v < vmin + 0.55 * (vmax - vmin) else "black",
                                  fontweight="bold" if bad is not None and bad[i, j] else "normal")
    ax.set_xticks(range(len(cl))); ax.set_xticklabels(cl, fontsize=9); ax.set_yticks(range(len(rl))); ax.set_yticklabels(rl, fontsize=10)
    ax.set_title(title, fontsize=11, pad=6); return im


# Fig A: K=46 paper archetypes
lab = [f"{l} (#{a})" for a, l in PAPER.items()]
fig = plt.figure(figsize=(17.5, 5.2)); gs = fig.add_gridspec(1, 4, width_ratios=[1, 1, 1, 0.9], wspace=0.28)
for k, (col, t) in enumerate([("map_r", "Spatial map r (matched)"), ("timecourse_r", "Time-course r (scan-wise, Fisher-z mean)"),
                               ("yeo7_profile_r", "Yeo-7 composition r")]):
    Mx = o.pivot(index="archetype_seed0", columns="seed", values=col).loc[list(PAPER)].to_numpy()
    bad = (o.pivot(index="archetype_seed0", columns="seed", values="map_r").loc[list(PAPER)].to_numpy() < REP_R)
    ax = fig.add_subplot(gs[k]); im = heat(ax, Mx, lab if k == 0 else [""] * 5, [str(s) for s in range(1, 10)], t, bad=bad)
    ax.set_xlabel("seed", fontsize=10)
ax = fig.add_subplot(gs[3]); x = np.arange(5)
for i, a in enumerate(PAPER):
    comp = A[A.archetype_seed0 == a][[f"yeo_{n}" for n in NSHORT]].to_numpy()  # 10 x 7
    left = 0
    for n in range(7):
        ax.barh(i, comp[0, n], left=left, color=NCOL[n], edgecolor="white", lw=0.5, height=0.62)
        left += comp[0, n]
    ax.text(1.01, i, f"{comp[:, 6].mean():.2f}±{comp[:, 6].std(ddof=1):.2f} DMN" if NETS[np.argmax(comp[0])] == "Default mode" else
            f"{comp[:, np.argmax(comp[0])].mean():.2f}±{comp[:, np.argmax(comp[0])].std(ddof=1):.2f} {NSHORT[np.argmax(comp[0])]}",
            va="center", fontsize=8.5, transform=ax.get_yaxis_transform())
ax.set_yticks(x); ax.set_yticklabels([""] * 5); ax.invert_yaxis(); ax.set_xlim(0, 1); ax.set_xlabel("share of |weight|", fontsize=10)
ax.set_title("Yeo-7 composition, seed 0\n(text: dominant net mean±SD over 10 seeds)", fontsize=10.5)
for n in range(7): ax.bar(0, 0, color=NCOL[n], label=NSHORT[n])
ax.legend(ncol=4, fontsize=7.5, loc="upper center", bbox_to_anchor=(0.5, -0.17), frameon=False)
fig.colorbar(im, ax=fig.axes[:3], fraction=0.015, pad=0.01)
fig.suptitle("K=46 paper archetypes across nine further initialisations (Hungarian-matched to seed 0; bold = map r < 0.9, i.e. not recovered)",
             fontsize=12.5, fontweight="bold", y=1.02)
for e in ("png", "pdf"): fig.savefig(OUT / f"K46_paper_archetypes_alignment.{e}", dpi=220 if e == "png" else None, bbox_inches="tight")
plt.close(fig)

# Fig B: top-5 overlap summary per (mode, K)
for mode, K in [("spatial", 46), ("spatial", 51), ("temporal", 14)]:
    ncol = 3 if mode == "spatial" else 2
    fig, axs = plt.subplots(3, ncol, figsize=(6.2 * ncol, 12), gridspec_kw=dict(width_ratios=[1.5, 1, 1.15][:ncol]), squeeze=False)
    for ri, c in enumerate(CONDS):
        st = STAB[(mode, K, c)]; cnt = st["cnt"]; ids = cnt.index.tolist()
        ids = sorted(ids, key=lambda a: (-cnt[a], -st["freq"][st["freq"].archetype0 == a].selection_frequency.mean()))
        P = st["freq"].pivot_table(index="archetype0", columns="seed", values="selection_frequency", fill_value=0).reindex(ids).fillna(0).to_numpy()
        ax = axs[ri, 0]; ax.imshow(P, vmin=0, vmax=1, cmap="Blues", aspect="auto")
        for i, a in enumerate(ids):
            for s in SEEDS:
                if a in st["top"][s]: ax.add_patch(plt.Rectangle((s - 0.5, i - 0.5), 1, 1, fill=False, ec="#E4572E", lw=2))
                if not st["rec"][s][a]: ax.text(s, i, "×", ha="center", va="center", color="0.4", fontsize=9)
        lab = [f"#{a}" + (f" {PAPER[a]}" if K == 46 and a in PAPER else "") + f"  ({cnt[a]}/10)" for a in ids]
        ax.set_yticks(range(len(ids))); ax.set_yticklabels(lab, fontsize=8.5); ax.set_xticks(SEEDS); ax.set_xlabel("seed")
        ax.set_title(f"{c}: split-level selection freq (blue); orange box = in seed's top-5", fontsize=10)
        ax = axs[ri, 1]; pj = st["J"][np.triu_indices(10, 1)]
        data = [pj, st["same"], st["within"]]
        bp = ax.boxplot(data, widths=0.55, patch_artist=True, showfliers=False)
        for b, col in zip(bp["boxes"], ["#E4572E", "#4C9BE8", "#9AA5B1"]): b.set_facecolor(col); b.set_alpha(0.6)
        for i, d in enumerate(data): ax.scatter(np.random.default_rng(i).normal(i + 1, 0.06, len(d)), d, s=4, color="k", alpha=0.25 if len(d) > 100 else 0.6, zorder=3)
        ax.axhline(st["null"], ls="--", color="0.3", lw=1); ax.text(3.45, st["null"], "random", fontsize=8, va="bottom", ha="right")
        ax.set_xticks([1, 2, 3]); ax.set_xticklabels(["seed top-5 sets\n(45 pairs)", "same split,\nacross seeds", "across splits,\nwithin seed"], fontsize=8.5)
        ax.set_ylim(-0.03, 1.05); ax.set_ylabel("Jaccard"); ax.set_title(f"{c}: overlap of top-5 sets (mean {pj.mean():.2f})", fontsize=10)
        if mode == "spatial":
            ax = axs[ri, 2]; comp = st["comp"]; left = np.zeros(10)
            for n in range(7): ax.bar(SEEDS, comp[:, n], bottom=left, color=NCOL[n], edgecolor="white", lw=0.4, label=NSHORT[n]); left += comp[:, n]
            py = st["Y"][np.triu_indices(10, 1)]
            ax.set_xticks(SEEDS); ax.set_xlabel("seed"); ax.set_ylim(0, 1); ax.set_ylabel("share of |weight| (mean of top-5)")
            ax.set_title(f"{c}: Yeo-7 of each seed's top-5\npairwise r mean {py.mean():.3f}, min {py.min():.3f}", fontsize=10)
            if ri == 0: ax.legend(ncol=4, fontsize=7.5, loc="upper center", bbox_to_anchor=(0.5, 1.42), frameon=False)
    fig.suptitle(f"{mode.capitalize()} K={K}: does each seed's own participant-held-out ranking pick the same top-5? (seed-0 IDs; × = seed's match has map r < 0.9)"
                 if mode == "spatial" else
                 f"Temporal K={K}: does each seed's own participant-held-out ranking pick the same top-5? (seed-0 IDs, matched on C; × = C r < 0.9)",
                 fontsize=12.5, fontweight="bold")
    fig.tight_layout(rect=(0, 0, 1, 0.965))
    for e in ("png", "pdf"): fig.savefig(OUT / f"top5_overlap_{mode}_K{K}.{e}", dpi=200 if e == "png" else None, bbox_inches="tight")
    plt.close(fig)
say(""); say("figures:", ", ".join(sorted(p.name for p in OUT.glob("*.png"))))
