"""Selected temporal archetypes on one figure, laid out as glassbrains_all_seeds.py: rows = seeds, columns = archetypes.
Motif = group-mean expression over the scan (row of Smean, K x T). y-limits are shared. Seed 0 comes from the full fit, with +- SEM over the 108
scans (slim files keep no per-scan S); seeds >= 1 from slim npz, seed-0 motif in grey behind. Matching:
Hungarian to seed 0 on the r of C columns (zero lag), so columns keep seed-0 IDs. Where the one-to-one match has r < MERGE_R
the seed has no separate copy of that archetype, so the panel is left blank instead of showing the unrelated leftover, and is
labelled with the seed-0 archetype that the seed's best-correlated archetype was matched to.
Usage: python temporal_motifs_by_seed.py K IDS(comma) [--fit FULLFIT_seed0.npz] [--out OUTDIR]
Defaults: the seed-0 fit in msaa_flexible_outputs_npz/, seeds 1-9 from seeds/ (paths.py),
OUTDIR outputs/multiseed_stability/temporal_motifs_by_seed.
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
ap.add_argument("--fit", default=None); ap.add_argument("--out", default=str(paths.OUT / "temporal_motifs_by_seed"))
a = ap.parse_args(); K = a.K; SEL = [int(x) for x in a.ids.split(",")]
FIT = Path(a.fit or paths.full_fit0("temporal", K)); OUT = Path(a.out); OUT.mkdir(parents=True, exist_ok=True)
MERGE_R = 0.5
f = np.load(FIT, allow_pickle=True); assert int(f["K"]) == K
S0 = np.stack([np.asarray(r["S"], float) for r in f["results_subj"].tolist()])          # scans x K x T
C = {0: np.asarray(f["C"].tolist(), float)}; S = {0: S0.mean(0)}; SE0 = S0.std(0, ddof=1) / np.sqrt(len(S0))
for s, p in paths.seed_fits("temporal", K, seeds=range(1, paths.N_SEEDS)).items():
    d = np.load(p); C[s] = d["C"].astype(float); S[s] = d["Smean"].astype(float)
seeds = sorted(C); print("seeds:", seeds); assert all(S[s].shape == (K, C[0].shape[0]) and C[s].shape == C[0].shape for s in seeds)

def corr(A, B):
    A = A - A.mean(0); B = B - B.mean(0); na, nb = np.linalg.norm(A, axis=0), np.linalg.norm(B, axis=0)
    return (A / np.where(na > 0, na, 1)).T @ (B / np.where(nb > 0, nb, 1))

rows = []; fig, axs = plt.subplots(len(seeds), len(SEL), figsize=(3.6 * len(SEL) + 1, 1.45 * len(seeds) + 1.1), sharex=True, sharey=True, squeeze=False)
t = np.arange(S[0].shape[1]); match = {}
for ri, s in enumerate(seeds):
    best = corr(C[0], C[s]); i, j = linear_sum_assignment(-best); m = dict(zip(i, j)); match[s] = m; inv = dict(zip(j, i))
    for ci, k in enumerate(SEL):
        b = m[k]; r = best[k, b]; ax = axs[ri, ci]; mrg = r < MERGE_R
        near = int(best[k].argmax()) if mrg else b; g = S[s][near]
        if s == 0: ax.fill_between(t, g - SE0[k], g + SE0[k], color="#1C75BC", alpha=0.25, lw=0)
        elif not mrg: ax.plot(t, S[0][k], color="0.75", lw=1.0)
        if not mrg: ax.plot(t, g, color="#1C75BC", lw=1.1)
        if mrg: ax.text(0.5, 0.5, f"no separate A{k} (merged with A{inv[near]})", transform=ax.transAxes, ha="center", va="center", fontsize=10.5, color="red")
        else: ax.set_title(f"#{b} · r={r:.2f}", fontsize=10, color="red" if r < 0.8 else "0.25", pad=3)
        ax.set_xlim(0, t[-1]); ax.set_ylim(0, 1.02); ax.set_yticks([0, 0.5, 1]); ax.tick_params(labelsize=9); ax.spines[["top", "right"]].set_visible(False)
        if ci == 0: ax.set_ylabel(f"seed {s}", rotation=0, labelpad=30, fontsize=12, fontweight="bold", va="center")
        if ri == len(seeds) - 1: ax.set_xlabel("TR", fontsize=10)
        rows.append(dict(seed=s, archetype0=k, matched=int(b), C_r=r, closest=near, closest_C_r=best[k, near], motif_r=np.corrcoef(S[0][k], g)[0, 1],
                         peak_TR=int(C[s][:, near].argmax()), merged_with=int(inv[near]) if mrg else -1))
for ci, k in enumerate(SEL):
    p = axs[0, ci].get_position(); axs[0, ci].annotate(f"A{k}", (0.5, 1.28), xycoords="axes fraction", ha="center", fontsize=15, fontweight="bold")
fig.suptitle(f"K={K} temporal archetypes across {len(seeds)} seeds: mean expression (seed-0 IDs)\n"
             "Panel text: seed's own # · C r with seed 0. Grey = seed-0 motif; band on seed 0 = ±SEM over 108 scans; blank = no separate archetype in that seed. "
             "The peak near 1 at the TRs that build each archetype is by construction.", fontsize=13, fontweight="bold")
fig.tight_layout(rect=(0, 0, 1, 0.985), h_pad=1.0)
for e in ("png", "pdf"): fig.savefig(OUT / f"K{K}_temporal_allseeds_motifs.{e}", dpi=250 if e == "png" else None, bbox_inches="tight")
plt.close(fig)

P = pd.DataFrame(rows); P.to_csv(OUT / f"K{K}_temporal_allseeds_motifs.csv", index=False)
print(P.pivot(index="archetype0", columns="seed", values="closest_C_r").loc[SEL].round(2).to_string())
print(P.pivot(index="archetype0", columns="seed", values="motif_r").loc[SEL].round(2).to_string())
