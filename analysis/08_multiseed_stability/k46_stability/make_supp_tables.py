"""Permanent M0-M4 crosswalk + compact K=46 supplementary table, from the topm_stability_neighbourhood.py outputs.

Usage: python make_supp_tables.py [--in TOPM_NEIGHBOURHOOD_OUTDIR] [--nodes NODES.csv] [--out OUTDIR]
Defaults (paths.py): outputs/multiseed_stability/topm_stability_neighbourhood (run it with K list including 46),
OUTDIR outputs/multiseed_stability/k46_supp_tables.
"""
import argparse, sys, numpy as np, pandas as pd
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))   # analysis/08_multiseed_stability, for paths.py
import paths  # noqa: E402
ap = argparse.ArgumentParser(); ap.add_argument("--in", dest="inp", default=str(paths.OUT / "topm_stability_neighbourhood"))
ap.add_argument("--nodes", default=str(paths.NODES)); ap.add_argument("--out", default=str(paths.OUT / "k46_supp_tables"))
a = ap.parse_args(); IN = Path(a.inp); OUT = Path(a.out); OUT.mkdir(parents=True, exist_ok=True)
MT = pd.read_csv(IN / "M0M4_selection_by_K.csv"); MR = pd.read_csv(IN / "M0M4_recovery_by_K.csv")
nodes = pd.read_csv(a.nodes)
nodes = nodes.drop_duplicates("node").set_index("node").sort_index()["network"].to_numpy()
S = np.load(paths.slim_fit("spatial", 46, 0))["Smean"].astype(float)
NETS = list(dict.fromkeys(["Visual", "Somatomotor", "Dorsal attention", "Ventral attention", "Limbic", "Frontoparietal", "Default mode"]))
r46 = MR[MR.K == 46].set_index("motif")
X = []
for m, row in r46.iterrows():
    a = np.abs(np.nan_to_num(S[int(row.k46_seed0_id)])); share = {n: a[nodes == n].sum() / a.sum() for n in NETS}
    top = sorted(share, key=share.get, reverse=True)
    X.append(dict(motif=m, k46_seed0_component_id_0based=int(row.k46_seed0_id), legacy_label=row.legacy,
                  in_submitted_paper_set=row.legacy in ("A4", "A5", "A6", "A12", "A8"),
                  dominant_yeo7=top[0], dominant_share=round(share[top[0]], 3), second_yeo7=top[1], second_share=round(share[top[1]], 3)))
X = pd.DataFrame(X); X.to_csv(OUT / "M0M4_crosswalk.csv", index=False); print(X.to_string(index=False))
d = MT[MT.K == 46]
rows = []
for (m, c), g in d.groupby(["motif", "condition"], sort=False):
    g5, g10 = g[g.top_m == 5].iloc[0], g[g.top_m == 10].iloc[0]; nrec = int(g5.n_recovered)
    rows.append(dict(motif=m, condition=c,
        recovered_seeds1to9=f"{nrec-1}/9", recovery_rate=round((nrec - 1) / 9, 2),
        min_map_r=r46.loc[m, "r_min"].round(3),
        top5_selected_and_recovered=f"{int(g5.n_selected_and_recovered)}/10",
        top10_selected_and_recovered=f"{int(g10.n_selected_and_recovered)}/10",
        selected_not_recovered_top5=int(g5.n_selected_forced_r_lt_090), selected_not_recovered_top10=int(g10.n_selected_forced_r_lt_090),
        P_top5_given_recovered=round(g5.frac_selected_given_recovered, 2), P_top10_given_recovered=round(g10.frac_selected_given_recovered, 2),
        split_freq_top5_given_recovered=round(g5.mean_split_selection_freq_given_recovered, 2),
        split_freq_top10_given_recovered=round(g10.mean_split_selection_freq_given_recovered, 2)))
T = pd.DataFrame(rows); T["condition"] = pd.Categorical(T.condition, ["intact", "word", "rest"])
T = T.sort_values(["condition", "motif"]); T.to_csv(OUT / "supp_table_K46_M0M4.csv", index=False)
pd.set_option("display.width", 250); print(T.to_string(index=False))
