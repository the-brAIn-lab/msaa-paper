#!/usr/bin/env python3
from __future__ import annotations

import argparse
import ast
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import linear_sum_assignment

FIT_NAME = "spatialAA_across_acrossCond_K46.npz"
OLD_SELECTED = [4, 5, 6, 12, 8]

def find_new_fit(root: Path) -> Path:
    candidates = [
        root / "outputs" / "mini_paper" / "msaa_flexible_outputs_npz" / FIT_NAME,
        root / "msaa_flexible_outputs_npz" / FIT_NAME,
    ]
    for p in candidates:
        if p.exists():
            return p
    hits = list(root.glob(f"**/{FIT_NAME}"))
    if hits:
        return hits[0]
    raise FileNotFoundError(f"Could not find {FIT_NAME} under {root}")

def find_old_fit(root: Path) -> Path:
    candidates = [
        root / "msaa_flexible_outputs_npz" / FIT_NAME,
        root / "outputs" / "mini_paper" / "msaa_flexible_outputs_npz" / FIT_NAME,
    ]
    for p in candidates:
        if p.exists():
            return p
    hits = list(root.glob(f"**/{FIT_NAME}"))
    if hits:
        return hits[0]
    raise FileNotFoundError(f"Could not find {FIT_NAME} under {root}")

def load_subject_results(path: Path):
    d = np.load(path, allow_pickle=True)
    results = d["results_subj"].tolist()
    if isinstance(results, np.ndarray):
        results = results.tolist()
    return list(results)

def mean_spatial_maps(path: Path) -> np.ndarray:
    results = load_subject_results(path)
    maps = []
    for sub in results:
        maps.append(np.asarray(sub["S"], dtype=float))
    stacked = np.stack(maps, axis=0)
    return np.nanmean(stacked, axis=0)

def row_corr(A: np.ndarray, B: np.ndarray) -> np.ndarray:
    A = A - np.nanmean(A, axis=1, keepdims=True)
    B = B - np.nanmean(B, axis=1, keepdims=True)
    A /= np.sqrt(np.nansum(A*A, axis=1, keepdims=True))
    B /= np.sqrt(np.nansum(B*B, axis=1, keepdims=True))
    return np.clip(A @ B.T, -1.0, 1.0)

def parse_set(x):
    try:
        v = ast.literal_eval(str(x))
        return set(int(z) for z in v)
    except Exception:
        return set()

def maybe_compare_heldout(new_root, old_root, new_to_old, outdir):
    new_candidates = [
        new_root / "outputs" / "mini_paper" / "msaa_heldout_topm_outputs_spatial_across" / "heldout_topm_crossfit_runs.csv",
        new_root / "msaa_heldout_topm_outputs_spatial_across" / "heldout_topm_crossfit_runs.csv",
    ]
    old_candidates = [
        old_root / "msaa_heldout_topm_outputs_spatial_across" / "heldout_topm_crossfit_runs.csv",
    ]
    new_path = next((p for p in new_candidates if p.exists()), None)
    old_path = next((p for p in old_candidates if p.exists()), None)
    if new_path is None or old_path is None:
        return {"status": "skipped"}

    new = pd.read_csv(new_path)
    old = pd.read_csv(old_path)
    new = new[new["K"].astype(int) == 46].copy()
    old = old[old["K"].astype(int) == 46].copy()
    keys = ["condition", "K", "outer_rep", "direction", "top_m"]
    m = new.merge(old, on=keys, suffixes=("_new", "_old"))

    m["set_new_raw"] = m["selected_archetypes_new"].apply(parse_set)
    m["set_new_aligned"] = m["set_new_raw"].apply(
        lambda s: {new_to_old.get(int(v), int(v)) for v in s}
    )
    m["set_old"] = m["selected_archetypes_old"].apply(parse_set)

    m["same_selection_raw"] = [a == b for a,b in zip(m["set_new_raw"], m["set_old"])]
    m["same_selection_aligned"] = [a == b for a,b in zip(m["set_new_aligned"], m["set_old"])]

    def jac(a,b):
        u = a | b
        return len(a & b) / len(u) if u else 1.0

    m["jaccard_raw"] = [jac(a,b) for a,b in zip(m["set_new_raw"], m["set_old"])]
    m["jaccard_aligned"] = [jac(a,b) for a,b in zip(m["set_new_aligned"], m["set_old"])]
    m["accuracy_diff"] = np.abs(m["mean_new"] - m["mean_old"])

    summary = (
        m.groupby("top_m", as_index=False)
        .agg(
            exact_agreement_raw=("same_selection_raw","mean"),
            exact_agreement_aligned=("same_selection_aligned","mean"),
            mean_jaccard_raw=("jaccard_raw","mean"),
            mean_jaccard_aligned=("jaccard_aligned","mean"),
            mean_abs_accuracy_diff=("accuracy_diff","mean"),
            max_abs_accuracy_diff=("accuracy_diff","max"),
            n_rows=("top_m","size"),
        )
    )
    summary.to_csv(outdir / "heldout_k46_overlap_after_alignment.csv", index=False)
    return {"status":"ok","by_top_m":summary.to_dict(orient="records")}

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--old-repo", required=True)
    ap.add_argument("--new-root", default=".")
    ap.add_argument("--new-fit")
    ap.add_argument("--old-fit")
    ap.add_argument("--outdir", default="outputs/mini_paper/k46_alignment")
    args = ap.parse_args()

    new_root = Path(args.new_root).resolve()
    old_root = Path(args.old_repo).resolve()
    outdir = Path(args.outdir)
    if not outdir.is_absolute():
        outdir = new_root / outdir
    outdir.mkdir(parents=True, exist_ok=True)

    new_fit = Path(args.new_fit).resolve() if args.new_fit else find_new_fit(new_root)
    old_fit = Path(args.old_fit).resolve() if args.old_fit else find_old_fit(old_root)

    old_maps = mean_spatial_maps(old_fit)
    new_maps = mean_spatial_maps(new_fit)
    if old_maps.shape != new_maps.shape:
        raise ValueError(f"Shape mismatch: {old_maps.shape} vs {new_maps.shape}")

    corr = row_corr(old_maps, new_maps)
    old_idx, new_idx = linear_sum_assignment(-np.nan_to_num(corr, nan=-1.0))

    alignment = pd.DataFrame({
        "old_archetype": old_idx.astype(int),
        "new_archetype": new_idx.astype(int),
        "correlation": corr[old_idx, new_idx],
    }).sort_values("old_archetype").reset_index(drop=True)
    alignment["abs_correlation"] = alignment["correlation"].abs()
    alignment.to_csv(outdir / "k46_archetype_alignment.csv", index=False)
    np.savetxt(outdir / "k46_old_vs_new_spatial_correlation_matrix.csv", corr, delimiter=",")

    old_to_new = dict(zip(alignment["old_archetype"].astype(int), alignment["new_archetype"].astype(int)))
    new_to_old = {v:k for k,v in old_to_new.items()}

    selected = alignment[alignment["old_archetype"].isin(OLD_SELECTED)].copy()
    selected["order"] = selected["old_archetype"].map({a:i for i,a in enumerate(OLD_SELECTED)})
    selected = selected.sort_values("order")
    selected.to_csv(outdir / "k46_manuscript_selected_mapping.csv", index=False)

    heldout = maybe_compare_heldout(new_root, old_root, new_to_old, outdir)

    summary = {
        "old_fit": str(old_fit),
        "new_fit": str(new_fit),
        "matched_correlation_mean": float(alignment["correlation"].mean()),
        "matched_correlation_median": float(alignment["correlation"].median()),
        "matched_correlation_min": float(alignment["correlation"].min()),
        "n_matched_corr_ge_0_90": int((alignment["correlation"] >= 0.90).sum()),
        "n_matched_corr_ge_0_95": int((alignment["correlation"] >= 0.95).sum()),
        "old_selected_mapping": [
            {
                "old_archetype": int(r.old_archetype),
                "new_archetype": int(r.new_archetype),
                "correlation": float(r.correlation),
            } for r in selected.itertuples()
        ],
        "heldout_after_alignment": heldout,
    }
    with open(outdir / "k46_alignment_summary.json","w") as f:
        json.dump(summary, f, indent=2)

    print("\nMatched-map correlation summary")
    print(alignment["correlation"].describe().to_string())
    print("\nManuscript-selected archetypes")
    print(selected[["old_archetype","new_archetype","correlation"]].to_string(index=False))
    if heldout.get("status") == "ok":
        print("\nHeld-out K=46 selection overlap after alignment")
        print(pd.DataFrame(heldout["by_top_m"]).to_string(index=False))
    print("\nSaved diagnostics to:", outdir)

if __name__ == "__main__":
    main()
