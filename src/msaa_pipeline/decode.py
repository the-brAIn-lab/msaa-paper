from __future__ import annotations
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.spatial.distance import cdist

CONDITIONS = ["intact", "word", "rest"]


def load_fit(path: Path):
    d = np.load(path, allow_pickle=True)
    results = d["results_subj"].tolist()
    if isinstance(results, np.ndarray):
        results = results.tolist()
    labels = d["condition_labels_str"].tolist() if "condition_labels_str" in d else None
    return {"results_subj": results, "condition_labels_str": labels, "K": int(d["K"])}


def fit_path(output_root: Path, analysis_type, fit_scope, K, condition=None):
    base = output_root / "msaa_flexible_outputs_npz"
    if fit_scope == "across":
        return base / f"{analysis_type}AA_across_acrossCond_K{K}.npz"
    return base / f"{analysis_type}AA_within_{condition}_K{K}.npz"


def get_subjects(output_root: Path, analysis_type, fit_scope, K, condition):
    p = fit_path(output_root, analysis_type, fit_scope, K, condition)
    if not p.exists():
        raise FileNotFoundError(p)
    fit = load_fit(p)
    if fit_scope == "across":
        labels = np.asarray(fit["condition_labels_str"], dtype=object)
        idx = np.where(labels == condition)[0]
        subjects = [fit["results_subj"][i] for i in idx]
    else:
        subjects = fit["results_subj"]
    return subjects, p


def reconstruct(sub, analysis_type, archetype_indices=None):
    sXC = np.asarray(sub["sXC"], dtype=float)
    S = np.asarray(sub["S"], dtype=float)
    idx = np.arange(S.shape[0]) if archetype_indices is None else np.asarray(list(archetype_indices), dtype=int)
    if analysis_type == "spatial":
        Xhat = sXC[:, idx] @ S[idx, :]
    elif analysis_type == "temporal":
        Xhat = (sXC[:, idx] @ S[idx, :]).T
    else:
        raise ValueError(analysis_type)
    return (Xhat - Xhat.mean(axis=0, keepdims=True)) / (Xhat.std(axis=0, keepdims=True) + 1e-8)


def _assignments(n, nfolds, rng):
    out = np.zeros(n, dtype=int)
    groupsize = int(np.ceil(n / nfolds))
    for i in range(1, nfolds):
        out[i * groupsize:min((i + 1) * groupsize, n)] = i
    rng.shuffle(out)
    return out


def _decode_corr(corrs):
    T = corrs.shape[0]
    accuracy = rank = error = 0.0
    for t in range(T):
        decoded = int(np.argmax(corrs[t]))
        error += abs(decoded - t) / T
        accuracy += (decoded == t)
        rank += np.mean((corrs[t] <= corrs[t, t]).astype(int))
    return dict(accuracy=accuracy / T, rank=rank / T, error=error / T)


def decode_subjects(subjects, analysis_type, archetype_indices=None, nfolds=2, nreps=20, seed=42):
    stack = np.stack([reconstruct(s, analysis_type, archetype_indices) for s in subjects])
    rng = np.random.default_rng(seed)
    rows = []
    for rep in range(int(nreps)):
        folds = _assignments(len(subjects), int(nfolds), rng)
        for fold in range(int(nfolds)):
            a = stack[folds == fold].mean(axis=0)
            b = stack[folds != fold].mean(axis=0)
            corrs = 1.0 - cdist(a, b, metric="correlation")
            corrs = np.nan_to_num(corrs, nan=0.0, posinf=0.0, neginf=0.0)
            row = _decode_corr(corrs)
            row.update(rep=rep, fold=fold)
            rows.append(row)
    return pd.DataFrame(rows)


def summarize(df):
    return {
        "mean": float(df.accuracy.mean()),
        "err": float(df.accuracy.std() / np.sqrt(max(len(df), 1))),
        "std": float(df.accuracy.std()),
        "count": int(len(df)),
        "rank_mean": float(df["rank"].mean()),
        "error_mean": float(df.error.mean()),
    }


def _write(df, path):
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, index=False)


def _load_csv(path: Path):
    return pd.read_csv(path) if path.exists() else pd.DataFrame()


def _keyset(df: pd.DataFrame, cols):
    if df.empty or not all(c in df.columns for c in cols):
        return set()
    return set(map(tuple, df[cols].astype(str).itertuples(index=False, name=None)))


def _append_unique(path: Path, rows, key_cols):
    if not rows:
        return _load_csv(path)
    new = pd.DataFrame(rows)
    old = _load_csv(path)
    if old.empty:
        combined = new
    else:
        keys = _keyset(old, key_cols)
        mask = [tuple(map(str, row)) not in keys for row in new[key_cols].itertuples(index=False, name=None)]
        combined = pd.concat([old, new.loc[mask]], ignore_index=True)
    sort_cols = [c for c in key_cols if c in combined.columns]
    if sort_cols:
        combined = combined.sort_values(sort_cols).reset_index(drop=True)
    _write(combined, path)
    return combined


def run_decode_job(output_root: Path, analysis_type: str, fit_scope: str, K_values,
                   *, top_m=(1,2,3,5,10), nfolds=2, nreps=20, seed=42):
    """Run/resume full, per-archetype, and top-m decoding.

    CSVs are checkpointed after each condition so a long full sweep can be
    resumed safely. Existing keyed rows are not recomputed.
    """
    outdir = output_root / f"msaa_condrank_decoding_outputs_{analysis_type}_{fit_scope}"
    outdir.mkdir(parents=True, exist_ok=True)
    full_path = outdir / "full_summary.csv"
    arch_path = outdir / "per_archetype_mean_accuracy.csv"
    top_path = outdir / "topm_summary.csv"
    rankings_path = outdir / "rankings.npy"

    full = _load_csv(full_path)
    arch = _load_csv(arch_path)
    top = _load_csv(top_path)
    full_keys = _keyset(full, ["analysis_type","fit_scope","condition","K"])
    arch_keys = _keyset(arch, ["analysis_type","fit_scope","condition","K","archetype"])
    top_keys = _keyset(top, ["analysis_type","fit_scope","condition","K","top_m"])

    for K in K_values:
        for condition in CONDITIONS:
            subjects, source = get_subjects(output_root, analysis_type, fit_scope, int(K), condition)

            fkey = tuple(map(str, [analysis_type, fit_scope, condition, int(K)]))
            if fkey not in full_keys:
                d = decode_subjects(subjects, analysis_type, None, nfolds, nreps, seed)
                rows = [dict(analysis_type=analysis_type, fit_scope=fit_scope, condition=condition,
                             K=int(K), fit_source_file=str(source), fit_source_basename=source.name,
                             **summarize(d))]
                full = _append_unique(full_path, rows, ["analysis_type","fit_scope","condition","K"])
                full_keys.add(fkey)

            n_arch = np.asarray(subjects[0]["S"]).shape[0]
            new_arch_rows = []
            for a in range(n_arch):
                akey = tuple(map(str, [analysis_type, fit_scope, condition, int(K), int(a)]))
                if akey in arch_keys:
                    continue
                dd = decode_subjects(subjects, analysis_type, [a], nfolds, nreps, seed)
                ss = summarize(dd)
                new_arch_rows.append(dict(
                    analysis_type=analysis_type, fit_scope=fit_scope, condition=condition,
                    K=int(K), archetype=int(a), fit_source_file=str(source), fit_source_basename=source.name,
                    mean_accuracy=ss["mean"], sem_accuracy=ss["err"], std_accuracy=ss["std"],
                    count=ss["count"], rank_mean=ss["rank_mean"], error_mean=ss["error_mean"]
                ))
                arch_keys.add(akey)
            if new_arch_rows:
                arch = _append_unique(arch_path, new_arch_rows,
                                      ["analysis_type","fit_scope","condition","K","archetype"])

            cur = arch[(arch["condition"].astype(str) == condition) & (arch["K"].astype(int) == int(K))]
            if len(cur) < n_arch:
                raise RuntimeError(f"Incomplete per-archetype table for {analysis_type}/{fit_scope}/{condition}/K={K}")
            ranked = (
                cur.groupby("archetype")["mean_accuracy"].mean()
                .sort_values(ascending=False).index.astype(int).tolist()
            )

            new_top_rows = []
            for m in top_m:
                if int(m) > len(ranked):
                    continue
                tkey = tuple(map(str, [analysis_type, fit_scope, condition, int(K), int(m)]))
                if tkey in top_keys:
                    continue
                selected = ranked[:int(m)]
                dd = decode_subjects(subjects, analysis_type, selected, nfolds, nreps, seed)
                new_top_rows.append(dict(
                    analysis_type=analysis_type, fit_scope=fit_scope, condition=condition,
                    K=int(K), top_m=int(m), selected_archetypes=str(selected),
                    fit_source_file=str(source), fit_source_basename=source.name,
                    **summarize(dd)
                ))
                top_keys.add(tkey)
            if new_top_rows:
                top = _append_unique(top_path, new_top_rows,
                                     ["analysis_type","fit_scope","condition","K","top_m"])

            rankings = np.load(rankings_path, allow_pickle=True).item() if rankings_path.exists() else {}
            rankings[(str(condition), int(K))] = ranked
            np.save(rankings_path, rankings, allow_pickle=True)

    return _load_csv(full_path), _load_csv(arch_path), _load_csv(top_path)

