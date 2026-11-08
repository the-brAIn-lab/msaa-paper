from __future__ import annotations
import os
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


def _row_corr(a, b):
    """Pearson correlation between rows of a and rows of b.

    Same values as 1 - cdist(a, b, "correlation") but via a BLAS matrix
    product, which is several times faster.
    """
    a = a - a.mean(axis=1, keepdims=True)
    b = b - b.mean(axis=1, keepdims=True)
    a = a / np.linalg.norm(a, axis=1, keepdims=True)
    b = b / np.linalg.norm(b, axis=1, keepdims=True)
    return a @ b.T


def _assignments(n, nfolds, rng):
    out = np.zeros(n, dtype=int)
    groupsize = int(np.ceil(n / nfolds))
    for i in range(1, nfolds):
        out[i * groupsize:min((i + 1) * groupsize, n)] = i
    rng.shuffle(out)
    return out


# Correlations within TIE_TOL of each other count as tied when scoring. 0.0 is
# the original exact comparison. Temporal S is mostly exact zeros, which makes
# many correlations mathematically equal; with exact comparison, whether they
# count as tied then depends on last-bit rounding (BLAS, summation order).
TIE_TOL = float(os.environ.get("MSAA_TIE_TOL", "0"))


def _decode_corr(corrs, tie_tol=None):
    """Accuracy, mean rank and mean error of time-point identification.

    Vectorised form of the original per-row loop. The per-row terms are the
    same, and np.cumsum adds them in the same left-to-right order as the loop's
    ``+=``, so with tie_tol=0 the result is bit-identical to the loop.
    With tie_tol > 0, the decoded time point is the first one within tie_tol
    of the row maximum, and the rank counts entries <= diagonal + tie_tol.
    """
    tol = TIE_TOL if tie_tol is None else float(tie_tol)
    T = corrs.shape[0]
    t = np.arange(T)
    if tol > 0:
        decoded = np.argmax(corrs >= corrs.max(axis=1, keepdims=True) - tol, axis=1)
    else:
        decoded = np.argmax(corrs, axis=1)
    accuracy = np.cumsum((decoded == t).astype(float))[-1]
    rank = np.cumsum((corrs <= corrs[t, t][:, None] + tol).astype(int).mean(axis=1))[-1]
    error = np.cumsum(np.abs(decoded - t) / T)[-1]
    return dict(accuracy=accuracy / T, rank=rank / T, error=error / T)


# Reconstructions with at most this many archetypes are decoded from their
# low-rank factors; larger ones (the full model) use the dense path.
LOWRANK_MAX_M = 10


def _lowrank_factors(sub, analysis_type, idx):
    """Factors (P, W) with reconstruct(sub, analysis_type, idx) == P @ W.T.

    The reconstruction is P0 @ Q.T with P0 (time x m) and Q (node x m), so it
    has rank m. Column-wise z-scoring subtracts the column mean (a rank-1 term
    built from the mean of P0) and divides each column by its std, so the
    z-scored matrix is still rank m: P = P0 - mean(P0), W = Q / (std + 1e-8).
    """
    sXC = np.asarray(sub["sXC"], dtype=float)[:, idx]
    S = np.asarray(sub["S"], dtype=float)[idx, :]
    if analysis_type == "spatial":
        P0, Q = sXC, S.T
    elif analysis_type == "temporal":
        P0, Q = S.T, sXC
    else:
        raise ValueError(analysis_type)
    P = P0 - P0.mean(axis=0, keepdims=True)
    if P.shape[1] == 1:
        sd = np.abs(Q[:, 0]) * P[:, 0].std()
    else:
        cov = P.T @ P / P.shape[0]
        sd = np.sqrt(np.clip(np.einsum("jm,mk,jk->j", Q, cov, Q), 0.0, None))
    W = Q / (sd + 1e-8)[:, None]
    return P, W


def _decode_subjects_lowrank(subjects, analysis_type, idx, nfolds, nreps, seed):
    """decode_subjects for a small archetype subset, without dense matrices.

    The fold average of the z-scored reconstructions is PA @ WA.T / n, where
    PA and WA stack the fold's subject factors side by side. The row
    correlations need only the row-centred W, so with G = Wc.T @ Wc
    (computed once for all subjects) the fold's correlation matrix is
    PA @ G[A, B] @ PB.T divided by the row norms sqrt(diag(PA G[A, A] PA.T)).
    The 1/n factors cancel in the correlation. Same values as the dense path
    up to rounding.
    """
    factors = [_lowrank_factors(s, analysis_type, idx) for s in subjects]
    m = factors[0][0].shape[1]
    P = np.concatenate([f[0] for f in factors], axis=1)  # time x (n*m)
    W = np.concatenate([f[1] for f in factors], axis=1)  # node x (n*m)
    Wc = W - W.mean(axis=0, keepdims=True)
    G = Wc.T @ Wc
    rng = np.random.default_rng(seed)
    rows = []
    for rep in range(int(nreps)):
        folds = _assignments(len(subjects), int(nfolds), rng)
        cols = np.repeat(folds, m)  # fold of each column of P / W
        for fold in range(int(nfolds)):
            ia = np.flatnonzero(cols == fold)
            ib = np.flatnonzero(cols != fold)
            PA, PB = P[:, ia], P[:, ib]
            num = PA @ (G[np.ix_(ia, ib)] @ PB.T)
            na = np.sqrt(np.einsum("tr,rt->t", PA, G[np.ix_(ia, ia)] @ PA.T))
            nb = np.sqrt(np.einsum("tr,rt->t", PB, G[np.ix_(ib, ib)] @ PB.T))
            with np.errstate(divide="ignore", invalid="ignore"):
                corrs = num / np.outer(na, nb)
            corrs = np.nan_to_num(corrs, nan=0.0, posinf=0.0, neginf=0.0)
            row = _decode_corr(corrs)
            row.update(rep=rep, fold=fold)
            rows.append(row)
    return pd.DataFrame(rows)


def decode_subjects(subjects, analysis_type, archetype_indices=None, nfolds=2, nreps=20, seed=42,
                    method="auto"):
    """Split-half time-point decoding, repeated over nreps random subject splits.

    method="auto" uses the low-rank path for subsets of up to LOWRANK_MAX_M
    archetypes and the dense path otherwise; "dense" and "lowrank" force one.
    The fold assignments (and so the seeds) are the same for both paths.
    """
    if method not in ("auto", "dense", "lowrank"):
        raise ValueError(method)
    if archetype_indices is not None and method != "dense":
        idx = np.asarray(list(archetype_indices), dtype=int)
        if method == "lowrank" or len(idx) <= LOWRANK_MAX_M:
            return _decode_subjects_lowrank(subjects, analysis_type, idx, nfolds, nreps, seed)
    stack = np.stack([reconstruct(s, analysis_type, archetype_indices) for s in subjects])
    rng = np.random.default_rng(seed)
    rows = []
    for rep in range(int(nreps)):
        folds = _assignments(len(subjects), int(nfolds), rng)
        for fold in range(int(nfolds)):
            a = stack[folds == fold].mean(axis=0)
            b = stack[folds != fold].mean(axis=0)
            corrs = _row_corr(a, b)
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
                   *, nfolds=2, nreps=20, seed=42):
    """Run/resume full and single-archetype decoding only.

    IMPORTANT: this function intentionally does *not* perform top-m selection.
    The canonical sparse analysis is the participant-held-out cross-fit in
    ``msaa_pipeline.heldout``. Keeping top-m out of this stage prevents the
    historical same-data ranking/evaluation procedure from being regenerated.

    CSVs are checkpointed after each condition so a long full sweep can be
    resumed safely. Existing keyed rows are not recomputed.
    """
    outdir = output_root / f"msaa_condrank_decoding_outputs_{analysis_type}_{fit_scope}"
    outdir.mkdir(parents=True, exist_ok=True)
    full_path = outdir / "full_summary.csv"
    arch_path = outdir / "per_archetype_mean_accuracy.csv"

    full = _load_csv(full_path)
    arch = _load_csv(arch_path)
    full_keys = _keyset(full, ["analysis_type","fit_scope","condition","K"])
    arch_keys = _keyset(arch, ["analysis_type","fit_scope","condition","K","archetype"])

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

    return _load_csv(full_path), _load_csv(arch_path)
