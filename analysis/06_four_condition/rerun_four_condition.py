"""Four-condition (paragraph-scrambled) sensitivity analysis with the current pipeline.

Redoes notebooks 01_fit_four_condition and 02_decode_four_condition (full-model part) with the
conventions of the October rerun:
  - fits: msaa_pipeline.fit.fit_one, FurthestSum picks computed up front (seed 42, the same picks
    multi_subject_aa would draw) and passed as init_idx, no post-fit sort (sequential numbering);
    the picks are stored in each fit and in fit_init_picks.csv;
  - decoding: the pipeline's dense split-half decoder (_assignments, _row_corr, _decode_corr) with a
    1e-8 tie tolerance, 2 folds x 20 repetitions, seed 42 (40 rows per summary, as in August).
Same data handling as August: paragraph source column 3 (duplicate) removed for the 700-factor data,
rest trimmed by 100 TRs, hypertools "across" normalisation within each condition; same 15-ratio grid
(spatial K = round(r x 700), temporal K = round(r x T_c)); PCA = node-space PCA fit on the pooled
condition data and reconstructed to 700 nodes. The same-data top-m part of notebook 02 is not redone.

Usage (repo root, BLAS threads set by the caller, e.g. OMP_NUM_THREADS=2):
    python analysis/06_four_condition/rerun_four_condition.py time  --out DIR [--fit spatial:intact:700]
    python analysis/06_four_condition/rerun_four_condition.py all   --out DIR [--workers 8]
Writes DIR/fits/*.npz and DIR/four_condition_ratio_decoding_outputs/{raw_htfa_summary,
pca_summary_with_retention, aa_full_summary_with_retention}.csv (same columns as August).
Pings ~/.discord_webhook (if present) at every 10% of each stage and at the end.
"""
import argparse
import json
import subprocess
import sys
import time
import traceback
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.io import loadmat

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
from msaa_pipeline.data import normalize_subject_list  # noqa: E402
from msaa_pipeline import decode as D  # noqa: E402
from msaa_pipeline.fit import fit_one, furthest_sum_picks, save_fit  # noqa: E402

CONDITIONS = ["intact", "paragraph", "word", "rest"]
PARAGRAPH_DUPLICATE_INDEX = 3  # 700-factor HTFA, as in August
TRIM_REST_BY = 100
RATIOS = [0.01, 0.02, 0.03, 0.05, 0.075, 0.10, 0.125, 0.15, 0.20, 0.25, 0.30, 0.40, 0.50, 0.75, 1.00]
NFOLDS, NREPS, SEED, TIE_TOL = 2, 20, 42, 1e-8
HOOK_FILE = Path.home() / ".discord_webhook"


def ping(msg, log=None):
    line = time.strftime("%F %T ") + msg
    print(line, flush=True)
    if log is not None:
        with open(log, "a") as f:
            f.write(line + "\n")
    if HOOK_FILE.exists():
        hook = HOOK_FILE.read_text().strip()
        subprocess.run(["curl", "-sS", "-o", "/dev/null", "-H", "Content-Type: application/json",
                        "-d", json.dumps({"content": msg}), hook], check=False)


_DATA = None


def load_data():
    global _DATA
    if _DATA is None:
        mat = loadmat(ROOT / "data" / "input" / "pieman_data.mat")
        out = {}
        for c in CONDITIONS:
            arr = mat[c]
            keep = [i for i in range(arr.shape[1]) if not (c == "paragraph" and i == PARAGRAPH_DUPLICATE_INDEX)]
            subj = [np.asarray(arr[:, i][0], dtype=float) for i in keep]
            if c == "rest":
                subj = [x[:-TRIM_REST_BY] for x in subj]
            out[c] = normalize_subject_list(subj, ROOT)
        _DATA = out
    return _DATA


def ratio_table():
    data = load_data()
    rows = []
    for c in CONDITIONS:
        T = data[c][0].shape[0]
        for r in RATIOS:
            for a, den in (("spatial", 700), ("temporal", T)):
                K = int(np.clip(np.rint(r * den), 1, den))
                rows.append(dict(condition=c, target_ratio=r, analysis_type=a, denominator=den, K=K,
                                 actual_ratio=K / den))
    t = pd.DataFrame(rows)
    return t.sort_values(["analysis_type", "condition", "target_ratio"]).drop_duplicates(
        ["analysis_type", "condition", "K"]).reset_index(drop=True)


def fit_path(fitdir, a, c, K):
    return fitdir / f"{a}AA_within_{c}_K{K}.npz"


def fit_task(fitdir, a, c, K):
    out = fit_path(fitdir, a, c, K)
    if out.exists():
        return dict(analysis_type=a, condition=c, K=K, status="skip")
    subj = load_data()[c]
    t = time.time()
    picks = furthest_sum_picks(subj, K, a, ROOT, seed=SEED)
    res = fit_one(subj, K, a, ROOT, maxiter=200, seed=SEED, verbose=False, init_idx=picks, sort_crit="none")
    tmp = out.with_name(out.stem + ".partial.npz")
    save_fit(tmp, res, [c] * len(subj), a, "within", np.zeros((0, 3)), np.zeros(0))
    tmp.rename(out)
    return dict(analysis_type=a, condition=c, K=K, status="fit", varexpl=float(res["varexpl"]),
                secs=float(res["secs"]), wall=time.time() - t, init_index=" ".join(map(str, picks)))


def decode_stack(stack):
    rng = np.random.default_rng(SEED)
    rows = []
    for rep in range(NREPS):
        folds = D._assignments(len(stack), NFOLDS, rng)
        for fold in range(NFOLDS):
            corrs = D._row_corr(stack[folds == fold].mean(axis=0), stack[folds != fold].mean(axis=0))
            corrs = np.nan_to_num(corrs, nan=0.0, posinf=0.0, neginf=0.0)
            rows.append(D._decode_corr(corrs, tie_tol=TIE_TOL))
    return D.summarize(pd.DataFrame(rows))


def zscore(X):
    X = np.asarray(X, dtype=float)
    return (X - X.mean(axis=0, keepdims=True)) / (X.std(axis=0, keepdims=True) + 1e-8)


def decode_task(kind, fitdir, a, c, K):
    subj = load_data()[c]
    T = subj[0].shape[0]
    base = dict(condition=c, K=K, T=T, chance=1.0 / T)
    if kind == "raw":
        return dict(base, method="raw_htfa", V=700, target_ratio=1.0, actual_ratio=1.0,
                    **decode_stack(np.stack([zscore(x) for x in subj])))
    if kind == "pca":
        from sklearn.decomposition import PCA
        model = PCA(n_components=700, svd_solver="full", whiten=False, random_state=SEED).fit(np.vstack(subj))
        stack = np.stack([zscore(model.transform(x)[:, :K] @ model.components_[:K] + model.mean_) for x in subj])
        return dict(base, method="PCA", actual_ratio=K / 700,
                    variance_explained=float(model.explained_variance_ratio_[:K].sum()), **decode_stack(stack))
    d = np.load(fit_path(fitdir, a, c, K), allow_pickle=True)
    subs = d["results_subj"].tolist()
    stack = np.stack([D.reconstruct(s, a) for s in subs])
    den = 700 if a == "spatial" else T
    return dict(base, analysis_type=a, ratio_denominator=den, actual_ratio=K / den,
                fit_source=f"fits/{fit_path(fitdir, a, c, K).name}", varexpl=float(d["varexpl"]),
                **decode_stack(stack))


def run_pool(stage, fn, jobs, workers, log):
    out, n, step = [], len(jobs), max(1, len(jobs) // 10)
    t0 = time.time()
    ping(f"[four-cond] {stage}: starting {n} jobs on {workers} workers", log)
    with ProcessPoolExecutor(max_workers=workers) as ex:
        futs = {ex.submit(fn, *j): j for j in jobs}
        for i, f in enumerate(as_completed(futs), 1):
            try:
                out.append(f.result())
            except Exception:
                ping(f"[four-cond] FAILED {stage} {futs[f]}: {traceback.format_exc().splitlines()[-1]}", log)
                raise
            if i % step == 0 or i == n:
                el = time.time() - t0
                ping(f"[four-cond] {stage} {i}/{n} ({100 * i // n}%) elapsed {el / 60:.1f} min, "
                     f"eta {el / i * (n - i) / 60:.1f} min", log)
    return out


def add_retention(df, raw):
    df = df.copy()
    rm = df.condition.map(raw.set_index("condition")["mean"])
    df["above_chance"] = df["mean"] - df["chance"]
    df["raw_mean"] = rm
    df["raw_above_chance"] = rm - df["chance"]
    df["fraction_raw_above_chance_retained"] = df["above_chance"] / df["raw_above_chance"]
    return df


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("cmd", choices=["time", "all"])
    ap.add_argument("--out", required=True)
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--fit", default="spatial:intact:700")
    a = ap.parse_args()
    out = Path(a.out).expanduser().resolve()
    fitdir = out / "fits"
    tabdir = out / "four_condition_ratio_decoding_outputs"
    fitdir.mkdir(parents=True, exist_ok=True); tabdir.mkdir(parents=True, exist_ok=True)
    log = out / "pings.log"

    if a.cmd == "time":
        at, c, K = a.fit.split(":")
        t = time.time(); r = fit_task(fitdir, at, c, int(K)); tf = time.time() - t
        t = time.time(); decode_task("aa", fitdir, at, c, int(K)); td = time.time() - t
        print(f"timed {a.fit}: fit {tf:.1f} s (status {r['status']}), decode {td:.1f} s")
        return

    t0 = time.time()
    try:
        rt = ratio_table()
        rt.to_csv(out / "ratio_K_map_unique_fits.csv", index=False)
        # largest fits first so the pool's tail is short
        jobs = sorted(((r.analysis_type, r.condition, int(r.K)) for r in rt.itertuples()),
                      key=lambda j: -j[2] * (700 if j[0] == "spatial" else 300))
        fits = run_pool("fit", fit_task, [(fitdir, *j) for j in jobs], a.workers, log)
        pd.DataFrame(fits).to_csv(out / "fit_log.csv", index=False)
        init = [dict(analysis_type=f["analysis_type"], condition=f["condition"], K=f["K"],
                     init_index=f.get("init_index", "")) for f in fits]
        pd.DataFrame(init).to_csv(out / "fit_init_picks.csv", index=False)

        sp = rt[rt.analysis_type == "spatial"]
        djobs = ([("raw", fitdir, None, c, 700) for c in CONDITIONS]
                 + [("pca", fitdir, None, r.condition, int(r.K)) for r in sp.itertuples()]
                 + [("aa", fitdir, r.analysis_type, r.condition, int(r.K)) for r in rt.itertuples()])
        res = run_pool("decode", decode_task, djobs, a.workers, log)
        raw = pd.DataFrame([r for r in res if r.get("method") == "raw_htfa"])
        pca = pd.DataFrame([r for r in res if r.get("method") == "PCA"])
        aa = pd.DataFrame([r for r in res if "analysis_type" in r])
        tr = rt.set_index(["analysis_type", "condition", "K"])["target_ratio"]
        pca["target_ratio"] = [tr[("spatial", c, k)] for c, k in zip(pca.condition, pca.K)]
        aa["target_ratio"] = [tr[(x, c, k)] for x, c, k in zip(aa.analysis_type, aa.condition, aa.K)]
        order = {c: i for i, c in enumerate(CONDITIONS)}
        raw.sort_values("condition", key=lambda s: s.map(order)).to_csv(tabdir / "raw_htfa_summary.csv", index=False)
        add_retention(pca, raw).sort_values(["condition", "K"]).to_csv(
            tabdir / "pca_summary_with_retention.csv", index=False)
        add_retention(aa, raw).sort_values(["analysis_type", "condition", "K"]).to_csv(
            tabdir / "aa_full_summary_with_retention.csv", index=False)
        ping(f"[four-cond] DONE in {(time.time() - t0) / 60:.1f} min: {len(fits)} fits, {len(res)} decodes -> {out}", log)
    except Exception as e:
        ping(f"[four-cond] FAILED after {(time.time() - t0) / 60:.1f} min: {type(e).__name__}: {e}", log)
        raise


if __name__ == "__main__":
    main()
