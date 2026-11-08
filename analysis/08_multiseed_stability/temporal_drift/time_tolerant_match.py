#!/usr/bin/env python3
"""Time-tolerant comparison of temporal archetypes between two runs (e.g. archbox vs Hellgate).

Why: each temporal archetype is built from a few time points (the columns of C, 300 x K).
fMRI is smooth in time, so neighbouring TRs are near-identical brain states and a fit can put
an archetype's weight on TR 101 instead of TR 100 at almost no cost. compare_hellgate.py
correlates C column j with column j time point by time point, so such a shift looks like a
different archetype. Here C is smoothed along time with a Gaussian (sigma in TRs) before
correlating, which forgives shifts of about one sigma but not moves to another moment.

Per temporal archetype j it reports:
  C_corr                raw corr of column j (what compare_hellgate.py calls C_corr)
  C_corr_s1/_s2/_s3     corr after smoothing both columns with sigma = 1, 2, 3 TRs
  com_shift_TR          shift of the weighted mean time point of column j, in TRs
  null_s{SIGMA}         best smoothed corr of A's j with any OTHER archetype k != j of B: the
                        baseline that smoothing alone gives; a tolerant match should beat it
  Smean_corr            raw corr of the subject-mean S row j (how the archetype is used over time)
  status                same     raw C_corr >= --close
                        shifted  raw < --close, but smoothed (sigma = --sigma) >= --tol and > null
                        changed  otherwise
Per fit: counts per status, and whether the one-to-one best match on smoothed C is the identity.

  python time_tolerant_match.py --a ~/aa/msaa-rerun-seqnum-2026-10-27 \
                                --b ~/aa/hellgate-setup/results/hellgate_fingerprint.npz
Each side is a run root (with msaa_flexible_outputs_npz/), a folder holding
hellgate_fingerprint.npz, or a fingerprint .npz from hellgate-setup/fingerprint.py.
Writes archetypes_tolerant.csv, fits_tolerant.csv and summary.txt to --out (default: results/
next to this script). numpy + pandas only; scipy is used for the one-to-one match if present.
"""
import argparse
import sys
from pathlib import Path
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))  # fingerprint.py sits next to this file
import fingerprint as fpmod  # noqa: E402

SIGMAS = (1, 2, 3)


def load(path):
    p = Path(path).expanduser()
    if p.is_file():
        fp = fpmod.load(p)
    elif (p / "hellgate_fingerprint.npz").exists():
        fp = fpmod.load(p / "hellgate_fingerprint.npz")
    else:
        names = {f.stem for f in (p / "msaa_flexible_outputs_npz").glob("temporal*.npz") if ".partial" not in f.name}
        print(f"reading {len(names)} temporal fits from {p} ...", flush=True)
        fp = fpmod.build(p, names)
    return {k: v for k, v in fp.items() if k.startswith("temporal")}


def smooth(C, sigma):
    """Gaussian smoothing of each column along time (rows), edges reflected."""
    h = int(4 * sigma)
    k = np.exp(-np.arange(-h, h + 1) ** 2 / (2 * sigma ** 2))
    k /= k.sum()
    P = np.pad(C, ((h, h), (0, 0)), mode="reflect")
    return np.stack([np.convolve(P[:, j], k, "valid") for j in range(C.shape[1])], axis=1)


def corr_matrix(A, B):
    A = A - A.mean(0); B = B - B.mean(0)
    na, nb = np.linalg.norm(A, axis=0), np.linalg.norm(B, axis=0)
    R = (A / np.where(na > 0, na, 1)).T @ (B / np.where(nb > 0, nb, 1))
    R[np.outer(na == 0, nb == 0)] = 1.0
    return R


def one_to_one(R):
    try:
        from scipy.optimize import linear_sum_assignment
        r, c = linear_sum_assignment(-R)
        m = np.empty(R.shape[0], int); m[r] = c
        return m
    except ImportError:  # greedy fallback: take the best remaining pair each time
        m = -np.ones(R.shape[0], int); R = R.copy()
        for _ in range(R.shape[0]):
            i, j = np.unravel_index(np.argmax(R), R.shape)
            m[i] = j; R[i, :] = -np.inf; R[:, j] = -np.inf
        return m


def compare_fit(name, a, b, sigma, close, tol):
    Ca, Cb = np.asarray(a["C"], float), np.asarray(b["C"], float)
    K = Ca.shape[1]
    t = np.arange(Ca.shape[0])[:, None]
    com = lambda C: (C * t).sum(0) / np.where(C.sum(0) > 0, C.sum(0), 1)
    out = dict(fit=name, archetype=np.arange(K), C_corr=np.diag(corr_matrix(Ca, Cb)))
    Rs = {}
    for s in SIGMAS:
        Rs[s] = corr_matrix(smooth(Ca, s), smooth(Cb, s))
        out[f"C_corr_s{s}"] = np.diag(Rs[s])
    R = Rs[sigma].copy()
    np.fill_diagonal(R, -np.inf)
    out[f"null_s{sigma}"] = R.max(1) if K > 1 else np.full(K, np.nan)
    out["com_shift_TR"] = np.abs(com(Ca) - com(Cb))
    out["Smean_corr"] = np.diag(corr_matrix(np.asarray(a["Smean"], float).T, np.asarray(b["Smean"], float).T))
    A = pd.DataFrame(out)
    tol_ok = (A[f"C_corr_s{sigma}"] >= tol) & (A[f"C_corr_s{sigma}"] > A[f"null_s{sigma}"])
    A["status"] = np.where(A.C_corr >= close, "same", np.where(tol_ok, "shifted", "changed"))
    m = one_to_one(Rs[sigma])
    F = dict(fit=name, K=K, n_same=int((A.status == "same").sum()), n_shifted=int((A.status == "shifted").sum()),
             n_changed=int((A.status == "changed").sum()), tolerant_numbering_same=bool((m == np.arange(K)).all()),
             n_moved_tolerant=int((m != np.arange(K)).sum()), min_C_corr=A.C_corr.min(),
             min_C_corr_tolerant=A[f"C_corr_s{sigma}"].min(), max_com_shift_TR=A.com_shift_TR.max())
    return A, F


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--a", default=str(Path.home() / "aa/msaa-rerun-seqnum-2026-10-27"))
    ap.add_argument("--b", default=str(Path.home() / "aa/hellgate-setup/results/hellgate_fingerprint.npz"))
    ap.add_argument("--out", default=str(HERE / "results"))
    ap.add_argument("--sigma", type=int, default=2, choices=SIGMAS, help="tolerance in TRs used for status")
    ap.add_argument("--close", type=float, default=0.99, help="raw corr counted as the same archetype")
    ap.add_argument("--tol", type=float, default=0.9, help="smoothed corr counted as a shifted copy")
    a = ap.parse_args()

    fa, fb = load(a.a), load(a.b)
    names = sorted(set(fa) & set(fb))
    print(f"temporal fits: {len(fa)} in A, {len(fb)} in B, {len(names)} compared", flush=True)
    arch, fits = [], []
    for n in names:
        if np.shape(fa[n]["C"]) != np.shape(fb[n]["C"]):
            print(f"skip {n}: C shape {np.shape(fa[n]['C'])} vs {np.shape(fb[n]['C'])}")
            continue
        A, F = compare_fit(n, fa[n], fb[n], a.sigma, a.close, a.tol)
        arch.append(A); fits.append(F)
    A, F = pd.concat(arch, ignore_index=True), pd.DataFrame(fits)
    A["scope"] = A.fit.str.extract(r"^temporalAA_(across|within)")[0]
    F["scope"] = F.fit.str.extract(r"^temporalAA_(across|within)")[0]

    out = Path(a.out).expanduser(); out.mkdir(parents=True, exist_ok=True)
    A.to_csv(out / "archetypes_tolerant.csv", index=False)
    F.to_csv(out / "fits_tolerant.csv", index=False)

    s = f"s{a.sigma}"
    L = [f"A: {a.a}", f"B: {a.b}",
         f"temporal fits compared: {len(F)}, archetypes: {len(A)}",
         f"status uses smoothing sigma = {a.sigma} TR, same if raw C corr >= {a.close}, "
         f"shifted if smoothed corr >= {a.tol} and above the best other archetype", "",
         "Archetypes by status:", A.status.value_counts().reindex(["same", "shifted", "changed"]).fillna(0)
         .astype(int).to_string(), "",
         "By scope and K band (counts):",
         A.assign(K_band=pd.cut(A.groupby("fit").archetype.transform("size"), [0, 40, 80, 150, 400]))
          .groupby(["scope", "K_band"], observed=True).status.value_counts().unstack(fill_value=0).to_string(), "",
         "Archetypes that are not 'same' (quantiles):",
         A[A.status != "same"][["C_corr", "C_corr_s1", "C_corr_s2", "C_corr_s3", f"null_{s}", "com_shift_TR",
                                "Smean_corr"]].describe(percentiles=[.1, .5, .9]).T
          [["count", "min", "10%", "50%", "90%", "max"]].to_string(), "",
         f"Numbering on smoothed C (sigma = {a.sigma}): identity in {int(F.tolerant_numbering_same.sum())}/{len(F)} fits", ""]
    ch = A[A.status == "changed"].sort_values(f"C_corr_{s}")
    L += [f"'changed' archetypes ({len(ch)}):",
          ch[["fit", "archetype", "C_corr", f"C_corr_{s}", f"null_{s}", "com_shift_TR"]].to_string(index=False)
          if len(ch) else "  none"]
    text = "\n".join(L)
    (out / "summary.txt").write_text(text + "\n")
    print(text)
    print(f"\nwrote {out}/archetypes_tolerant.csv, fits_tolerant.csv, summary.txt")


if __name__ == "__main__":
    main()
