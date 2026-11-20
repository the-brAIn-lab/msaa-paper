"""Parity of the low-rank decoding path with the dense path (synthetic data, seconds).

decode_subjects() decodes reconstructions built from a small archetype subset
(m <= LOWRANK_MAX_M) from their low-rank factors instead of averaging dense
time x node matrices (see the note in docs/REPRODUCE.md). These checks build
random MS-AA-shaped subjects and confirm, for both orientations:

  1. _lowrank_factors(): P @ W.T equals the dense z-scored reconstruction.
  2. decode_subjects(method="lowrank") gives the same accuracy, rank and error
     as decode_subjects(method="dense") on the same folds and seeds.

Run with  python tests/test_lowrank_parity.py  (or pytest, or make test).
No input data or fitted models are needed.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from msaa_pipeline import decode  # noqa: E402

T, NODES, K, N_SUBJECTS = 60, 80, 12, 10
SUBSETS = ([3], [0, 5, 7], [1, 2, 4, 8, 11], list(range(10)))


def _subjects(analysis_type, seed=0):
    """Random subjects with the shapes reconstruct() expects.

    spatial:  sXC is time x K, S is K x node (columns on the simplex).
    temporal: sXC is node x K, S is K x time, mostly exact zeros like real temporal S.
    """
    rng = np.random.default_rng(seed)
    rows, cols = (T, NODES) if analysis_type == "spatial" else (NODES, T)
    out = []
    for _ in range(N_SUBJECTS):
        sXC = rng.standard_normal((rows, K))
        S = rng.gamma(0.5, size=(K, cols))
        if analysis_type == "temporal":
            S[rng.random(S.shape) < 0.7] = 0.0
        S /= S.sum(axis=0, keepdims=True) + 1e-12
        out.append({"sXC": sXC, "S": S})
    return out


def test_factors_match_dense_reconstruction():
    for analysis_type in ("spatial", "temporal"):
        for sub in _subjects(analysis_type):
            for idx in SUBSETS:
                dense = decode.reconstruct(sub, analysis_type, idx)
                P, W = decode._lowrank_factors(sub, analysis_type, np.asarray(idx))
                assert P.shape[1] == len(idx)
                np.testing.assert_allclose(P @ W.T, dense, rtol=1e-10, atol=1e-10,
                                           err_msg=f"{analysis_type} m={len(idx)}")


def test_lowrank_decoding_matches_dense():
    for analysis_type in ("spatial", "temporal"):
        subjects = _subjects(analysis_type, seed=1)
        for idx in SUBSETS:
            kw = dict(nfolds=2, nreps=5, seed=42)
            lr = decode.decode_subjects(subjects, analysis_type, idx, method="lowrank", **kw)
            de = decode.decode_subjects(subjects, analysis_type, idx, method="dense", **kw)
            assert list(lr[["rep", "fold"]].itertuples(index=False)) == \
                list(de[["rep", "fold"]].itertuples(index=False))
            for col in ("accuracy", "rank", "error"):
                np.testing.assert_allclose(lr[col].to_numpy(), de[col].to_numpy(), rtol=0, atol=1e-12,
                                           err_msg=f"{analysis_type} m={len(idx)} {col}")


def test_auto_uses_lowrank_only_for_small_subsets():
    assert decode.LOWRANK_MAX_M == 10
    subjects = _subjects("spatial", seed=2)
    kw = dict(nfolds=2, nreps=2, seed=7)
    small = list(range(decode.LOWRANK_MAX_M))
    auto = decode.decode_subjects(subjects, "spatial", small, **kw)
    lr = decode.decode_subjects(subjects, "spatial", small, method="lowrank", **kw)
    assert auto.equals(lr)


if __name__ == "__main__":
    tests = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    for t in tests:
        t()
        print(f"PASS {t.__name__}")
    print(f"low-rank parity: {len(tests)} checks passed")
