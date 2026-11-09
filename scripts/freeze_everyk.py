"""freeze_init.py logic for the every-K grid, with the inputs loaded once per worker.
Uses msaa_pipeline.fit.furthest_sum_picks exactly as scripts/freeze_init.py does; existing rows kept.
Usage: python scripts/freeze_everyk.py REPO_ROOT OUT.csv.gz"""
import os, sys, time
for v in ["OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"]: os.environ[v] = "1"
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor
import pandas as pd
ROOT = Path(sys.argv[1]); OUT = Path(sys.argv[2])
sys.path.insert(0, str(ROOT / "src"))
KEY = ["analysis_type", "fit_scope", "condition", "K"]
GRID = {"spatial": range(2, 701), "temporal": range(2, 301)}
D = None
def init():
    global D
    from msaa_pipeline.data import load_inputs
    D, _, _ = load_inputs(ROOT, trim_rest_by=100, max_subjects_per_condition=None)
def picks(at, sc, cond, K):
    from msaa_pipeline.fit import furthest_sum_picks
    subjects = [x for c in ["intact", "word", "rest"] for x in D[c]] if sc == "across" else D[cond]
    idx = furthest_sum_picks(subjects, K, at, ROOT, seed=42)
    return pd.DataFrame(dict(analysis_type=at, fit_scope=sc, condition=cond or "acrossCond", K=K,
                             archetype=range(K), init_index=idx))
if __name__ == "__main__":
    old = pd.read_csv(ROOT / "configs" / "furthestsum_init.csv.gz")
    have = set(map(tuple, old[KEY].astype(str).drop_duplicates().itertuples(index=False, name=None)))
    todo = [(at, sc, c, K) for at, ks in GRID.items() for sc in ["across", "within"] for K in ks
            for c in ([None] if sc == "across" else ["intact", "word", "rest"])
            if (at, sc, c or "acrossCond", str(K)) not in have]
    todo.sort(key=lambda t: -t[3])
    print(len(todo), "fits need picks;", len(have), "frozen", flush=True)
    t = time.time(); new = []
    with ProcessPoolExecutor(int(sys.argv[3]), initializer=init) as ex:
        for i, df in enumerate(ex.map(picks, *zip(*todo), chunksize=4), 1):
            new.append(df)
            if i % 200 == 0: print(i, f"{time.time()-t:.0f}s", flush=True)
    df = pd.concat([old, *new], ignore_index=True).sort_values(KEY + ["archetype"]).reset_index(drop=True)
    df.to_csv(OUT, index=False)
    print("wrote", OUT, df.groupby(KEY).ngroups, "fits", len(df), "rows")
