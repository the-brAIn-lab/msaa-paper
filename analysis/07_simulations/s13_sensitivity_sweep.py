"""Simulation sensitivity sweep for Supplementary Figure S13.

Varies one simulation parameter at a time around the metastable regime of
Figure 3 and records PCA vs MS-AA recovery. The generator, model fits and
recovery metrics are executed directly from
``analysis/07_simulations/revision_simulation.ipynb`` (the notebook that
writes ``pca_msaa_extended_results_paper.csv``, which Figure 3 reads), so
nothing here re-implements the simulation.

Design
------
* Eight panels: observation noise, between-subject heterogeneity, fitted K,
  motif overlap (distinctiveness scale), and the four dynamical-mixing
  components (mean dwell, transition width, secondary-motif weight, AR(1) rho).
* Five levels per panel; all non-varied parameters stay at the notebook's
  metastable defaults (noise 0.25, heterogeneity 1.0, K 6, distinctiveness
  1.0, dwell 12, transition width 4, secondary weight 0.25, rho 0.85).
* Twenty replicates. Replicate r uses the Figure 3 metastable seed
  SeedSequence([RNG_SEED, 303, r]) at every level of every panel (common
  random numbers). The default setting is therefore one shared set of 20
  simulations, identical to Figure 3's metastable regime, and it is fit once.

Run from a directory that contains ``helpers/`` and ``hypertools/`` (the
repository root, or the paper run directory that symlinks them):

    python analysis/07_simulations/s13_sensitivity_sweep.py --check
    python analysis/07_simulations/s13_sensitivity_sweep.py --workers 14 \
        --out simulation_pca_msaa_extended_outputs/s13_sensitivity_results.csv

``--check`` fits only the shared baseline and compares it with the metastable
rows of the Figure 3 results CSV (``--reference``).
"""

from __future__ import annotations

import argparse
import contextlib
import io
import json
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")

import numpy as np
import pandas as pd

NOTEBOOK = Path(__file__).resolve().with_name("revision_simulation.ipynb")
# Config, MS-AA wrapper, motif construction, normalisation + metrics,
# evaluate_methods, metastable generator, metastable diagnostics.
CELLS = [1, 4, 6, 9, 11, 17, 19]
METASTABLE_REGIME_CODE = 303  # seed code of the Figure 3 metastable regime
N_REPS = 20

DEFAULTS = dict(
    noise_sd=0.25,
    heterogeneity_scale=1.0,
    K_fit=6,
    distinctiveness_scale=1.0,
    mean_dwell=12,
    transition_width=4,
    secondary_weight=0.25,
    ar1_rho=0.85,
)

# Levels are listed from easiest to hardest for motif recovery, so that
# "more noise / overlap / mixing" always plots to the right.
LEVELS = {
    "noise_sd": [0.05, 0.10, 0.25, 0.40, 0.60],
    "heterogeneity_scale": [0.0, 0.5, 1.0, 1.5, 2.0],
    "K_fit": [4, 5, 6, 7, 8],
    "distinctiveness_scale": [1.15, 1.0, 0.85, 0.70, 0.55],
    "mean_dwell": [16, 14, 12, 10, 8],
    "transition_width": [2, 3, 4, 6, 8],
    "secondary_weight": [0.10, 0.175, 0.25, 0.35, 0.45],
    "ar1_rho": [0.75, 0.80, 0.85, 0.875, 0.90],
}
INT_PARAMS = {"K_fit", "mean_dwell", "transition_width"}

# All four mixing components at their most-mixed level together (the
# notebook's "more_mixed" preset), to test for an interaction. Reported in
# the caption, not plotted.
COMBINED_MIXING = dict(mean_dwell=8, transition_width=8, secondary_weight=0.45, ar1_rho=0.90)

_NS: dict | None = None


def load_ns() -> dict:
    """Execute the needed notebook cells once per process."""
    global _NS
    if _NS is not None:
        return _NS
    nb = json.loads(NOTEBOOK.read_text())
    prelude = (
        "import os, sys, time\nimport numpy as np\nimport pandas as pd\n"
        "import matplotlib\nmatplotlib.use('Agg')\nimport matplotlib.pyplot as plt\n"
        "plt.show = lambda *a, **k: plt.close('all')\n"
        "from scipy.optimize import linear_sum_assignment\n"
        "from sklearn.decomposition import PCA\n"
        "sys.path.insert(0, 'helpers')\n"
        "from MultiSubject_AA import Subject, multi_subject_aa\n"
        "display = print\n"
    )
    ns: dict = {"__name__": "revision_simulation"}
    with contextlib.redirect_stdout(io.StringIO()):
        exec(prelude, ns)
        for i in CELLS:
            cell = nb["cells"][i]
            assert cell["cell_type"] == "code", i
            src = "".join(cell["source"]).splitlines()
            exec("\n".join(l for l in src if not l.lstrip().startswith(("%", "!"))), ns)
    assert ns["RUN_MODE"] == "paper"
    assert (ns["N_SUBJECTS"], ns["T_PER_SUBJECT"]) == (16, 180)
    for key, nb_name in [("mean_dwell", "MEAN_DWELL"), ("transition_width", "TRANSITION_WIDTH"),
                         ("secondary_weight", "SECONDARY_WEIGHT"), ("ar1_rho", "AR1_RHO"),
                         ("noise_sd", "BASE_NOISE_SD"), ("K_fit", "N_TRUE_MOTIFS")]:
        assert ns[nb_name] == DEFAULTS[key], (nb_name, ns[nb_name], DEFAULTS[key])
    _NS = ns
    return ns


def rep_seed(ns: dict, rep: int) -> int:
    return int(np.random.SeedSequence([ns["RNG_SEED"], METASTABLE_REGIME_CODE, rep])
               .generate_state(1)[0])


def run_task(task: dict) -> dict:
    ns = load_ns()
    rep = task["rep"]
    cfg = dict(DEFAULTS)
    if task["param"] == "mixing_combined":
        cfg.update(COMBINED_MIXING)
    elif task["param"] != "baseline":
        cfg[task["param"]] = task["level"]
    seed = rep_seed(ns, rep)
    motifs = ns["make_overlap_variant"](ns["TRUE_MOTIFS_RAW"], cfg["distinctiveness_scale"])

    t0 = time.time()
    sim = ns["simulate_metastable"](
        seed,
        noise_sd=cfg["noise_sd"],
        heterogeneity_scale=cfg["heterogeneity_scale"],
        motifs_raw=motifs,
        mean_dwell=int(cfg["mean_dwell"]),
        transition_width=int(cfg["transition_width"]),
        secondary_weight=float(cfg["secondary_weight"]),
        ar1_rho=float(cfg["ar1_rho"]),
    )
    diag = ns["metastable_diagnostics"](sim)
    row = ns["evaluate_methods"](
        sim["subjects"], sim["true_patterns"], K_fit=int(cfg["K_fit"]), seed=seed,
        regime="s13_sensitivity", rep=rep, true_activations=sim["activations"],
        truth_kind=sim["truth_kind"],
        extra={
            "param": task["param"],
            "level": task["level"],
            "condition": f"{task['param']}={task['level']}",
            **{k: v for k, v in cfg.items() if k != "K_fit"},
            "mean_pairwise_cosine": ns["mean_pairwise_cosine"](motifs),
            **{f"diag_{k}": v for k, v in diag.items()},
        },
    )
    row["task_sec"] = time.time() - t0
    return row


def build_tasks(check_only: bool) -> list[dict]:
    tasks = [dict(param="baseline", level=np.nan, rep=r) for r in range(N_REPS)]
    if check_only:
        return tasks
    for p, levels in LEVELS.items():
        for lv in levels:
            if lv == DEFAULTS[p]:
                continue  # shared baseline
            lv = int(lv) if p in INT_PARAMS else float(lv)
            tasks += [dict(param=p, level=lv, rep=r) for r in range(N_REPS)]
    tasks += [dict(param="mixing_combined", level=1.0, rep=r) for r in range(N_REPS)]
    return tasks


def task_key(param, level, rep):
    return (str(param), "nan" if pd.isna(level) else f"{float(level):g}", int(rep))


def compare_with_reference(res: pd.DataFrame, ref_path: Path) -> bool:
    ref = pd.read_csv(ref_path)
    ref = ref[(ref["regime"] == "metastable") & (ref["condition"] == "base")]
    base = res[res["param"] == "baseline"]
    m = base.merge(ref, on="rep", suffixes=("", "_ref"))
    ok = len(m) == N_REPS and (m["seed"] == m["seed_ref"]).all()
    for col in ["aa_motif_recovery", "pca_motif_recovery"]:
        d = (m[col] - m[f"{col}_ref"]).abs().max()
        print(f"{col}: new mean {m[col].mean():.4f} | Figure 3 mean "
              f"{m[col + '_ref'].mean():.4f} | max |diff| per rep {d:.2e}")
        ok &= bool(d < 1e-6)
    print("seeds identical:", bool((m["seed"] == m["seed_ref"]).all()), "| n =", len(m))
    print("BASELINE CHECK", "PASSED" if ok else "FAILED")
    return ok


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", default="simulation_pca_msaa_extended_outputs/s13_sensitivity_results.csv")
    ap.add_argument("--reference",
                    default="simulation_pca_msaa_extended_outputs/pca_msaa_extended_results_paper.csv")
    ap.add_argument("--workers", type=int, default=14)
    ap.add_argument("--check", action="store_true", help="fit the shared baseline only and compare with Figure 3")
    a = ap.parse_args()

    out = Path(a.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    tasks = build_tasks(a.check)
    done = set()
    if out.exists():
        old = pd.read_csv(out)
        done = {task_key(p, l, r) for p, l, r in zip(old["param"], old["level"], old["rep"])}
    todo = [t for t in tasks if task_key(t["param"], t["level"], t["rep"]) not in done]
    print(f"{len(todo)} fits to run ({len(tasks) - len(todo)} already in {out})", flush=True)

    t0 = time.time()
    with ProcessPoolExecutor(max_workers=a.workers) as ex:
        futs = [ex.submit(run_task, t) for t in todo]
        for i, f in enumerate(as_completed(futs), 1):
            pd.DataFrame([f.result()]).to_csv(out, mode="a", header=not out.exists(), index=False)
            if i % 50 == 0 or i == len(todo):
                print(f"{i}/{len(todo)} fits after {time.time() - t0:.0f}s", flush=True)

    ok = compare_with_reference(pd.read_csv(out), Path(a.reference))
    if not ok:
        sys.exit(1)


if __name__ == "__main__":
    main()
