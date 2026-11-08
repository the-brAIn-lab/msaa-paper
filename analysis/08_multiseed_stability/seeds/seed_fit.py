#!/usr/bin/env python3
"""One fit of the multi-seed run: line N of tasks/<STAGE>_fit.txt ("seed at scope K cond").
Runs the seqnum code's own fit_task inside that seed's directory, so the fit is exactly what
`rerun/run_parallel.py fit` does (no post-fit sort, 200 iterations), except that the frozen
FurthestSum picks it reads were drawn with this seed. Then writes the slim copy and meta/<fit>.json.

  python seed_fit.py STAGE N
"""
import json, os, platform, sys, time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from seeds_common import fit_name, seed_dir, write_slim  # noqa: E402

HERE = Path(__file__).resolve().parent


def main():
    stage, n = sys.argv[1], int(sys.argv[2])
    s, at, sc, K, cond = (HERE / "tasks" / f"{stage}_fit.txt").read_text().splitlines()[n - 1].split()
    s, K, cond = int(s), int(K), (None if cond == "-" else cond)
    d = seed_dir(stage, s)
    name = fit_name(at, sc, K, cond)
    slim = d / "slim" / name
    if slim.exists():
        print("skip", s, name); return
    threads = int(os.environ.get("MSAA_THREADS", "1"))
    os.environ.setdefault("MSAA_TIE_TOL", "1e-08")
    os.chdir(d)
    sys.path.insert(0, str(d / "rerun"))
    import run_parallel as rp          # this seed directory's copy
    rp._init(threads)
    os.chdir(d)
    t = time.time()
    msg = rp.fit_task(at, sc, K, cond, json.loads((d / "configs" / "profiles.json").read_text())["full"]["fit_maxiter"])
    rp._log(f"seed {s} {msg}")
    fit = d / "msaa_flexible_outputs_npz" / name
    if not fit.exists():
        raise SystemExit(f"fit_task did not write {fit}")
    write_slim(fit, slim)
    (d / "meta").mkdir(exist_ok=True)
    (d / "meta" / name.replace(".npz", ".json")).write_text(json.dumps(dict(
        fit=name, seed=s, host=platform.node(), threads=threads, wall_s=round(time.time() - t, 1))))
    print(f"done seed {s} {name} {time.time() - t:.0f}s", flush=True)


if __name__ == "__main__":
    main()
