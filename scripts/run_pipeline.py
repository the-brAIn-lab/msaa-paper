#!/usr/bin/env python3
from __future__ import annotations
import argparse, json, os, sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))

from msaa_pipeline.profiles import load_profile, parse_job_key
from msaa_pipeline.fit import run_fit_job
from msaa_pipeline.decode import run_decode_job
from msaa_pipeline.heldout import run_heldout_topm_job
from msaa_pipeline.qc import build_qc
import msaa_pipeline.decode as decode_mod


def main():
    ap = argparse.ArgumentParser(description="Run the cleaned MS-AA fit -> decode -> held-out top-m -> QC pipeline.")
    ap.add_argument("--profile", default="smoke", choices=["smoke", "mini-paper", "full"])
    ap.add_argument("--output-root", default=None, help="Output directory. Defaults to outputs/<profile>; for canonical full outputs use '.'.")
    ap.add_argument("--fit-only", action="store_true")
    ap.add_argument("--decode-only", action="store_true")
    ap.add_argument("--heldout-only", action="store_true", help="Run only the canonical held-out top-m stage; requires fits to exist.")
    ap.add_argument("--force", action="store_true", help="Refit existing model files.")
    ap.add_argument("--numbering", default="auto", choices=["auto", "frozen", "legacy"],
                    help="frozen: frozen FurthestSum picks + no post-fit sort (the paper's canonical numbering, "
                         "as scripts/run_parallel.py); legacy: recompute FurthestSum and sort after fitting "
                         "(the original outputs' numbering). auto = frozen for the full profile, legacy for "
                         "smoke and mini-paper, which are checked against the legacy outputs.")
    ap.add_argument("--tie-tol", type=float, default=None,
                    help="decode correlation tie tolerance; default 1e-8 for the full profile (as "
                         "scripts/run_parallel.py), 0 (exact comparison, legacy) otherwise.")
    args = ap.parse_args()

    profile = load_profile(REPO, args.profile)
    output_root = Path(args.output_root).resolve() if args.output_root else (REPO / "outputs" / args.profile.replace("-", "_"))
    output_root.mkdir(parents=True, exist_ok=True)
    print(f"Profile: {args.profile}\nOutput root: {output_root}\n{profile['description']}")
    full = args.profile == "full"
    frozen = args.numbering == "frozen" or (args.numbering == "auto" and full)
    if args.tie_tol is not None:
        decode_mod.TIE_TOL = float(args.tie_tol)
    elif "MSAA_TIE_TOL" not in os.environ:
        decode_mod.TIE_TOL = 1e-8 if full else 0.0
    print(f"Archetype numbering: {'frozen FurthestSum picks (canonical)' if frozen else 'legacy (FurthestSum + sort)'}; "
          f"decode tie tolerance {decode_mod.TIE_TOL:g}")

    manifests = []
    if not args.decode_only and not args.heldout_only:
        for key, ks in profile["jobs"].items():
            analysis_type, fit_scope = parse_job_key(key)
            print("\n=== FIT", key, ks, "===")
            m = run_fit_job(
                REPO, output_root, analysis_type, fit_scope, ks,
                max_subjects_per_condition=profile["max_subjects_per_condition"],
                maxiter=profile["fit_maxiter"], seed=42, skip_existing=not args.force,
                frozen_init=frozen,
            )
            manifests.append(m)
        if manifests:
            import pandas as pd
            pd.concat(manifests, ignore_index=True).to_csv(output_root / "fit_manifest.csv", index=False)

    if args.fit_only:
        return

    if not args.heldout_only:
        for key, ks in profile["jobs"].items():
            analysis_type, fit_scope = parse_job_key(key)
            print("\n=== FULL + SINGLE-ARCHETYPE DECODE", key, ks, "===")
            run_decode_job(
                output_root, analysis_type, fit_scope, ks,
                nfolds=profile["decode_nfolds"],
                nreps=profile["decode_nreps"], seed=42,
            )

    # Sparse top-m is intentionally *only* the participant-held-out cross-fit.
    for key, ks in profile.get("heldout_jobs", {}).items():
        analysis_type, fit_scope = parse_job_key(key)
        if fit_scope != "across":
            raise ValueError(f"Held-out top-m is only defined for across-condition fits, got {key}")
        print("\n=== HELD-OUT TOP-M", key, ks, "===")
        run_heldout_topm_job(
            output_root, analysis_type, ks,
            top_m=profile["top_m"],
            n_outer_reps=profile["heldout_outer_reps"],
            nfolds=profile["decode_nfolds"],
            nreps=profile["decode_nreps"], seed=42,
        )

    report = build_qc(output_root, args.profile)
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
