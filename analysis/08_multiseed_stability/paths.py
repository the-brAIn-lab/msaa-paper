"""Where the multi-seed analysis reads and writes, relative to the repo root.

The defaults match what run_local.py writes, so every script in this folder runs from a clone
without path arguments once the run (or the relevant steps of it) has finished:

  msaa_flexible_outputs_npz/                     seed-0 (canonical) fits, full .npz
  msaa_condrank_decoding_outputs_<analysis>_across/full_summary.csv   seed-0 decoding summary
  seeds/slim/seed00/, seeds/full/seedNN/slim/    compact fits of every seed (C, Smean, varexpl, ...)
  seeds/results/                                 cross-seed tables (per_K.csv, matching_*.csv, ...)
  heldout-everyk-<analysis>-<seed0|all>-results/ held-out top-m at every K
  heldout-seeds-results/                         held-out top-m at the 30 paper K
  outputs/multiseed_stability/                   default output folder of the scripts here

Set MSAA_ROOT to read from a run root other than the repo root.
"""
from __future__ import annotations

import os
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
ROOT = Path(os.environ.get("MSAA_ROOT", REPO)).expanduser().resolve()

FITS0 = ROOT / "msaa_flexible_outputs_npz"
SEEDS = ROOT / "seeds"
SEEDS_RESULTS = SEEDS / "results"
HELDOUT30 = ROOT / "heldout-seeds-results"
OUT = ROOT / "outputs" / "multiseed_stability"
DATA = REPO / "data" / "input" / "pieman_data.mat"
X_NPZ = OUT / "X.npz"                     # normalised scans, written by plots/loadx.py

# node table: node, x, y, z, width, network, code for the 700 HTFA nodes (scripts/build_node_yeo7.py)
NODES_CONFIG = REPO / "configs" / "node_yeo7.csv"
NODES_LEGACY = REPO / "reference_outputs" / "legacy_full" / "paper_figure_outputs" / "Figure5_spatial_maps_allnodes.csv"
NODES = NODES_CONFIG if NODES_CONFIG.exists() else NODES_LEGACY

N_SEEDS = 10


def decoding_summary(analysis: str) -> Path:
    return ROOT / f"msaa_condrank_decoding_outputs_{analysis}_across" / "full_summary.csv"


def heldout_everyk(analysis: str, which: str = "all") -> Path:
    """which = 'all' (seeds 0-9) or 'seed0'."""
    return ROOT / f"heldout-everyk-{analysis}-{which}-results"


def full_fit0(analysis: str, K: int) -> Path:
    return FITS0 / f"{analysis}AA_across_acrossCond_K{K}.npz"


def slim_fit(analysis: str, K: int, seed: int) -> Path:
    name = f"{analysis}AA_across_acrossCond_K{K}.npz"
    return SEEDS / "slim" / "seed00" / name if seed == 0 else SEEDS / "full" / f"seed{seed:02d}" / "slim" / name


def seed_fits(analysis: str, K: int, seeds=range(N_SEEDS), required: bool = True) -> dict[int, Path]:
    """{seed: compact across-condition fit} for one analysis and K."""
    out = {s: slim_fit(analysis, K, s) for s in seeds}
    missing = [s for s, p in out.items() if not p.exists()]
    if missing and required:
        raise SystemExit(f"missing {analysis} K={K} fits for seeds {missing} (expected e.g. {out[missing[0]]}); "
                         "run run_local.py --steps fits,merge for this K first")
    return {s: p for s, p in out.items() if p.exists()}


def outdir(*parts: str) -> Path:
    d = OUT.joinpath(*parts)
    d.mkdir(parents=True, exist_ok=True)
    return d
