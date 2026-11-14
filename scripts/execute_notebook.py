#!/usr/bin/env python3
"""Execute a notebook with its working directory set to the repository root (or --cwd).

The notebooks use paths relative to the run root (results and code side by side), so
they run unchanged from any results root that has src/, helpers/ and hypertools/.
"""
from __future__ import annotations
import argparse
from pathlib import Path
import nbformat
from nbconvert.preprocessors import ExecutePreprocessor

REPO = Path(__file__).resolve().parents[1]
ap = argparse.ArgumentParser()
ap.add_argument("notebook")
ap.add_argument("--output", default=None)
ap.add_argument("--timeout", type=int, default=-1)
ap.add_argument("--cwd", default=str(REPO), help="working directory (run root) for the notebook")
args = ap.parse_args()

src = (REPO / args.notebook).resolve() if not Path(args.notebook).is_absolute() else Path(args.notebook)
nb = nbformat.read(src, as_version=4)
ep = ExecutePreprocessor(timeout=args.timeout, kernel_name="python3")
ep.preprocess(nb, {"metadata": {"path": str(Path(args.cwd).resolve())}})
out = Path(args.output).resolve() if args.output else src.with_name(src.stem + "_executed.ipynb")
nbformat.write(nb, out)
print(out)
