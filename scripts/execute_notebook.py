#!/usr/bin/env python3
"""Execute a notebook while forcing its working directory to the repository root.

This keeps the historical notebooks' repo-relative paths working after they have
been organized into subdirectories.
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
args = ap.parse_args()

src = (REPO / args.notebook).resolve() if not Path(args.notebook).is_absolute() else Path(args.notebook)
nb = nbformat.read(src, as_version=4)
ep = ExecutePreprocessor(timeout=args.timeout, kernel_name="python3")
ep.preprocess(nb, {"metadata": {"path": str(REPO)}})
out = Path(args.output).resolve() if args.output else src.with_name(src.stem + "_executed.ipynb")
nbformat.write(nb, out)
print(out)
