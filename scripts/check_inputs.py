#!/usr/bin/env python3
from __future__ import annotations
import json, sys
from pathlib import Path
REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))
from msaa_pipeline.data import describe_inputs
try:
    print(json.dumps(describe_inputs(REPO), indent=2))
except Exception as e:
    print(f"INPUT CHECK FAILED: {e}", file=sys.stderr)
    raise
