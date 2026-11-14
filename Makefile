# Run with the project environment, e.g. `uv run make figures` (or PYTHON=.venv/bin/python make ...).
PYTHON ?= python
ROOT ?= .
OUT ?= $(ROOT)/figure_bundle

.PHONY: check atlas smoke results figures

# the two HTFA inputs are in data/input/
check:
	$(PYTHON) scripts/check_inputs.py

# Schaefer-1000 / Yeo-7 atlas into data/reference/ (Figure 5 and the node table)
atlas:
	$(PYTHON) scripts/fetch_reference_atlas.py

# small end-to-end run (4 participants, K = 5) into outputs/smoke/; minutes
smoke:
	$(PYTHON) scripts/run_pipeline.py --profile smoke

# every fit, decode and held-out run behind the paper (very large; see docs/REPRODUCE.md)
results:
	$(PYTHON) analysis/08_multiseed_stability/run_local.py

# draw every figure from finished results in ROOT and assemble the bundle in OUT
figures:
	$(PYTHON) scripts/build_figures.py --root $(ROOT) --out $(OUT)
