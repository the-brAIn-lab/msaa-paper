# Run with the project environment, e.g. `uv run make figures` (or PYTHON=.venv/bin/python make ...).
PYTHON ?= python
ROOT ?= .
OUT ?= $(ROOT)/figure_bundle

.PHONY: data check test atlas smoke results figures

# convert the Zenodo pieman2_htfa.pkl into the two HTFA inputs in data/input/:
#   make data PKL=/path/to/pieman2_htfa.pkl
data:
	@test -n "$(PKL)" || { echo "usage: make data PKL=/path/to/pieman2_htfa.pkl"; exit 2; }
	$(PYTHON) scripts/pkl_to_mat.py "$(PKL)"

# the two HTFA inputs are in data/input/, and the low-rank decoding path matches the dense one
check: test
	$(PYTHON) scripts/check_inputs.py

# fast-path parity checks on synthetic data (seconds; no inputs needed)
test:
	$(PYTHON) tests/test_lowrank_parity.py

# Schaefer-1000 / Yeo-7 atlas into data/reference/ (Figure 5 and the node table).
# Downloads from the internet: on a cluster, run it on a login node before offline batch jobs.
atlas:
	$(PYTHON) scripts/fetch_reference_atlas.py

# small end-to-end run (4 participants, K = 5, 20 iterations) into outputs/smoke/; seconds
smoke:
	$(PYTHON) scripts/run_pipeline.py --profile smoke

# every fit, decode and held-out run behind the paper (very large; see docs/REPRODUCE.md)
results:
	$(PYTHON) analysis/08_multiseed_stability/run_local.py

# draw every figure from finished results in ROOT and assemble the bundle in OUT
figures:
	$(PYTHON) scripts/build_figures.py --root $(ROOT) --out $(OUT)
