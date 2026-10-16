PYTHON ?= python
OLD_REPO ?=

.PHONY: check references smoke mini compare import-legacy validate-legacy validate-legacy-hash activate-legacy deactivate-legacy handoff-check full-fit full-core full-heldout paper-main paper-supp clean-smoke

check:
	$(PYTHON) scripts/check_inputs.py

references:
	$(PYTHON) scripts/fetch_reference_atlas.py

smoke:
	$(PYTHON) scripts/run_pipeline.py --profile smoke

mini:
	$(PYTHON) scripts/run_pipeline.py --profile mini-paper

compare:
	@if [ -z "$(OLD_REPO)" ]; then echo "Usage: make compare OLD_REPO=/path/to/original/MS_AA"; exit 2; fi
	$(PYTHON) scripts/compare_to_old_repo.py "$(OLD_REPO)" --profile mini-paper

# Copy the expensive original derived outputs into an ignored, portable store
# inside this project and write a tracked SHA-256 manifest. Resumable.
import-legacy:
	@if [ -z "$(OLD_REPO)" ]; then echo "Usage: make import-legacy OLD_REPO=/path/to/original/MS_AA"; exit 2; fi
	$(PYTHON) scripts/import_legacy_outputs.py "$(OLD_REPO)"

# Quick integrity check: existence + byte size.
validate-legacy:
	$(PYTHON) scripts/validate_legacy_outputs.py --mode quick

# Full integrity check: recompute SHA-256 for every imported file.
validate-legacy-hash:
	$(PYTHON) scripts/validate_legacy_outputs.py --mode hash

# Make the portable store visible at the historical top-level output paths
# expected by the original downstream notebooks. Links are relative and portable.
activate-legacy: validate-legacy
	$(PYTHON) scripts/activate_legacy_outputs.py

deactivate-legacy:
	$(PYTHON) scripts/activate_legacy_outputs.py --deactivate

handoff-check:
	$(PYTHON) scripts/check_handoff.py

# Full model sweep, fit stage only. This is computationally expensive.
full-fit:
	$(PYTHON) scripts/guard_no_legacy_links.py
	$(PYTHON) scripts/run_pipeline.py --profile full --output-root . --fit-only

# Full fit + full/single-archetype decoding + held-out spatial top-m. Extremely expensive.
full-core:
	$(PYTHON) scripts/guard_no_legacy_links.py
	$(PYTHON) scripts/run_pipeline.py --profile full --output-root .

# Run/resume only the canonical held-out sparse analysis after fits already exist.
full-heldout:
	$(PYTHON) scripts/guard_no_legacy_links.py
	$(PYTHON) scripts/run_pipeline.py --profile full --output-root . --heldout-only

# These assume the required upstream canonical result tables already exist.
paper-main:
	$(PYTHON) scripts/execute_notebook.py paper/main_figures/paper_figure02_empirical_PCA_MSAA_benchmark_main.ipynb
	$(PYTHON) scripts/execute_notebook.py paper/main_figures/paper_figure03_PCA_MSAA_simulation_benchmark.ipynb
	$(PYTHON) scripts/execute_notebook.py paper/main_figures/paper_figure04_heldout_sparse_MSAA_condensed_fixed.ipynb
	$(PYTHON) scripts/execute_notebook.py paper/main_figures/paper_figure05_selected_K46_archetypes.ipynb
	$(PYTHON) scripts/execute_notebook.py paper/main_figures/paper_figure06_support_and_robustness.ipynb

clean-smoke:
	rm -rf outputs/smoke/* outputs/mini_paper/*
	touch outputs/smoke/.gitkeep outputs/mini_paper/.gitkeep
