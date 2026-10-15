PYTHON ?= python
OLD_REPO ?=

.PHONY: check references smoke mini compare full-fit full-core full-heldout paper-main paper-supp clean-smoke

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

# Full model sweep, fit stage only. This is computationally expensive.
full-fit:
	$(PYTHON) scripts/run_pipeline.py --profile full --output-root . --fit-only

# Full fit + full/single-archetype decoding + held-out spatial top-m. Extremely expensive.
full-core:
	$(PYTHON) scripts/run_pipeline.py --profile full --output-root .

# Run/resume only the canonical held-out sparse analysis after fits already exist.
full-heldout:
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
