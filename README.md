# MS-AA paper

Code for the multi-subject archetypal analysis (MS-AA) paper: from the two HTFA inputs to every
figure and source table in the paper and supplement.

## Setup

```bash
uv sync
# put pieman_data.mat and pieman_posterior_K700.mat in data/input/ (see data/README.md)
uv run make check      # the inputs can be read
uv run make atlas      # Schaefer-1000 / Yeo-7 atlas into data/reference/
uv run make smoke      # small end-to-end run (4 participants, K = 5), a few minutes
```

Without uv (for example on a cluster), install the same pinned versions into a Python 3.12
environment with `pip install -r requirements.txt` (exported from `uv.lock` with
`uv export --frozen --no-hashes --no-emit-project -o requirements.txt`; regenerate it after changing
dependencies) and use that environment's `python` in place of `uv run`.

## Reproduce

```bash
uv run make results                      # every fit, decode and held-out run (very large)
uv run make figures OUT=figure_bundle    # every figure + source table, from those results
```

`docs/REPRODUCE.md` lists which step makes each result, which figure uses it, and how to run a
subset. `make figures ROOT=/path/to/results` draws from results kept outside the repo.

## Layout

```text
src/msaa_pipeline/    fitting, decoding, participant-held-out top-m (used by every analysis)
helpers/              the MS-AA implementation (MultiSubject_AA.py, FurthestSum)
hypertools/           the normalisation the inputs go through (kept for exact parity)
scripts/              run_parallel.py (fit/decode/held-out driver), build_figures.py, input checks,
                      atlas fetch, the frozen FurthestSum picks
analysis/             the analyses behind the figures, in pipeline order (see analysis/README.md)
paper/                figure notebooks/scripts that only draw, captions, bundle README
configs/              frozen FurthestSum picks, held-out outer-split seeds, node -> Yeo-7 table, profiles
data/                 input/ (the two HTFA .mat files) and reference/ (atlas); not committed
```

The HTFA inputs themselves (raw fMRI to HTFA) are made outside this repository.

## Source tables

`paper/source_tables/` holds the exact tables behind every figure in the paper and supplement, with
a README mapping each file to its figure.

## License

MIT; see `LICENSE`. The vendored `hypertools/` folder is from hypertools (Contextual Dynamics
Laboratory) and keeps its own MIT license in `hypertools/LICENSE`; `hypertools/srm.py` keeps its
Apache-2.0 header (Intel Corporation). `helpers/MultiSubject_AA.py` adapts the multisubject AA
implementation of Hinrich et al. (2016), doi:10.1109/JSTSP.2016.2595103.
