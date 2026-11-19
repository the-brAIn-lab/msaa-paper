#!/usr/bin/env python3
"""Draw every paper figure from finished results and assemble the figure bundle.

    uv run python scripts/build_figures.py [--root RUN_ROOT] [--out BUNDLE_DIR] [--only fig5,s10,...]

RUN_ROOT holds the results the figures are drawn from (default: the repo root, where the
pipeline writes them; see docs/REPRODUCE.md for the layout and which step makes each one).
If RUN_ROOT is not the repo, the code folders (src, helpers, hypertools, configs) are linked
into it, because the notebooks use paths relative to the run root.

Nothing here refits or re-decodes: each step only reads result tables (and, for Figure 5,
the K = 50 fit) and draws. Intermediate figure files go to RUN_ROOT/outputs/; the bundle
(main-png/, main-pdf/, supplement-png/, supplement-pdf/, tables/, README.txt)
goes to BUNDLE_DIR.
"""
from __future__ import annotations

import argparse
import os
import shutil
import subprocess
import sys
import time
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
MS = REPO / "analysis" / "08_multiseed_stability"
K_NEIGHBOURHOOD = "45,46,47,48,49,50,51,52,53,54,55"

# step name -> (command, figures it draws). Commands run with cwd = RUN_ROOT and MSAA_ROOT = RUN_ROOT;
# "{root}" is replaced by RUN_ROOT and "nb:" marks a notebook.
STEPS = {
    "fig2": [["nb:paper/main_figures/paper_figure02_empirical_PCA_MSAA_benchmark_main.ipynb"],
             ["paper/main_figures/figure02_plot.py", "--sharey"]],
    "fig3": [["nb:paper/main_figures/paper_figure03_PCA_MSAA_simulation_benchmark.ipynb"]],
    "fig4": [[f"{MS}/plots/fig4_heldout_panels.py"]],
    "fig4c": [[f"{MS}/plots/fig4c_matched_ratio_heatmaps.py"]],
    "fig5": [["nb:paper/main_figures/paper_figure05_selected_K50_archetypes.ipynb"]],
    "fig6": [["nb:paper/main_figures/paper_figure06_support_and_robustness.ipynb"]],
    "fig8": [["paper/main_figures/figure08_support_robustness.py"]],   # Fig 6 panels A-D plus E (retention) and F (edge trim)
    "s02": [[f"{MS}/plots/k14_supplement.py"]],
    "s03": [[f"{MS}/matched_fraction/k14_local_contrast.py",
             "{root}/outputs/multiseed_stability/matched_fraction/aggregate",
             "{root}/outputs/multiseed_stability/matched_fraction/k14_local_contrast"]],
    "s04": [[f"{MS}/k50_stability/m0m4_summary.py"]],
    "s05": [[f"{MS}/plots/topm_stability_neighbourhood.py", K_NEIGHBOURHOOD, "--plot-only"],
            [f"{MS}/k50_stability/make_supp_tables.py"]],
    "s07": [["paper/supplement/supp_atlas_condition_plot.py", "--root", "{root}"]],
    "s10": [["nb:paper/supplement/supp_05_edge_trim_diagnostics.ipynb"]],
    "s11": [["paper/supplement/supp_four_condition_plot.py",
             "--source", "{root}/four_condition_ratio_decoding_outputs"]],
    "s12": [[f"{MS}/plots/spatial_topm_landscape.py"]],
    "s13": [["paper/supplement/supp_13_simulation_sensitivity.py"]],
}
# S07 needs the Figure 5 notebook's audit tables; S04/S05 tables need nothing else drawn first.
ORDER = ["fig2", "fig3", "fig4", "fig4c", "fig5", "fig6", "s02", "s03", "s04", "s05", "s07", "s10", "s11",
         "s12", "s13", "fig8"]   # fig8 reads the S10 edge-trim summary

O, M = "outputs", "outputs/multiseed_stability"
# bundle file (without extension) <- figure written by the step (without extension), relative to RUN_ROOT
FIGURES = {
    "main/Figure2_PCA_vs_MSAA_benchmark": f"{O}/paper_figures/Figure2_empirical_PCA_MSAA_benchmark_sharey",
    "main/Figure3_simulation_motif_vs_span": f"{O}/paper_figures/Figure3_PCA_MSAA_simulation_benchmark",
    "main/Figure4A_spatial_fixed_topm": f"{M}/fig4_heldout_panels/Fig4A_spatial_heldout_topm",
    "main/Figure4B_temporal_fixed_topm": f"{M}/fig4_heldout_panels/Fig4B_temporal_heldout_topm",
    "main/Figure4C_spatial_minus_temporal_matched_ratio": f"{M}/fig4c_matched_ratio_heatmaps/Fig4C_matched_ratio_heatmaps_intact",
    "main/Figure4D_temporal_matched_retained_fraction": f"{M}/fig4_heldout_panels/Fig4D_temporal_matched_fraction",
    "main/Figure5_selected_spatial_motifs_M0-M4": "paper_figure_outputs/Figure5_selected_K50_archetypes",
    "main/Figure6_motif_support_and_robustness": f"{O}/paper_figures/Figure6_selected_motif_support_robustness",
    "main/Figure8_motif_support_and_robustness_A-F": f"{O}/paper_figures/Figure8_motif_support_and_robustness",
    "supplement/S01_spatial_minus_temporal_heatmaps_and_seed_consistency": f"{M}/fig4c_matched_ratio_heatmaps/FigS_matched_ratio_heatmaps_all_conditions",
    "supplement/S02_temporal_K14_accumulation_anchor_masking_higherK_splitting": f"{M}/k14_supplement/FigS_temporal_K14_accumulation_anchor_fragmentation",
    "supplement/S03_temporal_K12_14_16_seed_contrasts": f"{M}/matched_fraction/k14_local_contrast/per_seed_K12_14_16",
    "supplement/S04_M0-M4_cross_seed_stability": f"{M}/m0m4_summary/M0M4_stability_summary",
    "supplement/S05a_K45-55_selection_stability_top5": f"{M}/topm_stability_neighbourhood/topm_stability_neighbourhood_top5",
    "supplement/S05b_K45-55_selection_stability_top10": f"{M}/topm_stability_neighbourhood/topm_stability_neighbourhood_top10",
    "supplement/S06_M0-M4_selected_and_recovered_top5_top10": f"{M}/topm_stability_neighbourhood/M0M4_top5_vs_top10_K50",
    "supplement/S07_atlas_mapping_audit": f"{O}/paper_supplement/Supplement_M0M4_atlas_mapping_audit",
    "supplement/S08_condition_contributions": f"{O}/paper_supplement/Supplement_M0M4_condition_contributions",
    "supplement/S09_full_model_accuracy_vs_K": f"{M}/fig4_heldout_panels/FigS_full_model_accuracy_vs_K",
    "supplement/S10_edge_trim_diagnostics": f"{O}/paper_supplement/Supplement_edge_trim_diagnostics",
    "supplement/S11_four_condition_paragraph_sensitivity": f"{O}/paper_supplement/Supplement_four_condition_sensitivity",
    "supplement/S12_spatial_heldout_topm_landscape_K_by_m": f"{M}/spatial_topm_landscape/FigS_spatial_topm_landscape",
    "supplement/S13_simulation_sensitivity": f"{O}/paper_supplement/S13_simulation_sensitivity",
}
# bundle tables/ file <- source file, relative to RUN_ROOT
TABLES = {
    "Figure2_source.csv": f"{O}/paper_figures/Figure2_source_data.csv",
    "Figure4A_source.csv": f"{M}/fig4_heldout_panels/Fig4A_spatial_heldout_topm_source.csv",
    "Figure4B_source.csv": f"{M}/fig4_heldout_panels/Fig4B_temporal_heldout_topm_source.csv",
    "Figure4D_source.csv": f"{M}/fig4_heldout_panels/Fig4D_temporal_matched_fraction_source.csv",
    "spatial_vs_temporal_matched_ratio_per_seed.csv": f"{M}/fig4c_matched_ratio_heatmaps/matched_ratio_per_seed.csv",
    "spatial_vs_temporal_matched_ratio_summary.csv": f"{M}/fig4c_matched_ratio_heatmaps/matched_ratio_summary.csv",
    "Figure5_crossfit_selection_frequency.csv": "paper_figure_outputs/Figure5_crossfit_selection_frequency.csv",
    "Figure5_temporal_profiles.csv": "paper_figure_outputs/Figure5_temporal_profiles.csv",
    "Figure5_Yeo7_allnode_mass.csv": "paper_figure_outputs/Figure5_Yeo7_allnode_mass.csv",
    "Figure6A_cross_seed_map_r_source.csv": f"{O}/paper_figures/Figure6_cross_seed_map_r_source.csv",
    "Figure6B_C_support_source.csv": "revision_05_archetype_interpretation_and_support/selected_archetype_support_by_independent_participant.csv",
    "Figure6D_S08_LOCO_map_r_source.csv": f"{O}/paper_figures/Figure6_loco_map_r_source.csv",
    "Figure8E_K50_retention_source.csv": f"{O}/paper_figures/Figure8E_K50_retention_source.csv",
    "Figure8F_K50_edge_trim_source.csv": f"{O}/paper_figures/Figure8F_K50_edge_trim_source.csv",
    "temporal_K14_anchor_masked_vs_unmasked.csv": f"{O}/temporal_compact_dictionary_K14_local/K14_anchor_masked_vs_unmasked.csv",
    "temporal_K14_local_contrast_delta14.csv": f"{M}/matched_fraction/k14_local_contrast/delta14_contrast.csv",
    "temporal_matched_fraction_seed_level.csv": f"{M}/matched_fraction/k14_local_contrast/seed_level_table.csv",
    "M0-M4_cross_seed_stability_summary.csv": f"{M}/m0m4_summary/M0M4_stability_summary.csv",
    "K45-55_selection_stability.csv": f"{M}/topm_stability_neighbourhood/topm_stability_neighbourhood.csv",
    "M0-M4_selection_by_K.csv": f"{M}/topm_stability_neighbourhood/M0M4_selection_by_K.csv",
    "M0-M4_raw_ID_crosswalk.csv": f"{M}/k50_supp_tables/M0M4_crosswalk.csv",
    "K50_M0-M4_supplementary_table.csv": f"{M}/k50_supp_tables/supp_table_K50_M0M4.csv",
    "S07_atlas_mapping_audit_summary.csv": "paper_figure_outputs/Figure5_HTFA_Yeo_mapping_audit_summary.csv",
    "S09_full_model_accuracy_source.csv": f"{M}/fig4_heldout_panels/FigS_full_model_accuracy_vs_K_source.csv",
    "S10_edge_trim_source.csv": "msaa_edge_trim_sensitivity_K50_top5/edge_trim_accuracy_summary.csv",
    "S11_four_condition_msaa.csv": "four_condition_ratio_decoding_outputs/aa_full_summary_with_retention.csv",
    "S11_four_condition_pca.csv": "four_condition_ratio_decoding_outputs/pca_summary_with_retention.csv",
    "S11_four_condition_raw_htfa.csv": "four_condition_ratio_decoding_outputs/raw_htfa_summary.csv",
    "S11_new_vs_legacy_comparison.csv": "four_condition_comparison/comparison.csv",
    "S12_spatial_topm_landscape.csv": f"{M}/spatial_topm_landscape/spatial_topm_landscape.csv",
    "S12_landscape_run_validation.txt": "heldout-highm-landscape-results/summary.txt",
    "S13_simulation_sensitivity_source.csv": f"{O}/paper_supplement/S13_simulation_sensitivity_source.csv",
    "S13_simulation_sensitivity_runs.csv": "simulation_pca_msaa_extended_outputs/s13_sensitivity_results.csv",
}
CODE_LINKS = ["src", "helpers", "hypertools", "configs"]


def log(msg):
    print(time.strftime("%H:%M:%S"), msg, flush=True)


def link_code(root: Path):
    """The notebooks import from <run root>/src and read <run root>/hypertools; link the repo's copies in."""
    if root == REPO:
        return
    for name in CODE_LINKS:
        p = root / name
        if p.is_symlink() or not p.exists():
            if p.is_symlink():
                p.unlink()
            p.symlink_to(REPO / name)
        elif p.resolve() != (REPO / name).resolve():
            raise SystemExit(f"{p} exists and is not a link to the repo's {name}/: move it out of the run root")


def run_step(name: str, root: Path, env: dict):
    for cmd in STEPS[name]:
        cmd = [c.replace("{root}", str(root)) for c in cmd]
        if cmd[0].startswith("nb:"):
            nb = REPO / cmd[0][3:]
            out = root / "outputs" / "executed_notebooks" / nb.name
            out.parent.mkdir(parents=True, exist_ok=True)
            full = [sys.executable, str(REPO / "scripts" / "execute_notebook.py"), str(nb), "--cwd", str(root),
                    "--output", str(out)]
        else:
            script = Path(cmd[0]) if Path(cmd[0]).is_absolute() else REPO / cmd[0]
            full = [sys.executable, str(script), *cmd[1:]]
        log(f"[{name}] " + " ".join(Path(c).name if c.startswith(str(REPO)) else c for c in full[1:]))
        r = subprocess.run(full, cwd=root, env=env, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
        if r.returncode != 0:
            print(r.stdout[-4000:])
            raise SystemExit(f"step {name} failed (exit {r.returncode})")


def assemble(root: Path, out: Path):
    missing = []
    for dest, src in FIGURES.items():
        part, stem = dest.split("/")
        for ext in ("png", "pdf"):
            s = root / f"{src}.{ext}"
            if not s.exists():
                missing.append(str(s)); continue
            d = out / f"{part}-{ext}"; d.mkdir(parents=True, exist_ok=True)
            shutil.copy2(s, d / f"{stem}.{ext}")
    (out / "tables").mkdir(parents=True, exist_ok=True)
    for dest, src in TABLES.items():
        s = root / src
        if not s.exists():
            missing.append(str(s)); continue
        shutil.copy2(s, out / "tables" / dest)
    shutil.copy2(REPO / "paper" / "FIGURE_BUNDLE_README.txt", out / "README.txt")
    return missing


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--root", default=str(REPO), help="run root holding the results (default: the repo)")
    ap.add_argument("--out", default=None, help="bundle folder (default: RUN_ROOT/figure_bundle)")
    ap.add_argument("--only", default=None, help=f"comma list of steps: {','.join(ORDER)}")
    ap.add_argument("--no-draw", action="store_true", help="only assemble the bundle from figures already drawn")
    a = ap.parse_args()
    root = Path(a.root).expanduser().resolve()
    out = Path(a.out).expanduser().resolve() if a.out else root / "figure_bundle"
    steps = ORDER if not a.only else a.only.split(",")
    bad = [s for s in steps if s not in STEPS]
    if bad:
        raise SystemExit(f"unknown steps {bad}; choose from {ORDER}")
    env = dict(os.environ, MSAA_ROOT=str(root), MSAA_TIE_TOL=os.environ.get("MSAA_TIE_TOL", "1e-08"),
               PYTHONPATH=os.pathsep.join([str(REPO / "src"), os.environ.get("PYTHONPATH", "")]).rstrip(os.pathsep),
               MPLBACKEND="Agg")
    # Figure 6 reads the K = 50 cross-seed alignment from the run root
    env.setdefault("MSAA_SEED_ALIGNMENT",
                   str(root / "outputs/multiseed_stability/k50_cross_seed_stability/K50_paper_archetypes_alignment.csv"))
    link_code(root)
    if not a.no_draw:
        for s in steps:
            run_step(s, root, env)
    missing = assemble(root, out)
    if missing:
        print("not in the bundle (inputs missing or their step not run):")
        for m in missing:
            print("  ", m)
    log(f"bundle: {out}")
    sys.exit(1 if missing else 0)


if __name__ == "__main__":
    main()
