from __future__ import annotations
from pathlib import Path
import json
import pandas as pd
import matplotlib.pyplot as plt


def build_qc(output_root: Path, profile_name: str):
    rows = []
    for p in sorted(output_root.glob("msaa_condrank_decoding_outputs_*/full_summary.csv")):
        rows.append(pd.read_csv(p))

    qc_dir = output_root / "qc"
    qc_dir.mkdir(parents=True, exist_ok=True)

    full_rows = 0
    if rows:
        all_df = pd.concat(rows, ignore_index=True)
        full_rows = len(all_df)
        all_df.to_csv(qc_dir / "full_decoding_summary.csv", index=False)

        fig, ax = plt.subplots(figsize=(8, 4.5))
        labels = []
        values = []
        for _, r in all_df.sort_values(["analysis_type", "fit_scope", "K", "condition"]).iterrows():
            labels.append(f"{r.analysis_type[0].upper()}-{r.fit_scope[0].upper()}\nK{int(r.K)}\n{r.condition}")
            values.append(r["mean"])
        ax.bar(range(len(values)), values)
        ax.set_xticks(range(len(values)))
        ax.set_xticklabels(labels, rotation=90, fontsize=7)
        ax.set_ylabel("Timepoint decoding accuracy")
        ax.set_title(f"{profile_name} end-to-end QC")
        fig.tight_layout()
        fig.savefig(qc_dir / "pipeline_qc.png", dpi=180)
        plt.close(fig)

    heldout_rows = 0
    heldout_files = []
    for hp in sorted(output_root.glob("msaa_heldout_topm_outputs_*_across/heldout_topm_summary.csv")):
        heldout_files.append(str(hp.relative_to(output_root)))
        heldout_rows += len(pd.read_csv(hp))

    if not rows and heldout_rows == 0:
        raise FileNotFoundError("No full-decoding or held-out top-m summary files found")

    report = {
        "profile": profile_name,
        "n_full_decoding_summary_rows": int(full_rows),
        "n_heldout_topm_summary_rows": int(heldout_rows),
        "heldout_topm_summaries": heldout_files,
        "outputs": [str(p.relative_to(output_root)) for p in sorted(output_root.rglob("*.csv"))],
    }
    (qc_dir / "pipeline_qc.json").write_text(json.dumps(report, indent=2))
    return report
