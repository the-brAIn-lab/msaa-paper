from __future__ import annotations
from pathlib import Path
import json
import pandas as pd
import matplotlib.pyplot as plt


def build_qc(output_root: Path, profile_name: str):
    rows = []
    for p in sorted(output_root.glob("msaa_condrank_decoding_outputs_*/full_summary.csv")):
        df = pd.read_csv(p)
        rows.append(df)
    if not rows:
        raise FileNotFoundError("No decoding full_summary.csv files found")
    all_df = pd.concat(rows, ignore_index=True)
    qc_dir = output_root / "qc"
    qc_dir.mkdir(parents=True, exist_ok=True)
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

    report = {
        "profile": profile_name,
        "n_summary_rows": int(len(all_df)),
        "outputs": [str(p.relative_to(output_root)) for p in sorted(output_root.rglob("*.csv"))],
    }
    (qc_dir / "pipeline_qc.json").write_text(json.dumps(report, indent=2))
    return report
