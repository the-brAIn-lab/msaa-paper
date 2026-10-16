#!/usr/bin/env python3
from __future__ import annotations

from pathlib import Path
import io
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from PIL import Image, ImageChops

try:
    from nilearn import plotting as niplot
except Exception as exc:
    raise SystemExit(
        "Figure 5 rendering requires nilearn. Install/activate the paper environment first."
    ) from exc

SELECTED_ARCHETYPES = [4, 5, 6, 12, 8]
CONDITIONS = ["intact", "word", "rest"]
NETWORK_ORDER = [
    "Visual",
    "Somatomotor",
    "Dorsal attention",
    "Ventral attention",
    "Limbic",
    "Frontoparietal",
    "Default mode",
]
NETWORK_COLORS = {
    "Visual": "#D7DF23",
    "Somatomotor": "#39B54A",
    "Dorsal attention": "#00A79D",
    "Ventral attention": "#27AAE1",
    "Limbic": "#1C75BC",
    "Frontoparietal": "#92278F",
    "Default mode": "#EE2A7B",
}
DISPLAY_THR_FRAC = 0.30
DPI = 300

SOURCE_DIR = Path("reference_outputs/legacy_full/paper_figure_outputs")
OUTPUT_DIR = Path("outputs/paper_figures")
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

SPATIAL_CSV = SOURCE_DIR / "Figure5_spatial_maps_allnodes.csv"
TEMPORAL_CSV = SOURCE_DIR / "Figure5_temporal_profiles.csv"
NETWORK_CSV = SOURCE_DIR / "Figure5_Yeo7_allnode_mass.csv"
STABILITY_CSV = SOURCE_DIR / "Figure5_crossfit_selection_frequency.csv"

FIGURE_PNG = OUTPUT_DIR / "Figure5_selected_K46_archetypes.png"
FIGURE_PDF = OUTPUT_DIR / "Figure5_selected_K46_archetypes.pdf"

for p in [SPATIAL_CSV, TEMPORAL_CSV, NETWORK_CSV, STABILITY_CSV]:
    if not p.is_file():
        raise FileNotFoundError(f"Missing Figure 5 source table: {p}")

spatial = pd.read_csv(SPATIAL_CSV)
temporal = pd.read_csv(TEMPORAL_CSV)
network_mass = pd.read_csv(NETWORK_CSV)
stability = pd.read_csv(STABILITY_CSV)

# Normalize likely column names while remaining strict about the intended content.
required_spatial = {"archetype", "x", "y", "z", "network", "mean_spatial_weight"}
required_temporal = {"archetype", "timepoint", "mean_temporal_profile", "sem_temporal_profile"}
required_network = {"archetype", "network", "mass_fraction"}
required_stability = {"archetype", "condition", "selection_frequency"}

checks = [
    (spatial, required_spatial, SPATIAL_CSV),
    (temporal, required_temporal, TEMPORAL_CSV),
    (network_mass, required_network, NETWORK_CSV),
    (stability, required_stability, STABILITY_CSV),
]
for df, req, path in checks:
    missing = req - set(df.columns)
    if missing:
        raise ValueError(f"{path} is missing columns: {sorted(missing)}")


def scale_node_sizes(values, size_min=18, size_max=105):
    values = np.asarray(values, dtype=float)
    if len(values) == 0:
        return np.array([])
    vmin, vmax = values.min(), values.max()
    if np.isclose(vmax, vmin):
        return np.full(len(values), (size_min + size_max) / 2)
    scaled = (values - vmin) / (vmax - vmin)
    return size_min + scaled * (size_max - size_min)


def render_glass_brain_view(coords, networks, sizes, view):
    tmp_fig = plt.figure(figsize=(4.0, 3.4))
    display = niplot.plot_glass_brain(
        None,
        display_mode=view,
        figure=tmp_fig,
        black_bg=False,
        colorbar=False,
        plot_abs=False,
        annotate=True,
    )
    for network in NETWORK_ORDER:
        mask = networks == network
        if np.any(mask):
            display.add_markers(
                marker_coords=coords[mask],
                marker_color=NETWORK_COLORS[network],
                marker_size=sizes[mask].tolist(),
            )

    buf = io.BytesIO()
    tmp_fig.savefig(
        buf,
        format="png",
        dpi=220,
        facecolor="white",
        bbox_inches="tight",
        pad_inches=0.02,
    )
    plt.close(tmp_fig)
    buf.seek(0)
    img = Image.open(buf).convert("RGB")

    white = Image.new("RGB", img.size, "white")
    diff = ImageChops.difference(img, white).convert("L")
    diff = diff.point(lambda p: 255 if p > 8 else 0)
    bbox = diff.getbbox()
    if bbox is not None:
        img = img.crop(bbox)
    arr = np.asarray(img)
    buf.close()
    return arr


composition_table = (
    network_mass
    .pivot_table(index="archetype", columns="network", values="mass_fraction", aggfunc="first")
    .reindex(index=SELECTED_ARCHETYPES, columns=NETWORK_ORDER)
)

stability_heatmap = (
    stability
    .pivot_table(index="archetype", columns="condition", values="selection_frequency", aggfunc="mean")
    .reindex(index=SELECTED_ARCHETYPES, columns=CONDITIONS)
)

fig = plt.figure(figsize=(17, 22))
outer = fig.add_gridspec(2, 1, height_ratios=[6.4, 1.8], hspace=0.22)
top = outer[0].subgridspec(
    len(SELECTED_ARCHETYPES), 4,
    width_ratios=[1.55, 1.15, 0.90, 1.15],
    wspace=0.06,
    hspace=0.34,
)

for row_i, k in enumerate(SELECTED_ARCHETYPES):
    ax_t = fig.add_subplot(top[row_i, 0])
    cur_t = temporal[temporal["archetype"].astype(int) == int(k)].sort_values("timepoint")
    t = cur_t["timepoint"].to_numpy()
    mean_tc = cur_t["mean_temporal_profile"].to_numpy(dtype=float)
    sem_tc = cur_t["sem_temporal_profile"].to_numpy(dtype=float)

    ax_t.plot(t, mean_tc, linewidth=1.8)
    ax_t.fill_between(t, mean_tc - sem_tc, mean_tc + sem_tc, alpha=0.20)
    ax_t.axhline(0, linewidth=0.7, linestyle=":")
    ax_t.set_xlim(t.min(), t.max())
    ax_t.set_ylabel(f"A{k}", rotation=0, labelpad=28, fontsize=14, fontweight="bold", va="center")
    ax_t.tick_params(axis="both", labelsize=12)
    if row_i == len(SELECTED_ARCHETYPES) - 1:
        ax_t.set_xlabel("Timepoint", fontsize=14)
    else:
        ax_t.set_xticklabels([])
    if row_i == 0:
        ax_t.set_title("Temporal motif", fontweight="bold", fontsize=15, pad=10)

    cur_s = spatial[spatial["archetype"].astype(int) == int(k)].copy()
    weights = cur_s["mean_spatial_weight"].to_numpy(dtype=float)
    max_weight = np.nanmax(weights)
    keep = weights >= DISPLAY_THR_FRAC * max_weight
    shown = cur_s.loc[keep]
    coords = shown[["x", "y", "z"]].to_numpy(dtype=float)
    nets = shown["network"].astype(str).to_numpy()
    sizes = scale_node_sizes(shown["mean_spatial_weight"].to_numpy(dtype=float))

    brain_axes = []
    for col_i, view in enumerate(["l", "z", "r"], start=1):
        ax_b = fig.add_subplot(top[row_i, col_i])
        brain_axes.append(ax_b)
        brain_img = render_glass_brain_view(coords, nets, sizes, view)
        ax_b.imshow(brain_img, interpolation="antialiased")
        ax_b.set_axis_off()
        ax_b.set_anchor("C")
    if row_i == 0:
        brain_axes[1].set_title("Spatial deployment", fontweight="bold", fontsize=15, pad=10)

bottom = outer[1].subgridspec(1, 2, width_ratios=[1.70, 1.00], wspace=0.18)

ax_comp = fig.add_subplot(bottom[0, 0])
ypos = np.arange(len(SELECTED_ARCHETYPES))
left = np.zeros(len(SELECTED_ARCHETYPES))
for network in NETWORK_ORDER:
    vals = composition_table[network].to_numpy(dtype=float)
    ax_comp.barh(
        ypos, vals, left=left,
        color=NETWORK_COLORS[network],
        edgecolor="none", height=0.72,
        label=network,
    )
    left += vals
ax_comp.set_yticks(ypos)
ax_comp.set_yticklabels([f"A{k}" for k in SELECTED_ARCHETYPES], fontsize=15)
ax_comp.invert_yaxis()
ax_comp.set_xlim(0, 1)
ax_comp.set_xlabel("Fraction of |spatial weight|", fontsize=15, labelpad=8)
ax_comp.tick_params(axis="x", labelsize=15)
ax_comp.set_title("Yeo-7 composition across selected archetypes", fontweight="bold", fontsize=17)

ax_heat = fig.add_subplot(bottom[0, 1])
heat = stability_heatmap.to_numpy(dtype=float)
im = ax_heat.imshow(heat, aspect="auto", vmin=0, vmax=1)
ax_heat.set_xticks(np.arange(len(CONDITIONS)))
ax_heat.set_xticklabels(["Intact", "Word", "Rest"], fontsize=15)
ax_heat.set_yticks(np.arange(len(SELECTED_ARCHETYPES)))
ax_heat.set_yticklabels([f"A{k}" for k in SELECTED_ARCHETYPES], fontsize=15)
ax_heat.set_title("Cross-fitted top-5 selection stability", fontweight="bold", fontsize=17)
for r in range(heat.shape[0]):
    for c in range(heat.shape[1]):
        val = heat[r, c]
        txt_color = "white" if val >= 0.55 else "black"
        ax_heat.text(c, r, f"{val:.2f}", ha="center", va="center", fontsize=13, color=txt_color)

cbar = fig.colorbar(im, ax=ax_heat, fraction=0.046, pad=0.04)
cbar.set_label("Selection frequency", fontsize=15)
cbar.ax.tick_params(labelsize=14)

legend_handles = [plt.Rectangle((0, 0), 1, 1, color=NETWORK_COLORS[n]) for n in NETWORK_ORDER]
fig.legend(
    legend_handles,
    NETWORK_ORDER,
    loc="lower center",
    ncol=4,
    frameon=False,
    fontsize=14,
    bbox_to_anchor=(0.5, 0.018),
)
fig.suptitle("Held-out-selected spatial archetypes at K=46", fontsize=19, fontweight="bold", y=0.995)
fig.subplots_adjust(bottom=0.115, top=0.965)
fig.savefig(FIGURE_PNG, dpi=DPI, bbox_inches="tight")
fig.savefig(FIGURE_PDF, bbox_inches="tight")
plt.close(fig)

print("Loaded Figure 5 source tables from:", SOURCE_DIR)
print("Saved:", FIGURE_PNG)
print("Saved:", FIGURE_PDF)
