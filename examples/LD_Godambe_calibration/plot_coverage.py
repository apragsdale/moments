"""
Plots of 95% CI coverage and of standard-error calibration for FIM, GIM and bootstrap
refits, from the results in results/. Writes figures/coverage.png and
figures/sd_ratio.png.

Usage: python plot_coverage.py
"""

import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

import godambe_calibration as gc

# Reference palette (dataviz skill, light mode): categorical slots 1-4 in fixed order,
# with a distinct marker per series so identity never rests on colour alone.
SERIES = {
    "FIM": ("#2a78d6", "o"),
    "GIM": ("#eb6834", "s"),
    "Bootstrap SE": ("#1baf7a", "D"),
    "Bootstrap percentile": ("#eda100", "^"),
}
SURFACE = "#fcfcfb"
INK = "#0b0b0b"
INK2 = "#52514e"
MUTED = "#898781"
GRID = "#e1e0d9"
AXIS = "#c3c2b7"

PANELS = [("ld", "weak"), ("ld", "mig"), ("sfs", "weak"), ("sfs", "mig")]
TITLES = {
    ("ld", "weak"): "LD, weak migration (M = 0.2)",
    ("ld", "mig"): "LD, migration (M = 2)",
    ("sfs", "weak"): "SFS, weak migration (M = 0.2)",
    ("sfs", "mig"): "SFS, migration (M = 2)",
}
LABELS = ["nu0", "nu1", "T", "M", "Ne"]


def intervals(results):
    """Per-dataset 95% intervals and standard errors for each method."""
    out = {name: ([], [], []) for name in SERIES}
    for r in results:
        p = r["p_full"]
        boot_sd = r["p_boot"].std(axis=0, ddof=1)
        for name, sd in [("FIM", r["sd"]["FIM"]), ("GIM", r["sd"]["GIM"])]:
            out[name][0].append(p - 1.96 * sd)
            out[name][1].append(p + 1.96 * sd)
            out[name][2].append(sd)
        out["Bootstrap SE"][0].append(p - 1.96 * boot_sd)
        out["Bootstrap SE"][1].append(p + 1.96 * boot_sd)
        out["Bootstrap SE"][2].append(boot_sd)
        lo, hi = np.percentile(r["p_boot"], [2.5, 97.5], axis=0)
        out["Bootstrap percentile"][0].append(lo)
        out["Bootstrap percentile"][1].append(hi)
        out["Bootstrap percentile"][2].append(boot_sd)
    return {k: tuple(np.array(x) for x in v) for k, v in out.items()}


def style(ax):
    ax.set_facecolor(SURFACE)
    for side in ["top", "right", "left"]:
        ax.spines[side].set_visible(False)
    ax.spines["bottom"].set_color(AXIS)
    ax.tick_params(colors=INK2, labelsize=9, length=0)
    ax.grid(axis="y", color=GRID, linewidth=0.8)
    ax.set_axisbelow(True)


def legend(fig, names):
    handles = [
        plt.Line2D(
            [], [], ls="", marker=SERIES[n][1], color=SERIES[n][0], ms=7, mec=SURFACE
        )
        for n in names
    ]
    fig.legend(
        handles,
        names,
        loc="upper center",
        ncol=len(names),
        frameon=False,
        fontsize=9,
        labelcolor=INK2,
        bbox_to_anchor=(0.5, 0.94),
    )


def panel_grid(title):
    fig, axes = plt.subplots(2, 2, figsize=(10, 7), sharey=True, facecolor=SURFACE)
    fig.suptitle(title, x=0.06, ha="left", color=INK, fontsize=13)
    return fig, axes.ravel()


def plot_coverage(data, path):
    fig, axes = panel_grid("95% CI coverage by parameter (100 datasets per panel)")
    x = np.arange(len(LABELS))
    offsets = np.linspace(-0.24, 0.24, len(SERIES))
    for ax, key in zip(axes, PANELS):
        res, p_true = data[key]
        n = len(res["FIM"][0])
        half = 1.96 * np.sqrt(0.95 * 0.05 / n)
        ax.axhspan(0.95 - half, min(0.95 + half, 1), color=GRID, alpha=0.6, lw=0)
        ax.axhline(0.95, color=MUTED, lw=1, ls="--")
        for off, (name, (color, marker)) in zip(offsets, SERIES.items()):
            lo, hi, _ = res[name]
            cover = np.mean((lo <= p_true) & (p_true <= hi), axis=0)
            ax.plot(
                x + off, cover, ls="", marker=marker, color=color, ms=8, mec=SURFACE
            )
        style(ax)
        ax.set_title(TITLES[key], loc="left", color=INK, fontsize=10)
        ax.set_xticks(x, LABELS)
        ax.set_ylim(0.3, 1.02)
    for ax in axes[::2]:
        ax.set_ylabel("Coverage", color=INK2, fontsize=9)
    legend(fig, list(SERIES))
    fig.text(
        0.06,
        0.01,
        "Dashed line: nominal 0.95. Band: 95% binomial range for 100 datasets.",
        color=MUTED,
        fontsize=8,
    )
    fig.tight_layout(rect=(0, 0.03, 1, 0.9))
    fig.savefig(path, dpi=150, facecolor=SURFACE)


def plot_sd_ratio(data, path):
    names = ["FIM", "GIM", "Bootstrap SE"]
    fig, axes = panel_grid(
        "Mean standard error relative to the empirical SD of the estimates"
    )
    x = np.arange(len(LABELS))
    offsets = np.linspace(-0.18, 0.18, len(names))
    for ax, key in zip(axes, PANELS):
        res, p_true = data[key]
        emp = data[key + ("est",)].std(axis=0, ddof=1)
        ax.axhline(1, color=MUTED, lw=1, ls="--")
        for off, name in zip(offsets, names):
            color, marker = SERIES[name]
            ratio = res[name][2].mean(axis=0) / emp
            ax.plot(
                x + off, ratio, ls="", marker=marker, color=color, ms=8, mec=SURFACE
            )
        style(ax)
        ax.set_title(TITLES[key], loc="left", color=INK, fontsize=10)
        ax.set_xticks(x, LABELS)
        ax.set_ylim(0, 1.4)
    for ax in axes[::2]:
        ax.set_ylabel("SE / empirical SD", color=INK2, fontsize=9)
    legend(fig, names)
    fig.text(
        0.06,
        0.01,
        "Dashed line: 1 (standard error matches the spread of the estimates).",
        color=MUTED,
        fontsize=8,
    )
    fig.tight_layout(rect=(0, 0.03, 1, 0.9))
    fig.savefig(path, dpi=150, facecolor=SURFACE)


def main():
    data = {}
    for method, scenario in PANELS:
        results = gc.load_results(method, scenario)
        data[(method, scenario)] = (intervals(results), results[0]["p_true"])
        data[(method, scenario, "est")] = np.array([r["p_full"] for r in results])
    os.makedirs("figures", exist_ok=True)
    plot_coverage(data, "figures/coverage.png")
    plot_sd_ratio(data, "figures/sd_ratio.png")


if __name__ == "__main__":
    main()
