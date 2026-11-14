"""
Plots of 95% CI coverage and standard-error calibration for the second model
(growth_asym_calibration.py, scenario g_asym, the one with bootstrap refits), in the
same style as plot_coverage.py. Reads results_growth/ in the current directory (seed
9001, the cluster test, is excluded) and writes figures/coverage_growth.png and
figures/sd_ratio_growth.png.

Usage: python plot_coverage_growth.py
"""

import glob
import os
import pickle

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

import plot_coverage as pc

LABELS = ["nu0", "nu1_0", "nu1_F", "T", "M01", "M10", "Ne"]
PANELS = [("ld", "LD"), ("sfs", "SFS")]
SUBTITLE = "deme1 grows, M01 = 1, M10 = 0.2"


def load(method):
    out = []
    for path in sorted(glob.glob(f"results_growth/{method}_g_asym_*.pkl")):
        if not path.endswith("_9001.pkl"):
            with open(path, "rb") as fin:
                out.append(pickle.load(fin))
    return out


def intervals(results):
    """Per-dataset 95% intervals and standard errors for each method."""
    out = {name: ([], [], []) for name in pc.SERIES}
    for r in results:
        p = r["p_full"]
        boot_sd = np.nanstd(r["p_boot"], axis=0, ddof=1)
        for name in ["FIM", "GIM"]:
            sd = r["sd"][name]
            out[name][0].append(p - 1.96 * sd)
            out[name][1].append(p + 1.96 * sd)
            out[name][2].append(sd)
        out["Bootstrap SE"][0].append(p - 1.96 * boot_sd)
        out["Bootstrap SE"][1].append(p + 1.96 * boot_sd)
        out["Bootstrap SE"][2].append(boot_sd)
        lo, hi = np.nanpercentile(r["p_boot"], [2.5, 97.5], axis=0)
        out["Bootstrap percentile"][0].append(lo)
        out["Bootstrap percentile"][1].append(hi)
        out["Bootstrap percentile"][2].append(boot_sd)
    return {k: tuple(np.array(x) for x in v) for k, v in out.items()}


def panels(title):
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2), sharey=True, facecolor=pc.SURFACE)
    fig.suptitle(title, x=0.06, ha="left", color=pc.INK, fontsize=13)
    return fig, axes


def plot_coverage(data, path):
    fig, axes = panels("95% CI coverage by parameter, second model (100 datasets)")
    x = np.arange(len(LABELS))
    offsets = np.linspace(-0.27, 0.27, len(pc.SERIES))
    for ax, (method, label) in zip(axes, PANELS):
        res, p_true, _ = data[method]
        n = len(res["FIM"][0])
        half = 1.96 * np.sqrt(0.95 * 0.05 / n)
        ax.axhspan(0.95 - half, min(0.95 + half, 1), color=pc.GRID, alpha=0.6, lw=0)
        ax.axhline(0.95, color=pc.MUTED, lw=1, ls="--")
        for off, (name, (color, marker)) in zip(offsets, pc.SERIES.items()):
            lo, hi, _ = res[name]
            cover = np.mean((lo <= p_true) & (p_true <= hi), axis=0)
            ax.plot(
                x + off, cover, ls="", marker=marker, color=color, ms=7, mec=pc.SURFACE
            )
        pc.style(ax)
        ax.set_title(f"{label} ({SUBTITLE})", loc="left", color=pc.INK, fontsize=10)
        ax.set_xticks(x, LABELS)
        ax.set_ylim(0.2, 1.02)
    axes[0].set_ylabel("Coverage", color=pc.INK2, fontsize=9)
    pc.legend(fig, list(pc.SERIES))
    fig.text(
        0.06,
        0.01,
        "Dashed line: nominal 0.95. Band: 95% binomial range for 100 datasets.",
        color=pc.MUTED,
        fontsize=8,
    )
    fig.tight_layout(rect=(0, 0.04, 1, 0.86))
    fig.savefig(path, dpi=150, facecolor=pc.SURFACE)


def plot_sd_ratio(data, path):
    names = ["FIM", "GIM", "Bootstrap SE"]
    fig, axes = panels("Mean standard error relative to the empirical SD, second model")
    x = np.arange(len(LABELS))
    offsets = np.linspace(-0.2, 0.2, len(names))
    for ax, (method, label) in zip(axes, PANELS):
        res, p_true, est = data[method]
        emp = est.std(axis=0, ddof=1)
        ax.axhline(1, color=pc.MUTED, lw=1, ls="--")
        for off, name in zip(offsets, names):
            color, marker = pc.SERIES[name]
            ratio = np.nanmean(res[name][2], axis=0) / emp
            ax.plot(
                x + off, ratio, ls="", marker=marker, color=color, ms=7, mec=pc.SURFACE
            )
        pc.style(ax)
        ax.set_title(f"{label} ({SUBTITLE})", loc="left", color=pc.INK, fontsize=10)
        ax.set_xticks(x, LABELS)
        ax.set_ylim(0, 1.4)
    axes[0].set_ylabel("SE / empirical SD", color=pc.INK2, fontsize=9)
    pc.legend(fig, names)
    fig.text(
        0.06,
        0.01,
        "Dashed line: 1 (standard error matches the spread of the estimates).",
        color=pc.MUTED,
        fontsize=8,
    )
    fig.tight_layout(rect=(0, 0.04, 1, 0.86))
    fig.savefig(path, dpi=150, facecolor=pc.SURFACE)


def main():
    data = {}
    for method, _ in PANELS:
        results = load(method)
        est = np.array([r["p_full"] for r in results])
        data[method] = (intervals(results), results[0]["p_true"], est)
    os.makedirs("figures", exist_ok=True)
    plot_coverage(data, "figures/coverage_growth.png")
    plot_sd_ratio(data, "figures/sd_ratio_growth.png")


if __name__ == "__main__":
    main()
