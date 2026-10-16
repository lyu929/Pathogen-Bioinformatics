"""Figures for the report.

Colours follow a validated, colour-vision-deficiency-safe palette: categorical
slots are used in a fixed order, magnitude uses a single-hue blue ramp, and
polarity (up/down regulation) uses a blue <-> red diverging pair with a neutral
grey for "not significant". Every figure is rendered for a light and a dark
surface so the README can serve the right one via ``<picture>``.
"""

from __future__ import annotations

from collections.abc import Callable, Iterator, Mapping
from contextlib import contextmanager
from dataclasses import dataclass
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.colors import LinearSegmentedColormap, ListedColormap  # noqa: E402


@dataclass(frozen=True)
class Theme:
    name: str
    surface: str
    ink: str
    ink_secondary: str
    muted: str
    grid: str
    axis: str
    neutral: str
    series: tuple[str, ...]
    blue_ramp: tuple[str, ...]
    up: str
    down: str


LIGHT = Theme(
    name="light",
    surface="#fcfcfb",
    ink="#0b0b0b",
    ink_secondary="#52514e",
    muted="#898781",
    grid="#e1e0d9",
    axis="#c3c2b7",
    neutral="#c3c2b7",
    series=("#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6", "#e34948"),
    blue_ramp=(
        "#f0efec",
        "#cde2fb",
        "#9ec5f4",
        "#6da7ec",
        "#3987e5",
        "#256abf",
        "#184f95",
        "#0d366b",
    ),
    up="#e34948",
    down="#2a78d6",
)
DARK = Theme(
    name="dark",
    surface="#1a1a19",
    ink="#f0efec",
    ink_secondary="#c3c2b7",
    muted="#898781",
    grid="#2c2c2a",
    axis="#383835",
    neutral="#5c5b57",
    series=("#3987e5", "#d95926", "#199e70", "#c98500", "#d55181", "#008300", "#9085e9", "#e66767"),
    blue_ramp=(
        "#383835",
        "#104281",
        "#184f95",
        "#1c5cab",
        "#256abf",
        "#3987e5",
        "#6da7ec",
        "#b7d3f6",
    ),
    up="#e66767",
    down="#3987e5",
)
THEMES = (LIGHT, DARK)


@contextmanager
def themed(theme: Theme) -> Iterator[None]:
    rc = {
        "figure.facecolor": theme.surface,
        "axes.facecolor": theme.surface,
        "savefig.facecolor": theme.surface,
        "axes.edgecolor": theme.axis,
        "axes.labelcolor": theme.ink_secondary,
        "axes.titlecolor": theme.ink,
        "axes.titleweight": "bold",
        "axes.titlesize": 11,
        "axes.titlelocation": "left",
        "axes.labelsize": 9.5,
        "axes.grid": True,
        "axes.axisbelow": True,
        "grid.color": theme.grid,
        "grid.linewidth": 0.6,
        "grid.linestyle": "-",
        "axes.spines.top": False,
        "axes.spines.right": False,
        "xtick.color": theme.muted,
        "ytick.color": theme.muted,
        "xtick.labelcolor": theme.ink_secondary,
        "ytick.labelcolor": theme.ink_secondary,
        "xtick.labelsize": 8.5,
        "ytick.labelsize": 8.5,
        "legend.frameon": False,
        "legend.fontsize": 8.5,
        "legend.labelcolor": theme.ink_secondary,
        "text.color": theme.ink,
        "font.family": "sans-serif",
        "font.sans-serif": ["DejaVu Sans", "Helvetica", "Arial"],
        "lines.linewidth": 2.0,
        "axes.prop_cycle": matplotlib.cycler(color=list(theme.series)),
    }
    with matplotlib.rc_context(rc):
        yield


def sequential_cmap(theme: Theme) -> LinearSegmentedColormap:
    return LinearSegmentedColormap.from_list(f"blue_{theme.name}", list(theme.blue_ramp))


def save_both(
    draw: Callable[[Theme], plt.Figure], out_dir: str | Path, stem: str
) -> dict[str, Path]:
    """Render ``draw`` for each theme to ``<out_dir>/<theme>/<stem>.png``."""
    paths = {}
    for theme in THEMES:
        with themed(theme):
            fig = draw(theme)
            target = Path(out_dir) / theme.name / f"{stem}.png"
            target.parent.mkdir(parents=True, exist_ok=True)
            fig.savefig(target, dpi=150, bbox_inches="tight", metadata={"Software": None})
            plt.close(fig)
        paths[theme.name] = target
    return paths


# --------------------------------------------------------------------------- QC
def per_base_quality(
    profiles: Mapping[str, np.ndarray],
    counts: Mapping[str, np.ndarray],
    cutoff: int,
    min_fraction: float = 0.5,
) -> Callable[[Theme], plt.Figure]:
    """Mean quality per position; each curve stops where fewer than ``min_fraction`` of
    its reads still reach that position (later means describe a survivor subset)."""

    def draw(t: Theme) -> plt.Figure:
        fig, ax = plt.subplots(figsize=(6.4, 3.2))
        for i, (label, prof) in enumerate(profiles.items()):
            n = counts[label]
            keep = n >= min_fraction * n.max()
            x = np.arange(1, len(prof) + 1)[keep]
            ax.plot(x, prof[keep], color=t.series[i], label=label)
        ax.axhline(cutoff, color=t.axis, linewidth=1)
        ax.text(
            ax.get_xlim()[1],
            cutoff + 0.3,
            f"Q{cutoff} cutoff",
            color=t.muted,
            fontsize=8,
            va="bottom",
            ha="right",
        )
        ax.set_ylim(cutoff - 2, None)
        ax.set_xlabel("Position in read (bp)")
        ax.set_ylabel("Mean Phred quality")
        ax.set_title("Per-base quality of BR-Tn1 reads, raw vs Q20-trimmed")
        ax.legend(loc="upper right")
        return fig

    return draw


# --------------------------------------------------------------------------- assembly
def nx_curves(lengths: Mapping[str, np.ndarray]) -> Callable[[Theme], plt.Figure]:
    def draw(t: Theme) -> plt.Figure:
        fig, ax = plt.subplots(figsize=(6.4, 3.4))
        for i, (label, arr) in enumerate(lengths.items()):
            ordered = np.sort(arr)[::-1]
            frac = np.concatenate([[0], np.cumsum(ordered) / ordered.sum() * 100])
            y = np.concatenate([[ordered[0]], ordered]) / 1e3
            ax.step(frac, y, where="post", color=t.series[i], label=label)
        ax.axvline(50, color=t.axis, linewidth=1)
        ax.text(51, ax.get_ylim()[1] * 0.92, "N50", color=t.muted, fontsize=8)
        ax.set_yscale("log")
        ax.set_xlabel("Cumulative assembly length (%)")
        ax.set_ylabel("Contig length (kb, log)")
        ax.set_title("Nx curves of SPAdes assemblies (P. aeruginosa isolates)")
        ax.legend(loc="lower left", ncol=2)
        return fig

    return draw


# --------------------------------------------------------------------------- presence heatmaps
def presence_heatmap(
    presence: pd.DataFrame,
    title: str,
    xlabel: str = "",
    cell_labels: bool = False,
    max_width: float = 9.0,
) -> Callable[[Theme], plt.Figure]:
    def draw(t: Theme) -> plt.Figure:
        rows, cols = presence.shape
        width = min(max_width, max(5.5, 0.42 * cols + 2.5))
        fig, ax = plt.subplots(figsize=(width, max(2.2, 0.36 * rows + 1.2)))
        cmap = ListedColormap([t.blue_ramp[0], t.series[0]])
        ax.imshow(
            presence.to_numpy(dtype=float),
            cmap=cmap,
            aspect="auto",
            vmin=0,
            vmax=1,
            interpolation="nearest",
        )
        show_x = cols <= 30
        if show_x:
            ax.set_xticks(range(cols), presence.columns, rotation=45, ha="right")
            ax.set_xticks(np.arange(-0.5, cols), minor=True)
        else:
            ax.set_xticks([])
        ax.set_yticks(range(rows), presence.index)
        ax.set_yticks(np.arange(-0.5, rows), minor=True)
        ax.grid(which="minor", color=t.surface, linewidth=2)
        ax.grid(which="major", visible=False)
        ax.tick_params(which="both", length=0)
        for spine in ax.spines.values():
            spine.set_visible(False)
        if cell_labels:
            for r in range(rows):
                for c in range(cols):
                    if presence.iat[r, c]:
                        ax.text(
                            c,
                            r,
                            "+",
                            ha="center",
                            va="center",
                            color=t.surface,
                            fontsize=9,
                            fontweight="bold",
                        )
        ax.set_title(title)
        if xlabel:
            ax.set_xlabel(xlabel)
        return fig

    return draw


# --------------------------------------------------------------------------- pangenome
def rarefaction_plot(curve: pd.DataFrame) -> Callable[[Theme], plt.Figure]:
    def draw(t: Theme) -> plt.Figure:
        fig, ax = plt.subplots(figsize=(5.6, 3.3))
        x = curve["genomes"]
        for i, (key, label) in enumerate([("pan", "Pangenome"), ("core", "Core genome")]):
            mean, sd = curve[f"{key}_mean"], curve[f"{key}_sd"]
            ax.fill_between(x, mean - sd, mean + sd, color=t.series[i], alpha=0.18, linewidth=0)
            ax.plot(x, mean, color=t.series[i], marker="o", markersize=5, label=label)
            ax.annotate(
                f"{label}: {mean.iloc[-1]:,.0f}",
                (x.iloc[-1], mean.iloc[-1]),
                xytext=(6, 0),
                textcoords="offset points",
                va="center",
                fontsize=8.5,
                color=t.ink_secondary,
            )
        ax.set_xlim(0.8, x.max() + 1.4)
        ax.set_xticks(list(x))
        ax.set_xlabel("Genomes sampled")
        ax.set_ylabel("Gene clusters")
        ax.set_title("Pangenome rarefaction (mean +/- SD over permutations)")
        ax.legend(loc="center right", bbox_to_anchor=(1.0, 0.5))
        return fig

    return draw


# --------------------------------------------------------------------------- RNA-seq
def volcano(
    results: pd.DataFrame,
    title: str,
    alpha: float,
    min_lfc: float,
    labels: Mapping[str, str] | None = None,
) -> Callable[[Theme], plt.Figure]:
    def draw(t: Theme) -> plt.Figure:
        fig, ax = plt.subplots(figsize=(5.8, 4.2))
        res = results.dropna(subset=["log2FoldChange", "padj"])
        y = -np.log10(res["padj"].clip(lower=1e-300))
        sig = (res["padj"] < alpha) & (res["log2FoldChange"].abs() >= min_lfc)
        up = sig & (res["log2FoldChange"] > 0)
        down = sig & (res["log2FoldChange"] < 0)
        ax.scatter(
            res.loc[~sig, "log2FoldChange"],
            y[~sig],
            s=7,
            color=t.neutral,
            linewidths=0,
            label="Not significant",
        )
        ax.scatter(
            res.loc[up, "log2FoldChange"],
            y[up],
            s=14,
            color=t.up,
            edgecolors=t.surface,
            linewidths=0.6,
            label=f"Up ({int(up.sum())})",
        )
        ax.scatter(
            res.loc[down, "log2FoldChange"],
            y[down],
            s=14,
            color=t.down,
            edgecolors=t.surface,
            linewidths=0.6,
            label=f"Down ({int(down.sum())})",
        )
        for v in (-min_lfc, min_lfc):
            ax.axvline(v, color=t.axis, linewidth=1)
        ax.axhline(-np.log10(alpha), color=t.axis, linewidth=1)
        if labels:
            _label_points(
                ax, res[sig].assign(score=y[sig]).nlargest(8, "score"), labels, t, max_labels=5
            )
        ax.set_xlabel("log2 fold change")
        ax.set_ylabel("-log10 adjusted p")
        ax.set_title(title)
        ax.legend(loc="upper left", markerscale=1.4)
        return fig

    return draw


def _label_points(
    ax, top: pd.DataFrame, labels: Mapping[str, str], t: Theme, max_labels: int = 5
) -> None:
    """Greedy, collision-avoiding direct labels for the most significant points."""
    ax.figure.canvas.draw()
    placed: list = []
    renderer = ax.figure.canvas.get_renderer()
    for gid, row in top.iterrows():
        if len(placed) >= max_labels:
            break
        for dx, dy in ((5, 3), (5, -10), (-5, 3)):
            ann = ax.annotate(
                labels.get(gid) or gid,
                (row["log2FoldChange"], row["score"]),
                xytext=(dx, dy),
                textcoords="offset points",
                ha="left" if dx > 0 else "right",
                fontsize=7.5,
                color=t.ink_secondary,
            )
            box = ann.get_window_extent(renderer).expanded(1.05, 1.15)
            inside = ax.get_window_extent(renderer).contains(
                box.x0, box.y0
            ) and ax.get_window_extent(renderer).contains(box.x1, box.y1)
            if inside and not any(box.overlaps(b) for b in placed):
                placed.append(box)
                break
            ann.remove()


def lfc_concordance(
    ours: pd.DataFrame, reference: pd.DataFrame, title: str, r: float
) -> Callable[[Theme], plt.Figure]:
    def draw(t: Theme) -> plt.Figure:
        common = ours.index.intersection(reference.index)
        a = ours.loc[common, "log2FoldChange"]
        b = reference.loc[common, "log2FoldChange"]
        ok = a.notna() & b.notna()
        fig, ax = plt.subplots(figsize=(4.4, 4.2))
        lim = float(np.nanmax(np.abs(np.concatenate([a[ok], b[ok]])))) * 1.05
        ax.plot([-lim, lim], [-lim, lim], color=t.axis, linewidth=1)
        ax.scatter(b[ok], a[ok], s=6, color=t.series[0], alpha=0.55, linewidths=0)
        ax.set_xlim(-lim, lim)
        ax.set_ylim(-lim, lim)
        ax.set_aspect("equal")
        ax.set_xlabel("log2FC, DESeq2 in Galaxy (R)")
        ax.set_ylabel("log2FC, PyDESeq2 (this repo)")
        ax.set_title(title)
        ax.text(
            0.04,
            0.95,
            f"Pearson r = {r:.3f}",
            transform=ax.transAxes,
            va="top",
            fontsize=8.5,
            color=t.ink_secondary,
        )
        return fig

    return draw


def pca_plot(
    coords: pd.DataFrame, groups: pd.Series, explained: np.ndarray, title: str
) -> Callable[[Theme], plt.Figure]:
    def draw(t: Theme) -> plt.Figure:
        fig, ax = plt.subplots(figsize=(4.6, 3.6))
        for i, level in enumerate(sorted(groups.unique())):
            idx = groups.index[groups == level]
            ax.scatter(
                coords.loc[idx, "PC1"],
                coords.loc[idx, "PC2"],
                s=60,
                color=t.series[i],
                edgecolors=t.surface,
                linewidths=2,
                label=level,
            )
            for s in idx:
                ax.annotate(
                    s,
                    (coords.at[s, "PC1"], coords.at[s, "PC2"]),
                    xytext=(5, 4),
                    textcoords="offset points",
                    fontsize=7,
                    color=t.muted,
                )
        ax.set_xlabel(f"PC1 ({explained[0] * 100:.0f}% variance)")
        ax.set_ylabel(f"PC2 ({explained[1] * 100:.0f}% variance)")
        ax.set_title(title)
        ax.legend(loc="best")
        ax.margins(0.25)
        return fig

    return draw


# --------------------------------------------------------------------------- structure
def plddt_and_pae(
    resnums: np.ndarray, plddt: np.ndarray, pae: np.ndarray, title: str
) -> Callable[[Theme], plt.Figure]:
    def draw(t: Theme) -> plt.Figure:
        fig, (ax1, ax2) = plt.subplots(
            1, 2, figsize=(9.2, 3.6), gridspec_kw={"width_ratios": [1.5, 1]}
        )
        ax1.plot(resnums, plddt, color=t.series[0])
        for level, name in ((90, "very high"), (70, "confident"), (50, "low")):
            ax1.axhline(level, color=t.axis, linewidth=1)
            ax1.text(resnums.max() + 1, level, name, color=t.muted, fontsize=7.5, va="center")
        ax1.set_ylim(0, 100)
        ax1.set_xlabel("Residue")
        ax1.set_ylabel("pLDDT")
        ax1.set_title(title)
        cmap = LinearSegmentedColormap.from_list("pae", list(reversed(t.blue_ramp)))
        im = ax2.imshow(pae, cmap=cmap, vmin=0, vmax=max(30.0, float(pae.max())), origin="upper")
        ax2.grid(False)
        ax2.set_xlabel("Scored residue")
        ax2.set_ylabel("Aligned residue")
        ax2.set_title("Predicted aligned error")
        cb = fig.colorbar(im, ax=ax2, fraction=0.046, pad=0.04)
        cb.set_label("PAE (Angstrom)", color=t.ink_secondary, fontsize=8.5)
        cb.outline.set_visible(False)
        cb.ax.tick_params(colors=t.muted, labelcolor=t.ink_secondary, labelsize=8)
        fig.tight_layout()
        return fig

    return draw


def plddt_profiles(
    profiles: Mapping[str, tuple[np.ndarray, np.ndarray]], title: str
) -> Callable[[Theme], plt.Figure]:
    def draw(t: Theme) -> plt.Figure:
        fig, ax = plt.subplots(figsize=(6.4, 3.0))
        for i, (label, (x, y)) in enumerate(profiles.items()):
            ax.plot(x, y, color=t.series[i], label=label)
        ax.axhline(70, color=t.axis, linewidth=1)
        ax.set_ylim(0, 100)
        ax.set_xlabel("Residue")
        ax.set_ylabel("pLDDT")
        ax.set_title(title)
        ax.legend(loc="lower center", ncol=2)
        return fig

    return draw
