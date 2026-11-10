#!/usr/bin/env python3
"""Runs persistence2D on several maps and draws them side by side: map, barcode, persistence diagram and
Betti curves, one column per map.

  python3 scripts/compare_maps.py examples/porous_phi0.40.pgm::porous examples/random_phi0.40.pgm::random \\
      examples/hyperuniform_phi0.40.pgm:light:hyperuniform --out comparison.png

Each map is given as PATH[:OBSTACLE[:LABEL]], OBSTACLE = dark (default) or light. The filtration options are
common to all maps (default: signed filtration, Euclidean widths, periodic, 4-connectivity).

Rows: (1) the map, free sites shaded by corridor width w and obstacles by wall thickness d; (2) barcode of the
most persistent bars, H0 then H1, essential bars as arrows; (3) persistence diagram, marker area proportional
to the number of bars at the same point; (4) Betti curves, with the percolating width r* marked.
Needs NumPy and Matplotlib.
"""
from __future__ import annotations

import argparse
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import p2d  # noqa: E402

H0, H1, WRAP = "#2a78d6", "#eb6834", "#1baf7a"  # categorical slots 1-3 (CVD-checked)
INK, MUTED, GRID = "#1f2328", "#6e7781", "#d0d7de"


def parse_spec(spec: str) -> tuple[Path, str, str]:
    parts = spec.split(":")
    path = Path(parts[0])
    obstacle = parts[1] if len(parts) > 1 and parts[1] else "dark"
    label = parts[2] if len(parts) > 2 and parts[2] else path.stem
    if obstacle not in ("dark", "light"):
        sys.exit(f"{spec}: obstacle must be dark or light")
    return path, obstacle, label


def analyse(binary: Path, path: Path, obstacle: str, a: argparse.Namespace, tmp: Path):
    params = tmp / f"{path.stem}.in"
    params.write_text("\n".join([
        f"obstacle: {obstacle}", f"filtration: {a.filtration}", f"metric: {a.metric}",
        f"periodic: {str(not a.open).lower()}", f"connectivity: {a.connectivity}",
        f"bars_file: {tmp / 'bars.csv'}", f"betti_file: {tmp / 'betti.csv'}", f"summary_file: {tmp / 'summary.csv'}",
        "num_threads: 1", ""]))
    subprocess.run([str(binary), str(params), str(path)], check=True, capture_output=True, text=True)
    bars = next(iter(p2d.read_bars(tmp / "bars.csv").values()))
    betti = next(iter(p2d.read_betti(tmp / "betti.csv").values()))
    summary = p2d.read_summary(tmp / "summary.csv")[0]
    return bars, betti, summary


def draw_map(ax, grid, w, label, summary):
    from matplotlib.colors import ListedColormap  # noqa: F401
    free = np.ma.masked_where(grid, w)
    walls = np.ma.masked_where(~grid, -w)
    ax.imshow(walls, cmap="Greys", vmin=-0.5 * max(walls.max(), 1), vmax=max(walls.max(), 1), interpolation="nearest")
    ax.imshow(free, cmap="Blues", vmin=0, vmax=max(free.max(), 1), interpolation="nearest")
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_title(f"{label}\nφ = {summary['phi']:.3f},  $w_{{\\max}}$ = {summary['w_max']:.3g},  "
                 f"$d_{{\\max}}$ = {summary['d_max']:.3g}", color=INK, fontsize=9)


def draw_barcode(ax, bars, max_bars, r_lo, r_hi):
    y = 0
    ticks, labels = [], []
    for dim, color in ((0, H0), (1, H1)):
        m = bars.dim == dim
        births, deaths, ess = bars.birth[m], bars.death[m], bars.essential[m]
        pers = np.where(ess, np.inf, births - deaths)
        order = np.lexsort((-births, -pers))[:max_bars]  # most persistent first, then latest born
        order = order[np.argsort(-births[order], kind="stable")]  # draw in order of birth
        for i in order:
            end = r_lo if ess[i] else deaths[i]
            ax.plot([births[i], end], [y, y], color=color, lw=1.6, solid_capstyle="butt")
            if ess[i]:
                ax.annotate("", xy=(r_lo, y), xytext=(r_lo + 0.04 * (r_hi - r_lo), y),
                            arrowprops=dict(arrowstyle="->", color=color, lw=1.2))
            y += 1
        shown, total = len(order), int(m.sum())
        ticks.append(y - shown / 2)
        labels.append(f"$H_{dim}$\n{shown} of {total}")
        y += 2
    ax.set_xlim(r_hi, r_lo)
    ax.set_ylim(-1, y - 1)
    ax.set_yticks(ticks, labels, fontsize=8)
    for t, color in zip(ax.get_yticklabels(), (H0, H1)):
        t.set_color(color)
    ax.tick_params(axis="y", length=0)
    ax.spines["left"].set_visible(False)
    ax.set_xlabel("width $r$  (filtration →)")


def draw_diagram(ax, bars, r_lo, r_hi):
    pad = 0.05 * (r_hi - r_lo)
    top = r_hi + 2.5 * pad
    ax.plot([r_lo, r_hi], [r_lo, r_hi], color=MUTED, lw=1, zorder=1)
    ax.axhline(r_hi + 1.5 * pad, color=GRID, lw=0.8, ls=":")
    for dim, color, marker in ((0, H0, "o"), (1, H1, "^")):
        for essential in (False, True):
            d = bars.diagram(dim, essential=essential)
            if not len(d):
                continue
            if essential:
                d = np.column_stack([d[:, 0], np.full(len(d), top)])
            pts, counts = np.unique(np.round(d, 9), axis=0, return_counts=True)
            ax.scatter(pts[:, 0], pts[:, 1], s=14 + 10 * np.sqrt(counts), c=color, marker=marker, alpha=0.6,
                       edgecolors="white", linewidths=0.6, zorder=3,
                       label=None if essential else f"$H_{dim}$ ({int((bars.dim == dim).sum())} bars)")
    ax.text(r_lo, top, " essential", va="center", ha="left", color=MUTED, fontsize=7)
    ax.set_xlim(r_lo - pad, r_hi + pad)
    ax.set_ylim(r_lo - pad, top + pad)
    ax.set_xlabel("birth $r_b$")
    ax.set_ylabel("death $r_d$")
    ax.legend(loc="center left", frameon=False, fontsize=8, markerscale=0.8)


def draw_betti(ax, betti, summary, signed):
    r = betti["r"]
    step = r.max() - r.min()
    rx = np.append(r, r.min() - 0.05 * (step if step > 0 else 1))
    for key, color, name in (("beta0", H0, r"$\beta_0$"), ("beta1", H1, r"$\beta_1$")):
        ax.step(rx, np.append(betti[key], betti[key][-1]), where="post", color=color, lw=1.8, label=name)
    if signed:
        ax.axvspan(rx.min(), 0, color="#eef1f4", zorder=0)
    rs = summary["r_star"]
    if np.isfinite(rs):
        ax.axvline(rs, color=WRAP, lw=1.5, ls="--")
        ax.text(rs, 0.97, f" $r^*$ = {rs:.3g}", transform=ax.get_xaxis_transform(), color=INK, fontsize=8,
                ha="left", va="top")
    ax.set_xlim(rx.max(), rx.min())
    ax.set_xlabel("width $r$  (filtration →)")
    ax.set_ylabel("Betti number")
    ax.legend(loc="best", frameon=False, fontsize=8)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("maps", nargs="+", help="PATH[:OBSTACLE[:LABEL]]")
    ap.add_argument("--out", type=Path, default=Path("comparison.png"))
    ap.add_argument("--binary", type=Path, default=HERE.parent / "build" / "persistence2D")
    ap.add_argument("--filtration", choices=("free", "signed"), default="signed")
    ap.add_argument("--metric", choices=("l1", "l2"), default="l2")
    ap.add_argument("--open", action="store_true", help="open boundaries (default: periodic)")
    ap.add_argument("--connectivity", type=int, choices=(4, 8), default=4)
    ap.add_argument("--max-bars", type=int, default=40, help="bars per dimension in the barcode")
    a = ap.parse_args()

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    plt.rcParams.update({"font.size": 9, "axes.edgecolor": MUTED, "axes.labelcolor": INK, "xtick.color": MUTED,
                         "ytick.color": MUTED, "axes.spines.top": False, "axes.spines.right": False})
    specs = [parse_spec(s) for s in a.maps]
    n = len(specs)
    fig, ax = plt.subplots(4, n, figsize=(3.6 * n, 13), squeeze=False, constrained_layout=True,
                           gridspec_kw={"height_ratios": [1.1, 1, 1, 0.8]})
    with tempfile.TemporaryDirectory(prefix="p2d_cmp_") as tmp:
        for j, (path, obstacle, label) in enumerate(specs):
            bars, betti, summary = analyse(a.binary, path, obstacle, a, Path(tmp))
            grid = p2d.read_map(path, obstacle)
            w = p2d.width_field(grid, not a.open, a.metric, signed=a.filtration == "signed")
            r_lo, r_hi = float(betti["r"].min()), float(betti["r"].max())
            draw_map(ax[0, j], grid, w, label, summary)
            draw_barcode(ax[1, j], bars, a.max_bars, r_lo, r_hi)
            draw_diagram(ax[2, j], bars, r_lo, r_hi)
            draw_betti(ax[3, j], betti, summary, a.filtration == "signed")
    for j in range(n):
        for i, title in ((1, "barcode"), (2, "persistence diagram"), (3, "Betti curves")):
            ax[i, j].set_title(title, loc="left", color=INK, fontsize=9)
    fig.suptitle(f"{a.filtration} filtration, {a.metric.upper()} widths, {'open' if a.open else 'periodic'}, "
                 f"{a.connectivity}-connectivity; barcodes show the {a.max_bars} most persistent bars per dimension" + ("  (shaded: inside walls, r < 0)" if a.filtration == "signed" else ""),
                 color=INK)
    fig.savefig(a.out, dpi=130, facecolor="white")
    print(a.out)
    return 0


if __name__ == "__main__":
    sys.exit(main())
