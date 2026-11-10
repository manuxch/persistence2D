#!/usr/bin/env python3
"""Runs persistence2D on one map and draws a four-panel summary figure.

  python3 scripts/plot_persistence.py examples/random_phi0.40.pgm --out random.png
  python3 scripts/plot_persistence.py examples/hyperuniform_phi0.40.pgm --obstacle light --filtration signed

Panels: (a) the map, free sites shaded by corridor width w; (b) persistence diagram, birth against death
(superlevel filtration, so the points lie below the diagonal; essential classes on the top strip);
(c) Betti curves beta0(r), beta1(r); (d) wrap dimension of V_r. Needs NumPy and Matplotlib.
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
INK, MUTED = "#1f2328", "#6e7781"


def run(binary: Path, map_file: Path, a: argparse.Namespace, out: Path) -> None:
    params = out / "plot.in"
    params.write_text("\n".join([
        f"obstacle: {a.obstacle}", f"filtration: {a.filtration}", f"metric: {a.metric}",
        f"periodic: {str(not a.open).lower()}", f"boundary_obstacle: {str(a.boundary_obstacle).lower()}",
        f"connectivity: {a.connectivity}", f"bars_file: {out / 'bars.csv'}", f"betti_file: {out / 'betti.csv'}",
        f"summary_file: {out / 'summary.csv'}", "num_threads: 1", ""]))
    subprocess.run([str(binary), str(params), str(map_file)], check=True, capture_output=True, text=True)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("map", type=Path)
    ap.add_argument("--out", type=Path, help="figure file (default: <map stem>_persistence.png)")
    ap.add_argument("--binary", type=Path, default=HERE.parent / "build" / "persistence2D")
    ap.add_argument("--obstacle", choices=("dark", "light"), default="dark")
    ap.add_argument("--filtration", choices=("free", "signed"), default="free")
    ap.add_argument("--metric", choices=("l1", "l2"), default="l1")
    ap.add_argument("--open", action="store_true", help="open boundaries (default: periodic)")
    ap.add_argument("--boundary-obstacle", action="store_true")
    ap.add_argument("--connectivity", type=int, choices=(4, 8), default=4)
    ap.add_argument("--title", default=None)
    a = ap.parse_args()

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import ListedColormap

    with tempfile.TemporaryDirectory(prefix="p2d_plot_") as tmp:
        tmp = Path(tmp)
        run(a.binary, a.map, a, tmp)
        bars = next(iter(p2d.read_bars(tmp / "bars.csv").values()))
        betti = next(iter(p2d.read_betti(tmp / "betti.csv").values()))
        summary = p2d.read_summary(tmp / "summary.csv")[0]

    grid = p2d.read_map(a.map, a.obstacle)
    w = p2d.width_field(grid, not a.open, a.metric) if grid.size <= 256 * 256 else None

    plt.rcParams.update({"font.size": 9, "axes.edgecolor": MUTED, "axes.labelcolor": INK, "xtick.color": MUTED,
                         "ytick.color": MUTED, "axes.spines.top": False, "axes.spines.right": False})
    fig, ax = plt.subplots(2, 2, figsize=(9, 8), constrained_layout=True)

    # (a) map
    if w is not None:
        shade = np.ma.masked_where(grid, w)
        ax[0, 0].imshow(np.zeros_like(w), cmap=ListedColormap(["#1f2328"]), interpolation="nearest")
        im = ax[0, 0].imshow(shade, cmap="Blues", vmin=0, interpolation="nearest")
        fig.colorbar(im, ax=ax[0, 0], shrink=0.8, label="corridor width $w$")
    else:
        ax[0, 0].imshow(grid, cmap="gray_r", interpolation="nearest")
    ax[0, 0].set_title(f"(a) map {grid.shape[1]} × {grid.shape[0]}, φ = {summary['phi']:.3f} (obstacles black)",
                       loc="left", color=INK)
    ax[0, 0].set_xticks([])
    ax[0, 0].set_yticks([])

    # (b) persistence diagram
    b = ax[0, 1]
    finite_all = np.concatenate([bars.diagram(0), bars.diagram(1)])
    lo = min(betti["r"].min(), finite_all[:, 1].min() if len(finite_all) else betti["r"].min())
    hi = betti["r"].max()
    pad = 0.06 * (hi - lo + 1)
    top = hi + 3 * pad
    b.plot([lo - pad, hi + pad], [lo - pad, hi + pad], color=MUTED, lw=1, zorder=1)
    b.axhline(hi + 2 * pad, color=MUTED, lw=0.6, ls=":")
    rng = np.random.default_rng(0)
    for dim, color, marker in ((0, H0, "o"), (1, H1, "^")):
        d = bars.diagram(dim)
        e = bars.diagram(dim, essential=True)
        jit = lambda n: rng.uniform(-0.12, 0.12, n) * (1 if a.metric == "l1" else 0)  # noqa: E731
        b.scatter(d[:, 0] + jit(len(d)), d[:, 1] + jit(len(d)), s=22, c=color, marker=marker, alpha=0.55,
                  edgecolors="white", linewidths=0.6, label=f"$H_{dim}$: {len(d)} finite, {len(e)} essential", zorder=3)
        b.scatter(e[:, 0] + jit(len(e)), np.full(len(e), top) + jit(len(e)), s=22, c=color, marker=marker,
                  alpha=0.55, edgecolors="white", linewidths=0.6, zorder=3)
    b.text(lo - pad, top, "essential", va="center", ha="left", color=MUTED, fontsize=8)
    b.set_xlim(lo - pad, hi + pad)
    b.set_ylim(lo - pad, top + pad)
    b.set_xlabel("birth width $r_b$")
    b.set_ylabel("death width $r_d$")
    b.set_title("(b) persistence diagram" + (" (points jittered)" if a.metric == "l1" else ""), loc="left", color=INK)
    b.legend(loc="center left", frameon=False)

    # (c) Betti curves and (d) wrap dimension, filtration progressing to the right (r decreasing)
    # Extend the last level to the right so that its plateau is drawn.
    step = (r_all := betti["r"]).max() - r_all.min()
    r = np.append(r_all, r_all.min() - 0.08 * (step if step > 0 else 1))
    ext = {k: np.append(betti[k], betti[k][-1]) for k in ("beta0", "beta1", "wraps")}
    shaded = " (shaded: inside walls, r < 0)" if a.filtration == "signed" else ""
    for axis in (ax[1, 0], ax[1, 1]):
        axis.invert_xaxis()
        axis.set_xlabel("width $r$  (filtration →)")
        if a.filtration == "signed":
            axis.axvspan(r.min(), 0, color="#eaeef2", zorder=0)
    c = ax[1, 0]
    c.step(r, ext["beta0"], where="post", color=H0, lw=2, label=r"$\beta_0$ (components)")
    c.step(r, ext["beta1"], where="post", color=H1, lw=2, label=r"$\beta_1$ (loops)")
    c.scatter(r_all, betti["beta0"], color=H0, s=14, zorder=3)
    c.scatter(r_all, betti["beta1"], color=H1, s=14, marker="^", zorder=3)
    c.set_ylim(top=1.22 * max(ext["beta0"].max(), ext["beta1"].max(), 1))  # headroom for the legend
    c.set_ylabel("Betti number of $V_r$")
    c.set_title("(c) Betti curves" + shaded, loc="left", color=INK)
    c.legend(frameon=False, loc="upper left")
    d = ax[1, 1]
    d.step(r, ext["wraps"], where="post", color=WRAP, lw=2)
    d.scatter(r_all, betti["wraps"], color=WRAP, s=14, zorder=3)
    d.set_yticks([0, 1, 2])
    d.set_ylim(-0.2, 2.3)
    d.set_ylabel("independent wrapping loops")
    rs = summary["r_star"]
    d.set_title(f"(d) wrap dimension;  $r^*$ = {rs:g}" if np.isfinite(rs) else "(d) wrap dimension;  $r^*$ undefined",
                loc="left", color=INK)

    fig.suptitle(a.title or f"{a.map.name}: {a.filtration} filtration, {a.metric.upper()}, "
                 f"{'open' if a.open else 'periodic'}, {a.connectivity}-connectivity", color=INK)
    out = a.out or Path(f"{a.map.stem}_persistence.png")
    fig.savefig(out, dpi=150, facecolor="white")
    print(out)
    return 0


if __name__ == "__main__":
    sys.exit(main())
