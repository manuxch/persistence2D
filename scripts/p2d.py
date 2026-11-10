"""Python helpers for persistence2D: read maps and the CSV outputs, and vectorize persistence diagrams.

Needs only NumPy. Import it from another directory with
    import sys; sys.path.insert(0, "/path/to/persistence2D/scripts"); import p2d

    grid = p2d.read_map("examples/random_phi0.40.pgm")           # bool array [y, x], True = obstacle
    bars = p2d.read_bars("results/bars.csv")                       # {file: Bars}
    b = bars["examples/random_phi0.40.pgm"]
    b.diagram(1)                                                   # (n, 2) finite (birth, death) pairs of H1
    img = p2d.persistence_image(b.diagram(0), grid=p2d.image_grid(1, 4, 0, 3, 20), sigma=0.25)
"""
from __future__ import annotations

import csv
import math
from dataclasses import dataclass
from pathlib import Path

import numpy as np

# ---------------------------------------------------------------------------------------------------- maps


def _pnm_tokens(data: bytes, start: int, count: int) -> tuple[list[int], int]:
    """`count` integers of a PNM header from byte `start`, skipping whitespace and '#' comments."""
    out, i = [], start
    while len(out) < count:
        while i < len(data) and (data[i:i + 1].isspace() or data[i:i + 1] == b"#"):
            if data[i:i + 1] == b"#":
                while i < len(data) and data[i:i + 1] != b"\n":
                    i += 1
            else:
                i += 1
        j = i
        while j < len(data) and data[j:j + 1].isdigit():
            j += 1
        if j == i:
            raise ValueError("malformed PNM header")
        out.append(int(data[i:j]))
        i = j
    return out, i


def read_map(path: str | Path, obstacle: str = "dark") -> np.ndarray:
    """Reads any persistence2D input format; returns a bool array [y, x], True = obstacle.

    `obstacle` ("dark" or "light") says which PBM/PGM pixels are obstacles, as the key of the same name.
    """
    data = Path(path).read_bytes()
    magic = data[:2]
    if magic in (b"P1", b"P2", b"P4", b"P5"):
        pbm = magic in (b"P1", b"P4")
        (nx, ny), i = _pnm_tokens(data, 2, 2)
        maxval = 1
        if not pbm:
            (maxval,), i = _pnm_tokens(data, i, 1)
        if magic == b"P1":
            digits = [c for c in data[i:].decode() if c in "01"]
            dark = np.array(digits[:nx * ny], dtype=int).reshape(ny, nx) == 1
        elif magic == b"P2":
            vals = np.array(data[i:].split()[:nx * ny], dtype=int).reshape(ny, nx)
            dark = 2 * vals <= maxval
        elif magic == b"P4":
            row = (nx + 7) // 8
            raw = np.frombuffer(data, np.uint8, row * ny, i + 1).reshape(ny, row)
            dark = np.unpackbits(raw, axis=1)[:, :nx] == 1
        else:
            dtype = np.uint8 if maxval < 256 else ">u2"
            vals = np.frombuffer(data, dtype, nx * ny, i + 1).reshape(ny, nx)
            dark = 2 * vals.astype(int) <= maxval
        return dark if obstacle == "dark" else ~dark
    text = data.decode()
    if text.startswith("# chasers-map"):
        text = text.split("\n---\n", 1)[1]
    rows = []
    for line in text.splitlines():
        if line.startswith("%"):
            continue
        row = [c for c in line if not c.isspace()]
        if row:
            rows.append([c in "#1" for c in row])
    return np.array(rows, dtype=bool)


# ---------------------------------------------------------------------------------------------- CSV files


def _rows(path: str | Path):
    with open(path, newline="") as fh:
        yield from csv.DictReader(line for line in fh if not line.startswith("#"))


@dataclass
class Bars:
    """Persistence intervals of one map (columns of bars.csv)."""
    dim: np.ndarray
    birth: np.ndarray
    death: np.ndarray        # inf for essential classes
    essential: np.ndarray    # bool
    wraps: np.ndarray        # bool

    def diagram(self, dim: int, essential: bool = False) -> np.ndarray:
        """(n, 2) array of (birth, death) of dimension `dim`: finite bars, or essential ones (death = inf)."""
        m = (self.dim == dim) & (self.essential == essential)
        return np.column_stack([self.birth[m], self.death[m]])


def read_bars(path: str | Path) -> dict[str, Bars]:
    cols: dict[str, dict[str, list]] = {}
    for r in _rows(path):
        c = cols.setdefault(r["file"], {k: [] for k in ("dim", "birth", "death", "essential", "wraps")})
        c["dim"].append(int(r["dim"]))
        c["birth"].append(float(r["birth"]))
        c["death"].append(float(r["death"]))
        c["essential"].append(r["essential"] == "1")
        c["wraps"].append(r["wraps"] == "1")
    return {f: Bars(**{k: np.array(v) for k, v in c.items()}) for f, c in cols.items()}


def read_betti(path: str | Path) -> dict[str, dict[str, np.ndarray]]:
    """{file: {"r", "beta0", "beta1", "wraps"}}, r decreasing."""
    out: dict[str, dict[str, list]] = {}
    for r in _rows(path):
        c = out.setdefault(r["file"], {"r": [], "beta0": [], "beta1": [], "wraps": []})
        c["r"].append(float(r["r"]))
        for k in ("beta0", "beta1", "wraps"):
            c[k].append(int(r[k]))
    return {f: {k: np.array(v) for k, v in c.items()} for f, c in out.items()}


def read_summary(path: str | Path) -> list[dict]:
    """One dict per map; numeric columns as float (nan where empty), file/format/error as strings."""
    rows = []
    for r in _rows(path):
        d = {}
        for k, v in r.items():
            if k in ("file", "format", "error"):
                d[k] = v
            else:
                d[k] = float(v) if v not in ("", None) else math.nan
        rows.append(d)
    return rows


# --------------------------------------------------------------------------------------- vectorization


def image_grid(b_min: float, b_max: float, p_min: float, p_max: float, n: int) -> tuple[np.ndarray, np.ndarray]:
    """Pixel centres of an n x n persistence image over birth in [b_min, b_max], persistence in [p_min, p_max]."""
    return np.linspace(b_min, b_max, n), np.linspace(p_min, p_max, n)


def persistence_image(diagram: np.ndarray, grid: tuple[np.ndarray, np.ndarray], sigma: float,
                      weight: str = "linear") -> np.ndarray:
    """Persistence image (Adams et al., JMLR 18(8), 2017) of finite (birth, death) pairs.

    In this superlevel filtration birth > death, so each pair maps to (birth, persistence = birth - death).
    Each point adds a Gaussian of width `sigma` weighted by its persistence (`linear`, the usual choice,
    which makes the image continuous in the diagram) or by 1 (`none`). The result is sampled at the pixel
    centres of `grid` (not integrated over pixels), shape (n_persistence, n_birth).
    """
    births, pers = grid
    if len(diagram) == 0:
        return np.zeros((len(pers), len(births)))
    b = diagram[:, 0]
    p = diagram[:, 0] - diagram[:, 1]
    w = p if weight == "linear" else np.ones_like(p)
    gb = np.exp(-((births[None, :] - b[:, None]) ** 2) / (2 * sigma**2))   # (points, n_birth)
    gp = np.exp(-((pers[None, :] - p[:, None]) ** 2) / (2 * sigma**2))     # (points, n_pers)
    return np.einsum("i,ip,ib->pb", w, gp, gb) / (2 * math.pi * sigma**2)



# ------------------------------------------------------------------------------------------ width field


def _distance_to(sources: np.ndarray, periodic: bool, metric: str) -> np.ndarray:
    """Distance from every site to the nearest True site of `sources` (inf if there is none)."""
    ny, nx = sources.shape
    if metric == "l1":  # breadth-first layers of the 4-neighbour lattice
        d = np.where(sources, 0.0, np.inf)
        front, step = sources.copy(), 0
        while front.any():
            step += 1
            grow = np.zeros_like(front)
            for axis, shift in ((0, 1), (0, -1), (1, 1), (1, -1)):
                rolled = np.roll(front, shift, axis=axis)
                if not periodic:  # drop the row or column that wrapped around
                    idx = [slice(None)] * 2
                    idx[axis] = 0 if shift == 1 else -1
                    rolled[tuple(idx)] = False
                grow |= rolled
            front = grow & np.isinf(d)
            d[front] = step
        return d
    sy, sx = np.nonzero(sources)
    yy, xx = np.mgrid[0:ny, 0:nx]
    d2 = np.full((ny, nx), np.inf)
    for start in range(0, len(sy), 256):  # brute force in chunks of sources
        dx = np.abs(xx[None] - sx[start:start + 256, None, None])
        dy = np.abs(yy[None] - sy[start:start + 256, None, None])
        if periodic:
            dx, dy = np.minimum(dx, nx - dx), np.minimum(dy, ny - dy)
        d2 = np.minimum(d2, (dx**2 + dy**2).min(axis=0))
    return np.sqrt(d2)


def width_field(obstacle: np.ndarray, periodic: bool = True, metric: str = "l1", signed: bool = False) -> np.ndarray:
    """Corridor width w on the free sites, as persistence2D computes it (without boundary_obstacle).

    On obstacles: 0, or minus the wall thickness (distance to the nearest free site) if `signed`.
    The L2 transform is brute force, meant for display of maps up to a few hundred sites a side.
    """
    w = _distance_to(obstacle, periodic, metric)
    w[np.isinf(w)] = 1.0
    if signed:
        d = _distance_to(~obstacle, periodic, metric)
        w = np.where(obstacle, -np.where(np.isinf(d), 1.0, d), w)
    else:
        w = np.where(obstacle, 0.0, w)
    return w
