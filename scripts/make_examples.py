#!/usr/bin/env python3
"""Regenerates the example maps of examples/ (deterministic; needs only NumPy).

  python3 scripts/make_examples.py            # from the repository root

  random_phi0.40.pgm   128 x 128, independent obstacles with probability 0.40, plain PGM, obstacles black
  blocks_q16.txt       64 x 64, 4 x 4 square obstacles every 16 sites (channels 12 wide), text format
  channels.map         48 x 48 chasers .map: two horizontal walls 3 sites thick, with one 2-site gap in one wall
  porous_phi0.40.pgm   128 x 128 porous medium: Gaussian random field (correlation length 4) thresholded at
                       obstacle fraction 0.40, periodic, plain PGM, obstacles black

examples/hyperuniform_phi0.40.pgm is not made here: it comes from hyperuniform2D
(https://github.com/manuxch/hyperuniform2D) with examples/hyperuniform_phi0.40.hyper2D.in.
"""
from pathlib import Path

import numpy as np

OUT = Path(__file__).resolve().parent.parent / "examples"


def write_pgm(path: Path, obstacle: np.ndarray, comment: str) -> None:
    """Plain PGM, maxval 1, obstacles black (0)."""
    ny, nx = obstacle.shape
    rows = "\n".join(" ".join("0" if o else "1" for o in row) for row in obstacle)
    path.write_text(f"P2\n# {comment}\n{nx} {ny}\n1\n{rows}\n")


def write_text(path: Path, obstacle: np.ndarray, comment: str) -> None:
    rows = "\n".join("".join("#" if o else "." for o in row) for row in obstacle)
    path.write_text(f"% {comment}\n{rows}\n")


def write_chasers_map(path: Path, obstacle: np.ndarray, kind: str) -> None:
    L = obstacle.shape[0]
    n = int(obstacle.sum())
    rows = "\n".join("".join("#" if o else "." for o in row) for row in obstacle)
    path.write_text(f"# chasers-map v1\nL: {L}\nkind: {kind}\nn_obs: {n}\nphi: {n / L**2}\n---\n{rows}\n")


def main() -> None:
    OUT.mkdir(exist_ok=True)
    rng = np.random.default_rng(20261110)
    random = rng.random((128, 128)) < 0.40
    write_pgm(OUT / "random_phi0.40.pgm", random, "random obstacles, p = 0.40, seed 20261110 (make_examples.py)")

    y, x = np.mgrid[0:64, 0:64]
    blocks = (x % 16 < 4) & (y % 16 < 4)
    write_text(OUT / "blocks_q16.txt", blocks, "4 x 4 blocks every 16 sites (make_examples.py)")

    walls = np.zeros((48, 48), dtype=bool)
    walls[10:13, :] = True
    walls[34:37, :] = True
    walls[34:37, 20:22] = False  # a 2-site gap in the lower wall
    write_chasers_map(OUT / "channels.map", walls, "periodic")

    # Periodic Gaussian random field: white noise filtered in Fourier space, thresholded at the 0.60 quantile.
    noise = np.random.default_rng(20261111).standard_normal((128, 128))
    k = np.fft.fftfreq(128)
    k2 = k[:, None] ** 2 + k[None, :] ** 2
    field = np.real(np.fft.ifft2(np.fft.fft2(noise) * np.exp(-2 * np.pi**2 * 4.0**2 * k2)))
    porous = field > np.quantile(field, 0.60)
    write_pgm(OUT / "porous_phi0.40.pgm", porous,
              "Gaussian random field, correlation length 4, thresholded at phi = 0.40, seed 20261111 (make_examples.py)")

    for p in sorted(OUT.glob("*")):
        print(p.relative_to(OUT.parent))


if __name__ == "__main__":
    main()
