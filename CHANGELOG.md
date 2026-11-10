# Changelog

All notable changes to this project are documented here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/) and the project uses
[semantic versioning](https://semver.org/).

## [0.1.0] - 2026-11-10

First version.

### Added

- `persistence2D`: persistent homology of the free space of 2D binary maps, filtered by corridor width
  (distance to the nearest obstacle), on GUDHI periodic bitmap cubical complexes with Z2 coefficients.
- Filtrations `free` (free space only) and `signed` (continued into the obstacles by wall thickness).
- L1 (breadth-first search) and exact Euclidean (Felzenszwalb–Huttenlocher) distance transforms, periodic
  (minimum image) or open, optionally with the frame as an obstacle.
- 4-connectivity (vertex construction) and 8-connectivity (top-cell construction).
- Wrapping classes on the torus: wrap dimension of every sublevel and the percolating widths `r_wrap1`
  and `r_star`.
- Inputs: PBM (P1, P4), PGM (P2, P5, 8 and 16 bit) with dark or light obstacles, chasers `.map`, text.
- Outputs: `bars.csv`, `betti.csv`, `summary.csv` (Betti numbers, wrap dimension, bar counts, persistence
  quantiles, total persistence and persistence entropy per dimension); per-map errors do not stop a run.
- Maps processed in parallel; CMake options `P2D_NATIVE` and `P2D_STATIC` (static binary for clusters).
- Python: `scripts/p2d.py` (map and CSV readers, width field, persistence images),
  `scripts/plot_persistence.py` (four-panel figure of one map), `scripts/compare_maps.py` (maps side by
  side: map, barcode, persistence diagram, Betti curves), `scripts/make_examples.py`.
- Examples (blocks, channels, random, porous and hyperuniform maps) and their parameter files; analysis
  example in the README comparing three maps at the same obstacle fraction.
- Tests: hand-computable topology, Euler characteristic at every level of random maps, distance
  transforms against brute force, readers, and an end-to-end test of the executable and scripts.

[0.1.0]: https://github.com/manuxch/persistence2D/releases/tag/v0.1.0
