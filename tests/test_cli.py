#!/usr/bin/env python3
"""End-to-end test of the persistence2D executable and of the Python scripts (run by ctest).

  python3 tests/test_cli.py build/persistence2D scripts
"""
import math
import subprocess
import sys
import tempfile
from pathlib import Path

binary = Path(sys.argv[1]).resolve()
scripts = Path(sys.argv[2]).resolve()
root = scripts.parent
examples = root / "examples"
sys.path.insert(0, str(scripts))
import p2d  # noqa: E402

failures = 0


def check(ok: bool, what: str) -> None:
    global failures
    print(("  ok   " if ok else "  FAIL ") + what)
    failures += not ok


def run(tmp: Path, params: str, *maps: Path) -> subprocess.CompletedProcess:
    (tmp / "p.in").write_text(params)
    return subprocess.run([str(binary), str(tmp / "p.in"), *map(str, maps)], cwd=tmp, capture_output=True, text=True)


with tempfile.TemporaryDirectory(prefix="p2d_cli_") as t:
    tmp = Path(t)
    print("command line")
    r = subprocess.run([str(binary), "--version"], capture_output=True, text=True)
    check(r.returncode == 0 and r.stdout.startswith("persistence2D "), "--version")
    r = subprocess.run([str(binary)], capture_output=True, text=True)
    check(r.returncode == 2 and "usage" in r.stderr, "no arguments: usage and exit code 2")
    r = run(tmp, "colour: dark\n", examples / "blocks_q16.txt")
    check(r.returncode == 2 and "unknown key" in r.stderr, "unknown key: exit code 2")

    print("free filtration on the example maps")
    maps = [examples / "blocks_q16.txt", examples / "channels.map", examples / "random_phi0.40.pgm"]
    r = run(tmp, "bars_file: out/bars.csv\nbetti_file: out/betti.csv\nsummary_file: out/summary.csv\n", *maps)
    check(r.returncode == 0, f"exit code 0 ({r.stderr.strip()})")
    s = {Path(row["file"]).name: row for row in p2d.read_summary(tmp / "out/summary.csv")}
    check(len(s) == 3 and all(not row["error"] for row in s.values()), "three summary rows, no errors")
    b = s["blocks_q16.txt"]
    check(b["beta0"] == 1 and b["beta1"] == 17 and b["wrap_dim"] == 2 and b["r_star"] == 6,
          "blocks: beta0 = 1, beta1 = 17, wraps twice, r* = 6 (channels 12 wide)")
    c = s["channels.map"]
    check(c["wrap_dim"] == 1 and c["r_wrap1"] == 11 and math.isnan(c["r_star"]),
          "channels: wraps once at r = 11, never in both directions")
    bars = p2d.read_bars(tmp / "out/bars.csv")
    betti = p2d.read_betti(tmp / "out/betti.csv")
    rb = bars[str(examples / "random_phi0.40.pgm")]
    rt = betti[str(examples / "random_phi0.40.pgm")]
    rr = s["random_phi0.40.pgm"]
    check(rt["beta0"][-1] == rr["beta0"] and rt["beta1"][-1] == rr["beta1"], "Betti curve ends at the summary values")
    check(int((rb.dim == 0).sum() - rb.essential[rb.dim == 0].sum()) == rr["h0_n"], "finite H0 bars = h0_n")

    print("signed filtration")
    r = run(tmp, "filtration: signed\nbars_file:\nbetti_file:\nsummary_file: signed.csv\n", examples / "channels.map")
    sc = p2d.read_summary(tmp / "signed.csv")[0]
    check(r.returncode == 0 and sc["r_star"] == -2 and sc["d_max"] == 2,
          "channels: crossing a 3-site wall gives r* = -2 (d_max = 2)")
    check(sc["r_wrap1"] == c["r_wrap1"] and sc["beta1"] == c["beta1"], "signed agrees with free for r >= 1")

    print("obstacle colour")
    r = run(tmp, "obstacle: light\nbars_file:\nbetti_file:\nsummary_file: hu.csv\n", examples / "hyperuniform_phi0.40.pgm")
    hu = p2d.read_summary(tmp / "hu.csv")[0]
    check(r.returncode == 0 and abs(hu["phi"] - 0.4) < 1e-3, "hyperuniform example read with obstacle: light (phi = 0.4)")

    print("errors")
    r = run(tmp, "bars_file:\nbetti_file:\nsummary_file: err.csv\n", examples / "blocks_q16.txt", tmp / "missing.pgm")
    rows = p2d.read_summary(tmp / "err.csv")
    check(r.returncode == 1 and len(rows) == 2 and rows[1]["error"] and not rows[0]["error"],
          "a missing map: exit code 1, error in its summary row, the other map analysed")

    print("Python scripts")
    for path, obstacle in ((examples / "random_phi0.40.pgm", "dark"), (examples / "hyperuniform_phi0.40.pgm", "light")):
        g = p2d.read_map(path, obstacle)
        check(g.shape == (128, 128) and abs(g.mean() - 0.4) < 0.01, f"read_map {path.name}")
    porous = p2d.read_map(examples / "porous_phi0.40.pgm")
    w = p2d.width_field(porous, True, "l2", signed=True)
    r = run(tmp, "metric: l2\nfiltration: signed\nbars_file:\nbetti_file:\nsummary_file: po.csv\n", examples / "porous_phi0.40.pgm")
    po = p2d.read_summary(tmp / "po.csv")[0]
    check(abs(w.max() - po["w_max"]) < 1e-9 and abs(-w.min() - po["d_max"]) < 1e-9,
          "width_field (Python) matches w_max and d_max of the executable")
    img = p2d.persistence_image(rb.diagram(0), p2d.image_grid(1, 3, 0, 2, 8), sigma=0.3)
    check(img.shape == (8, 8) and img.max() > 0, "persistence_image")
    try:
        import matplotlib  # noqa: F401
        out = tmp / "fig.png"
        r = subprocess.run([sys.executable, str(scripts / "plot_persistence.py"), str(examples / "blocks_q16.txt"),
                            "--binary", str(binary), "--out", str(out)], capture_output=True, text=True)
        check(r.returncode == 0 and out.stat().st_size > 10000, f"plot_persistence.py ({r.stderr.strip()[-200:]})")
        out = tmp / "cmp.png"
        r = subprocess.run([sys.executable, str(scripts / "compare_maps.py"), f"{examples / 'blocks_q16.txt'}::blocks",
                            f"{examples / 'hyperuniform_phi0.40.pgm'}:light:hu", "--binary", str(binary),
                            "--out", str(out)], capture_output=True, text=True)
        check(r.returncode == 0 and out.stat().st_size > 10000, f"compare_maps.py ({r.stderr.strip()[-200:]})")
    except ImportError:
        print("  skip plot_persistence.py and compare_maps.py (no Matplotlib)")

    print("continuum filtrations on an .obs file")
    obs = tmp / "pair.obs"
    obs.write_text("# format: 1\n# Lx: 10\n# Ly: 6\n# R: 0.2\n# boundary: periodic_x\n# n_obstacles: 2\n4.5 3\n5.5 3\n")
    r = run(tmp, "filtration: height\npixel: 0.05\nbars_file: h/bars.csv\nbetti_file:\nsummary_file: h/summary.csv\n", obs)
    check(r.returncode == 0, f"height on a .obs file ({r.stderr.strip()})")
    import csv
    rows = list(csv.DictReader(l for l in (tmp / "h/bars.csv").read_text().splitlines() if not l.startswith("#")))
    trap = [x for x in rows if x["dim"] == "0" and x["essential"] == "0"]
    check(len(trap) == 1 and abs(float(trap[0]["persistence"]) - 0.21) < 0.06 and "xs" in rows[0],
          "one V trap with barrier ~ 0.21 and positions in bars.csv")
    r = run(tmp, "filtration: height_exact\nbars_file: e/bars.csv\nbetti_file:\nsummary_file:\n", obs)
    rows = list(csv.DictReader(l for l in (tmp / "e/bars.csv").read_text().splitlines() if not l.startswith("#")))
    check(r.returncode == 0 and len(rows) == 1 and abs(float(rows[0]["persistence"]) - (0.7 - 0.7 ** 2 * 0 - (0.49 - 0.25) ** 0.5)) < 1e-9,
          "height_exact: one trap, exact barrier 0.7 - sqrt(0.49 - 0.25)")
    r = run(tmp, "filtration: height\n", examples / "blocks_q16.txt")
    check(r.returncode == 1 and "needs a .obs" in r.stderr, "height on a text map is an error")
    r = run(tmp, "filtration: free\n", obs)
    check(r.returncode == 1 and "needs filtration height" in r.stderr, "an .obs file with filtration free is an error")
    r = run(tmp, "filtration: rips\nmax_radius: 1\nbars_file: r/bars.csv\nbetti_file:\nsummary_file: r/summary.csv\n", obs)
    rows = list(csv.DictReader(l for l in (tmp / "r/bars.csv").read_text().splitlines() if not l.startswith("#")))
    fin = [x for x in rows if x["essential"] == "0"]
    check(r.returncode == 0 and len(fin) == 1 and abs(float(fin[0]["persistence"]) - 0.5) < 1e-12 and abs(float(fin[0]["xs"]) - 5.0) < 1e-12,
          f"rips: one H0 bar of persistence 1/2 at the midpoint ({r.stderr.strip()})")
    summ = list(csv.DictReader(l for l in (tmp / "r/summary.csv").read_text().splitlines() if not l.startswith("#")))
    check(summ[0]["format"] == "obs-rips" and summ[0]["h0_n"] == "1" and summ[0]["h0_essential"] == "1", "rips: summary row (format obs-rips, 1 finite + 1 essential H0)")
    r = run(tmp, "filtration: rips\nmax_radius: 3\n", obs)
    check(r.returncode == 1 and "too large" in r.stderr, "rips: max_radius too large for the box is an error")
    r = run(tmp, "filtration: rips\n", examples / "blocks_q16.txt")
    check(r.returncode == 1 and "needs a .obs" in r.stderr, "rips on a text map is an error")


print("all checks passed" if not failures else f"{failures} check(s) failed")
sys.exit(1 if failures else 0)
