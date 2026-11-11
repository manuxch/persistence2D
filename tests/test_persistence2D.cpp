// Unit tests for persistence2D. No framework: each check prints its result and the exit code counts failures
// (run with ctest, or ./build/test_persistence2D).
#include "continuum.hpp"
#include "descriptors.hpp"
#include "distance.hpp"
#include "grid.hpp"
#include "output.hpp"
#include "params.hpp"
#include "persistence.hpp"
#include "rips.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <format>
#include <fstream>
#include <iostream>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

#include <unistd.h>

using namespace p2d;

namespace {

int failures = 0;

void check(bool ok, const std::string& what) {
    std::cout << (ok ? "  ok   " : "  FAIL ") << what << '\n';
    if (!ok) ++failures;
}

Grid empty_grid(std::size_t nx, std::size_t ny) {
    Grid g;
    g.nx = nx;
    g.ny = ny;
    g.obstacle.assign(nx * ny, 0);
    return g;
}

Grid random_grid(std::size_t nx, std::size_t ny, double phi, unsigned seed) {
    Grid g = empty_grid(nx, ny);
    std::mt19937 rng(seed);
    std::bernoulli_distribution occ(phi);
    for (auto& c : g.obstacle) c = occ(rng) ? 1 : 0;
    return g;
}

// q x q cells, each with a b x b obstacle block in its corner: free channels of width q - b.
Grid blocks_grid(std::size_t L, std::size_t q, std::size_t b) {
    Grid g = empty_grid(L, L);
    for (std::size_t y = 0; y < L; ++y)
        for (std::size_t x = 0; x < L; ++x)
            if (x % q < b && y % q < b) g.obstacle[g.index(x, y)] = 1;
    return g;
}

PersistenceOptions opts(Metric m = Metric::l1, bool periodic = true, int conn = 4, bool boundary = false,
                        Filtration filt = Filtration::free) {
    PersistenceOptions o;
    o.filtration = filt;
    o.distance.metric = m;
    o.distance.periodic = periodic;
    o.distance.boundary_obstacle = boundary;
    o.connectivity = conn;
    return o;
}

// Euler characteristic of the cubical complex of the sites with w >= r, counted cell by cell.
long euler_direct(const Grid& g, const PersistenceResult& res, double r, bool periodic, int conn) {
    const auto& w = res.w;
    const long nx = static_cast<long>(g.nx), ny = static_cast<long>(g.ny);
    auto in = [&](long x, long y) {
        if (periodic) { x = (x % nx + nx) % nx; y = (y % ny + ny) % ny; }
        else if (x < 0 || y < 0 || x >= nx || y >= ny) return false;
        const auto i = static_cast<std::size_t>(x + nx * y);
        return res.in_filtration(g, i) && w[i] >= r;
    };
    long v = 0, e = 0, f = 0;
    if (conn == 4) {  // sites are vertices; edges join 4-neighbours; squares need their four corners
        const long ex = periodic ? nx : nx - 1, ey = periodic ? ny : ny - 1;
        for (long y = 0; y < ny; ++y)
            for (long x = 0; x < nx; ++x) v += in(x, y);
        for (long y = 0; y < ny; ++y)
            for (long x = 0; x < ex; ++x) e += in(x, y) && in(x + 1, y);
        for (long y = 0; y < ey; ++y)
            for (long x = 0; x < nx; ++x) e += in(x, y) && in(x, y + 1);
        for (long y = 0; y < ey; ++y)
            for (long x = 0; x < ex; ++x) f += in(x, y) && in(x + 1, y) && in(x, y + 1) && in(x + 1, y + 1);
    } else {  // sites are pixels; a corner or an edge is present if some pixel containing it is
        const long cx = periodic ? nx : nx + 1, cy = periodic ? ny : ny + 1;
        for (long y = 0; y < ny; ++y)
            for (long x = 0; x < nx; ++x) f += in(x, y);
        for (long y = 0; y < cy; ++y)  // corner (x, y) touches pixels (x-1..x, y-1..y)
            for (long x = 0; x < cx; ++x) v += in(x - 1, y - 1) || in(x, y - 1) || in(x - 1, y) || in(x, y);
        for (long y = 0; y < cy; ++y)  // horizontal edge from corner (x, y): pixels (x, y-1) and (x, y)
            for (long x = 0; x < nx; ++x) e += in(x, y - 1) || in(x, y);
        for (long y = 0; y < ny; ++y)  // vertical edge from corner (x, y): pixels (x-1, y) and (x, y)
            for (long x = 0; x < cx; ++x) e += in(x - 1, y) || in(x, y);
    }
    return v - e + f;
}

long euler_from_bars(const std::vector<Bar>& bars, double r) {
    return betti_at(bars, 0, r) - betti_at(bars, 1, r) + betti_at(bars, 2, r);
}

void test_empty_torus() {
    std::cout << "empty torus\n";
    const Grid g = empty_grid(8, 8);
    const auto res = compute_persistence(g, opts());
    check(betti_at(res.bars, 0, 1) == 1 && betti_at(res.bars, 1, 1) == 2 && betti_at(res.bars, 2, 1) == 1,
          "8 x 8 torus without obstacles: beta = (1, 2, 1)");
    check(wraps_at(res.bars, 1) == 2, "both H1 classes wrap");
    const auto open = compute_persistence(g, opts(Metric::l1, false));
    check(betti_at(open.bars, 0, 1) == 1 && betti_at(open.bars, 1, 1) == 0 && betti_at(open.bars, 2, 1) == 0,
          "open square without obstacles: beta = (1, 0, 0)");
}

void test_single_obstacle() {
    std::cout << "single obstacle on a 10 x 10 torus\n";
    Grid g = empty_grid(10, 10);
    g.obstacle[0] = 1;
    const auto res = compute_persistence(g, opts());
    const Summary s = summarize(g, res);
    check(s.w_max == 10.0, std::format("w_max = L = 10 (got {})", s.w_max));
    check(s.beta0 == 1 && s.beta1 == 2 && s.wrap_dim == 2, "punctured torus: beta0 = 1, beta1 = 2, both wrap");
    check(s.r_star == 5.0 && s.r_wrap1 == 5.0, std::format("r* = r_wrap1 = L/2 = 5 (got {}, {})", s.r_star, s.r_wrap1));
    check(s.h0.n_finite == 0 && s.h0.n_essential == 1, "one H0 class, born at the farthest site");
}

void test_blocks() {
    std::cout << "blocks: L = 24, q = 6, block 2 x 2 (16 blocks, channels of width 4)\n";
    const Grid g = blocks_grid(24, 6, 2);
    const auto res = compute_persistence(g, opts());
    const Summary s = summarize(g, res);
    check(s.beta0 == 1, "free space connected");
    check(s.beta1 == 16 + 1, std::format("beta1 = blocks + 1 = 17 on the torus (got {})", s.beta1));
    check(s.wrap_dim == 2, "free space wraps in both directions");
    check(s.r_star == 2.0, std::format("r* = floor((channel + 1) / 2) = 2 (got {})", s.r_star));
    const auto open = compute_persistence(g, opts(Metric::l1, false));
    const Summary so = summarize(g, open);
    check(so.wrap_dim == 0 && std::isnan(so.r_star), "open boundary: nothing wraps, r* undefined");
    check(so.beta1 != s.beta1, std::format("periodic and open complexes differ (beta1 {} vs {})", s.beta1, so.beta1));
}

void test_connectivity() {
    std::cout << "connectivity on a 4 x 4 checkerboard\n";
    Grid g = empty_grid(4, 4);
    for (std::size_t y = 0; y < 4; ++y)
        for (std::size_t x = 0; x < 4; ++x) g.obstacle[g.index(x, y)] = (x + y) % 2;
    const auto r4 = compute_persistence(g, opts(Metric::l1, true, 4));
    const auto r8 = compute_persistence(g, opts(Metric::l1, true, 8));
    check(betti_at(r4.bars, 0, 1) == 8, "4-connectivity: 8 isolated free sites");
    check(betti_at(r8.bars, 0, 1) == 1, "8-connectivity: diagonal contacts join them into one component");
}

void test_euler_random() {
    std::cout << "Euler characteristic at every level, random maps\n";
    int bad = 0, cases = 0;
    for (unsigned seed = 1; seed <= 6; ++seed)
        for (double phi : {0.2, 0.45, 0.6})
            for (Metric m : {Metric::l1, Metric::l2})
                for (bool periodic : {true, false})
                    for (int conn : {4, 8})
                        for (Filtration filt : {Filtration::free, Filtration::signed_distance}) {
                            const Grid g = random_grid(17 + seed, 13 + 2 * seed, phi, seed * 97 + 1);
                            const auto res = compute_persistence(g, opts(m, periodic, conn, false, filt));
                            for (double r : levels(g, res)) {
                                ++cases;
                                if (euler_direct(g, res, r, periodic, conn) != euler_from_bars(res.bars, r)) ++bad;
                            }
                        }
    check(bad == 0, std::format("chi(V_r) from the cells equals beta0 - beta1 + beta2 ({} levels, {} mismatches)",
                                cases, bad));
}

void test_signed() {
    std::cout << "signed filtration\n";
    // It continues the free one: the same Betti numbers and wrap dimension at every free width.
    int bad = 0;
    for (unsigned seed = 1; seed <= 6; ++seed)
        for (double phi : {0.3, 0.5})
            for (int conn : {4, 8}) {
                const Grid g = random_grid(20, 16, phi, seed);
                const auto rf = compute_persistence(g, opts(Metric::l1, true, conn));
                const auto rs = compute_persistence(g, opts(Metric::l1, true, conn, false, Filtration::signed_distance));
                for (double r : levels(g, rf))
                    for (int d : {0, 1})
                        if (betti_at(rf.bars, d, r) != betti_at(rs.bars, d, r) || wraps_at(rf.bars, r) != wraps_at(rs.bars, r))
                            ++bad;
            }
    check(bad == 0, std::format("signed and free agree for r >= 1 ({} mismatches)", bad));

    // A vertical wall three sites thick on a 12 x 12 torus: the free space wraps only vertically; crossing the
    // wall needs its middle column (L1 thickness d = 2), so r* = -2.
    Grid g = empty_grid(12, 12);
    for (std::size_t y = 0; y < 12; ++y)
        for (std::size_t x = 0; x < 3; ++x) g.obstacle[g.index(x, y)] = 1;
    const auto rs = compute_persistence(g, opts(Metric::l1, true, 4, false, Filtration::signed_distance));
    const Summary s = summarize(g, rs);
    check(s.wrap_dim == 1 && s.beta0 == 1 && s.beta1 == 1, "free band: one wrapping loop");
    check(s.r_wrap1 == 5.0, std::format("r_wrap1 = widest free row = 5 (got {})", s.r_wrap1));
    check(s.r_star == -2.0 && s.d_max == 2.0, std::format("r* = -2, d_max = 2 (got {}, {})", s.r_star, s.d_max));
    check(betti_at(rs.bars, 2, -2) == 1, "the whole torus is reached: beta2 = 1");
}

void test_distances() {
    std::cout << "distance transforms against brute force\n";
    int bad = 0;
    for (unsigned seed = 1; seed <= 5; ++seed)
        for (bool periodic : {true, false})
            for (bool boundary : {false, true}) {
                if (periodic && boundary) continue;
                const Grid g = random_grid(15 + seed, 11 + seed, 0.08, seed);
                for (Metric m : {Metric::l1, Metric::l2}) {
                    DistanceOptions o{m, periodic, boundary};
                    const auto w = distance_to_obstacles(g, o);
                    const long nx = static_cast<long>(g.nx), ny = static_cast<long>(g.ny);
                    for (long y = 0; y < ny; ++y)
                        for (long x = 0; x < nx; ++x) {
                            double best = std::numeric_limits<double>::infinity();
                            auto consider = [&](double dx, double dy) {
                                const double d = m == Metric::l1 ? std::abs(dx) + std::abs(dy) : std::hypot(dx, dy);
                                best = std::min(best, d);
                            };
                            for (long oy = 0; oy < ny; ++oy)
                                for (long ox = 0; ox < nx; ++ox) {
                                    if (!g.obstacle[g.index(ox, oy)]) continue;
                                    double dx = std::abs(x - ox), dy = std::abs(y - oy);
                                    if (periodic) { dx = std::min(dx, nx - dx); dy = std::min(dy, ny - dy); }
                                    consider(dx, dy);
                                }
                            if (boundary) {  // obstacles at x = -1, x = nx, y = -1, y = ny
                                for (long ox = -1; ox <= nx; ++ox) { consider(x - ox, y + 1); consider(x - ox, ny - y); }
                                for (long oy = -1; oy <= ny; ++oy) { consider(x + 1, y - oy); consider(nx - x, y - oy); }
                            }
                            if (std::abs(w[g.index(x, y)] - best) > 1e-9) ++bad;
                        }
                }
            }
    check(bad == 0, std::format("L1 and L2, periodic, open and boundary_obstacle ({} mismatches)", bad));
}

void test_readers() {
    std::cout << "readers: the same map in every format\n";
    const auto dir = std::filesystem::temp_directory_path() / std::format("p2d_test_{}", ::getpid());
    std::filesystem::create_directories(dir);
    const Grid g = random_grid(13, 9, 0.4, 7);
    const auto path = [&](const char* name) { return (dir / name).string(); };

    write_pgm(path("dark.pgm"), g, ObstacleColor::dark);
    write_pgm(path("light.pgm"), g, ObstacleColor::light);
    write_text(path("map.txt"), g);
    {
        std::ofstream out(path("map.map"));
        out << "# chasers-map v1\nL: 13\nkind: random  # comment\n---\n";  // not square: L check must fail
    }
    {
        std::ofstream out(path("raw.pgm"), std::ios::binary);  // P5, maxval 255, obstacles dark (0)
        out << "P5\n# comment\n13 9\n255\n";
        for (auto c : g.obstacle) out.put(c ? char(0) : char(255));
    }
    {
        std::ofstream out(path("raw.pbm"), std::ios::binary);  // P4, bit 1 = black = obstacle
        out << "P4\n13 9\n";
        for (std::size_t y = 0; y < g.ny; ++y) {
            unsigned char byte = 0;
            for (std::size_t x = 0; x < g.nx; ++x) {
                if (g.obstacle[g.index(x, y)]) byte |= static_cast<unsigned char>(1u << (7 - x % 8));
                if (x % 8 == 7 || x + 1 == g.nx) { out.put(static_cast<char>(byte)); byte = 0; }
            }
        }
    }
    {
        std::ofstream out(path("plain.pbm"));
        out << "P1\n13 9\n";
        for (std::size_t i = 0; i < g.size(); ++i) out << int(g.obstacle[i]);  // packed digits are legal
        out << '\n';
    }
    auto same = [&](const Grid& h) { return h.nx == g.nx && h.ny == g.ny && h.obstacle == g.obstacle; };
    check(same(read_grid(path("dark.pgm"), ObstacleColor::dark)), "P2, obstacles dark");
    check(same(read_grid(path("light.pgm"), ObstacleColor::light)), "P2, obstacles light");
    check(same(read_grid(path("raw.pgm"), ObstacleColor::dark)), "P5");
    check(same(read_grid(path("raw.pbm"), ObstacleColor::dark)), "P4");
    check(same(read_grid(path("plain.pbm"), ObstacleColor::dark)), "P1 with packed digits");
    check(same(read_grid(path("map.txt"))), "text");
    check(detect_format(path("map.map")) == FileFormat::chasers_map, "chasers .map detected");

    // A square chasers map, compared with its PGM version through the whole pipeline.
    const Grid sq = random_grid(16, 16, 0.35, 11);
    {
        std::ofstream out(path("sq.map"));
        out << std::format("# chasers-map v1\nL: 16\nkind: random\nn_obs: {}\n---\n", sq.n_obstacles());
        for (std::size_t y = 0; y < 16; ++y) {
            for (std::size_t x = 0; x < 16; ++x) out << (sq.obstacle[sq.index(x, y)] ? '#' : '.');
            out << '\n';
        }
    }
    write_pgm(path("sq.pgm"), sq, ObstacleColor::dark);
    const Grid a = read_grid(path("sq.map")), b = read_grid(path("sq.pgm"));
    check(a.header.at("kind") == "random", "header of the .map is parsed");
    const auto ra = compute_persistence(a, opts()), rb = compute_persistence(b, opts());
    bool equal = ra.bars.size() == rb.bars.size();
    for (std::size_t i = 0; equal && i < ra.bars.size(); ++i)
        equal = ra.bars[i].dim == rb.bars[i].dim && ra.bars[i].birth == rb.bars[i].birth &&
                ra.bars[i].death == rb.bars[i].death && ra.bars[i].essential == rb.bars[i].essential;
    check(equal, ".map and PGM of the same map give identical diagrams");

    bool threw = false;
    try { read_grid(path("map.map")); } catch (const std::runtime_error&) { threw = true; }
    check(threw, "a .map whose size disagrees with L is rejected");
    std::filesystem::remove_all(dir);
}

void test_quantile_and_params() {
    std::cout << "quantile and parameters\n";
    check(std::abs(quantile({4, 1, 3, 2}, 0.9) - 3.7) < 1e-12, "quantile matches numpy (linear)");
    check(std::abs(quantile({1, 2, 3, 4}, 0.5) - 2.5) < 1e-12, "median of an even sample");
    const Params p = parse_params(
        "metric: l2  # comment\nperiodic: false\nconnectivity: 8\ninput: a.pgm b.map\nfiltration: signed\n");
    check(p.persistence.distance.metric == Metric::l2 && !p.persistence.distance.periodic &&
              p.persistence.filtration == Filtration::signed_distance &&
              p.persistence.connectivity == 8 && p.inputs.size() == 2,
          "keys parsed");
    auto throws = [](const std::string& text) {
        try { parse_params(text); } catch (const std::runtime_error&) { return true; }
        return false;
    };
    check(throws("colour: dark\n"), "unknown key rejected");
    check(throws("connectivity: 6\n"), "invalid connectivity rejected");
    check(throws("metric: l1\nmetric: l2\n"), "duplicate key rejected");
    check(throws("boundary_obstacle: true\n"), "boundary_obstacle on a periodic lattice rejected");
}

}  // namespace


// --- Continuum filtrations (.obs input) -------------------------------------------------------------------

ObsFile make_obs(double Lx, double Ly, double R, std::vector<std::array<double, 2>> c, bool periodic = true) {
    ObsFile o;
    o.Lx = Lx;
    o.Ly = Ly;
    o.R = R;
    o.periodic_x = periodic;
    o.centres = std::move(c);
    return o;
}

void test_continuum() {
    std::cout << "continuum\n";
    // Reader: header, count check, errors.
    const auto dir = std::filesystem::temp_directory_path() / std::format("p2d_obs_{}", getpid());
    std::filesystem::create_directories(dir);
    auto write = [&](const std::string& name, const std::string& body) {
        std::ofstream(dir / name) << body;
        return (dir / name).string();
    };
    const std::string hdr = "# format: 1\n# Lx: 10\n# Ly: 6\n# R: 0.2\n# boundary: periodic_x\n";
    const auto good = write("good.obs", hdr + "# n_obstacles: 2\n# columns: x y\n4.5 3\n5.5 3\n");
    check(is_obs_file(good), "is_obs_file recognises '# format: 1'");
    const ObsFile o = read_obs(good);
    check(o.centres.size() == 2 && o.Lx == 10 && o.R == 0.2 && o.periodic_x, "read_obs: header and centres");
    check(!is_obs_file(write("map.txt", "..#\n")), "is_obs_file rejects a text map");
    auto throws = [&](const std::string& path) {
        try { read_obs(path); } catch (const std::runtime_error&) { return true; }
        return false;
    };
    check(throws(write("count.obs", hdr + "# n_obstacles: 3\n4.5 3\n")), "read_obs: wrong obstacle count");
    check(throws(write("nohdr.obs", "# format: 1\n# Lx: 10\n1 1\n")), "read_obs: missing header keys");
    check(throws(write("bnd.obs", "# format: 1\n# Lx: 10\n# Ly: 6\n# R: 0.2\n# boundary: periodic_xy\n1 1\n")),
          "read_obs: unsupported boundary");
    std::filesystem::remove_all(dir);

    // Distance to the nearest centre against brute force, with minimum image in x.
    std::mt19937 rng(7);
    std::uniform_real_distribution<double> ux(0, 12), uy(0, 5);
    std::vector<std::array<double, 2>> pts;
    for (int i = 0; i < 40; ++i) pts.push_back({ux(rng), uy(rng)});
    for (bool per : {true, false}) {
        const ObsFile r = make_obs(12, 5, 0.2, pts, per);
        const std::size_t nx = 60, ny = 25;
        const auto d = centre_distance(r, nx, ny);
        double worst = 0;
        for (std::size_t j = 0; j < ny; ++j)
            for (std::size_t i = 0; i < nx; ++i) {
                const double x = (i + 0.5) * 12.0 / nx, y = (j + 0.5) * 5.0 / ny;
                double best = 1e9;
                for (const auto& p : pts) {
                    double dx = std::abs(p[0] - x);
                    if (per) dx = std::min(dx, 12 - dx);
                    best = std::min(best, std::hypot(dx, p[1] - y));
                }
                worst = std::max(worst, std::abs(best - d[i + nx * j]));
            }
        check(worst < 1e-12, std::format("centre_distance matches brute force (periodic={})", per));
    }

    // One V trap with a known barrier: two obstacles 2*delta apart side by side at height y0. A grain (rho = R + 1/2)
    // rests at the cusp, y0 + sqrt(rho^2 - delta^2), and escapes over the top of the discs, y0 + rho.
    {
        const double R = 0.2, rho = R + 0.5, delta = 0.5, y0 = 3.0;
        const ObsFile v = make_obs(10, 6, R, {{5 - delta, y0}, {5 + delta, y0}});
        ContinuumOptions co;
        co.pixel = 0.02;
        const auto res = compute_continuum(v, co);
        const Bar* trap = nullptr;
        int finite0 = 0, essential0 = 0;
        for (const auto& b : res.bars) {
            if (b.dim != 0) continue;
            (b.essential ? essential0 : finite0)++;
            if (!b.essential) trap = &b;
        }
        const double cusp = y0 + std::sqrt(rho * rho - delta * delta);
        check(finite0 == 1 && essential0 == 1, "height: one V trap and the global minimum");
        if (trap) {
            check(std::abs(trap->birth - cusp) < 0.03 && std::abs(trap->x - 5.0) < 0.03,
                  std::format("height: trap at the cusp ({:.3f}, {:.3f}) vs ({:.3f}, 5)", trap->y, trap->x, cusp));
            check(std::abs(trap->death - (y0 + rho)) < 0.03 && std::abs(trap->persistence() - (rho - std::sqrt(rho * rho - delta * delta))) < 0.04,
                  std::format("height: barrier {:.3f} vs {:.3f}", trap->persistence(), rho - std::sqrt(rho * rho - delta * delta)));
            check(std::abs(trap->ys - (y0 + rho)) < 0.03, "height: the saddle is at the height of the obstacle tops");
        }
        // The pair far enough apart is no trap: only the global minimum remains.
        const ObsFile w = make_obs(10, 6, R, {{5 - 0.8, y0}, {5 + 0.8, y0}});
        int f = 0;
        for (const auto& b : compute_continuum(w, co).bars)
            if (b.dim == 0 && !b.essential) ++f;
        check(f == 0, "height: obstacles farther apart than 2 rho make no trap");
    }

    // Distance filtration of a triangular lattice: every triangle is a chamber of the free space, born at the
    // circumradius a/sqrt(3) and dying at the bottleneck a/2.
    {
        const double a = 4.0;
        std::vector<std::array<double, 2>> c;
        const int ncol = 6, nrow = 8;
        for (int j = 0; j < nrow; ++j)
            for (int i = 0; i < ncol; ++i) c.push_back({a * (i + 0.5 * (j % 2)), a * std::sqrt(3.0) / 2 * j});
        const ObsFile lat = make_obs(ncol * a, nrow * a * std::sqrt(3.0) / 2, 0.2, c);
        ContinuumOptions co;
        co.mode = ContinuumMode::centres;
        co.pixel = 0.02;
        co.min_persistence = 0.04;   // two pixels: the distance field has discretisation noise
        const auto res = compute_continuum(lat, co);
        int inner = 0;
        double worst = 0;
        for (const auto& b : res.bars) {
            if (b.dim != 0 || b.essential || b.y < a || b.y > lat.Ly - a) continue;
            ++inner;
            worst = std::max({worst, std::abs(b.birth - a / std::sqrt(3.0)), std::abs(b.death - a / 2)});
        }
        // The band between y = a and Ly - a holds 2 triangles per obstacle of the rows in it.
        check(inner >= 2 * ncol * 2 && worst < 0.04, std::format("centres: {} chambers of the lattice, worst error {:.3f}", inner, worst));
    }
    // height_exact: the same pair is a trap with the barrier given by the formula, with no pixel error.
    {
        const double R = 0.2, rho = R + 0.5, delta = 0.5, y0 = 3.0;
        const ObsFile v = make_obs(10, 6, R, {{5 - delta, y0}, {5 + delta, y0}});
        const auto bars = height_exact(v, rho);
        const double cusp = y0 + std::sqrt(rho * rho - delta * delta);
        check(bars.size() == 1 && std::abs(bars[0].birth - cusp) < 1e-12 && std::abs(bars[0].death - (y0 + rho)) < 1e-12,
              "height_exact: a pair is one trap, cusp and barrier exact");
        const ObsFile w = make_obs(10, 6, R, {{5 - 0.8, y0}, {5 + 0.8, y0}});
        check(height_exact(w, rho).empty(), "height_exact: obstacles farther apart than 2 rho make no trap");
        // Three stacked pairs: the lower trap is covered by nothing, the upper pair sits above a disc top.
        const ObsFile st = make_obs(10, 12, R, {{4.5, 3}, {5.5, 3}, {4.5, 7}, {5.5, 7}});
        const auto b2 = height_exact(st, rho);
        int fin = 0;
        for (const auto& b : b2) fin += !b.essential;
        check(b2.size() == 2 && fin == 2, "height_exact: two separate pairs, two traps that both drain");
    }
    // Random arrangement: exact barriers against the pixel grid (grid error is about one pixel).
    {
        std::mt19937 g(11);
        std::uniform_real_distribution<double> ux(0, 8), uy(2.5, 6.5);
        std::vector<std::array<double, 2>> c;
        for (int i = 0; i < 28; ++i) c.push_back({ux(g), uy(g)});
        const ObsFile r = make_obs(8, 9, 0.2, c);
        const auto exact = height_exact(r, 0.7);
        ContinuumOptions co;
        co.pixel = 0.01;
        const auto grid = compute_continuum(r, co).bars;
        int checked = 0, bad = 0;
        for (const auto& e : exact) {
            if (e.essential || e.persistence() < 0.12) continue;
            ++checked;
            bool ok = false;
            for (const auto& b : grid) {
                if (b.dim != 0 || b.essential) continue;
                double dx = std::abs(b.x - e.x);
                dx = std::min(dx, 8 - dx);
                if (std::hypot(dx, b.y - e.y) < 0.05 && std::abs(b.persistence() - e.persistence()) < 0.05) ok = true;
            }
            bad += !ok;
        }
        check(checked > 3 && bad == 0, std::format("height_exact agrees with the fine grid on {} traps ({} mismatches)", checked, bad));
    }
}

void test_rips() {
    std::cout << "rips\n";
    auto make = [](double Lx, double Ly, bool periodic, std::vector<std::array<double, 2>> c) {
        ObsFile o;
        o.Lx = Lx;
        o.Ly = Ly;
        o.R = 0.2;
        o.periodic_x = periodic;
        o.centres = std::move(c);
        return o;
    };
    auto finite = [](const std::vector<Bar>& b, int dim) {
        std::vector<Bar> v;
        for (const auto& x : b) if (x.dim == dim && !x.essential) v.push_back(x);
        return v;
    };
    auto n_ess = [](const std::vector<Bar>& b, int dim) {
        int n = 0;
        for (const auto& x : b) n += x.dim == dim && x.essential;
        return n;
    };

    // Two points at distance 1: one H0 bar born at 0, dead at the radius 1/2, at the midpoint of the edge.
    auto b = rips_persistence(make(10, 6, true, {{4.5, 3}, {5.5, 3}}), 1.0);
    auto f0 = finite(b, 0);
    check(f0.size() == 1 && n_ess(b, 0) == 1 && n_ess(b, 1) == 0 && finite(b, 1).empty(), "pair: one finite and one essential H0, no H1");
    check(f0.size() == 1 && std::abs(f0[0].birth) < 1e-15 && std::abs(f0[0].death - 0.5) < 1e-12 && std::abs(f0[0].xs - 5.0) < 1e-12 &&
              std::abs(f0[0].ys - 3.0) < 1e-12,
          "pair: bar [0, 1/2) with the death cell at the midpoint of the edge");

    // Minimum image in x: 0.2 and 9.8 in a box of width 10 are 0.4 apart; with walls they are 9.6 apart.
    b = rips_persistence(make(10, 6, true, {{0.2, 3}, {9.8, 3}}), 0.5);
    f0 = finite(b, 0);
    check(f0.size() == 1 && std::abs(f0[0].death - 0.2) < 1e-12 && std::abs(f0[0].xs) < 1e-12, "periodic: minimum image across x = 0 (midpoint at x = 0)");
    b = rips_persistence(make(10, 6, false, {{0.2, 3}, {9.8, 3}}), 0.5);
    check(finite(b, 0).empty() && n_ess(b, 0) == 2, "walls_x: no wrap, two components");

    // Unit square: four edges of length 1 close a loop at 1/2; the diagonals (sqrt 2) fill it at sqrt(2)/2.
    b = rips_persistence(make(20, 20, true, {{5, 5}, {6, 5}, {6, 6}, {5, 6}}), 1.0);
    auto f1 = finite(b, 1);
    check(finite(b, 0).size() == 3 && n_ess(b, 0) == 1, "square: three finite H0 bars and one essential");
    check(f1.size() == 1 && std::abs(f1[0].birth - 0.5) < 1e-12 && std::abs(f1[0].death - std::sqrt(0.5)) < 1e-12,
          "square: H1 bar [1/2, sqrt(2)/2)");
    check(f1.size() == 1 && std::abs(f1[0].x - 5.5) < 1e-12 + 0.5 && f1[0].xs > 5.0 && f1[0].xs < 6.0 && f1[0].ys > 5.0 && f1[0].ys < 6.0,
          "square: killing triangle inside the square");

    // Truncation: with max_radius below sqrt(2)/2 the loop is still alive and the bar is essential.
    b = rips_persistence(make(20, 20, true, {{5, 5}, {6, 5}, {6, 6}, {5, 6}}), 0.6);
    check(n_ess(b, 1) == 1 && finite(b, 1).empty(), "square, max_radius 0.6: the loop is essential");

    // Triangular lattice, a = 1, on a cylinder: every H0 bar dies at 1/2 and the only loop is the one around x.
    {
        std::vector<std::array<double, 2>> c;
        for (int j = 0; j < 6; ++j)
            for (int i = 0; i < 8; ++i) c.push_back({i + (j % 2 ? 0.5 : 0.0), 1.0 + j * std::sqrt(3.0) / 2});
        b = rips_persistence(make(8, 8, true, c), 1.0);
        bool all_half = true;
        for (const auto& x : finite(b, 0)) all_half = all_half && std::abs(x.death - 0.5) < 1e-12;
        check(finite(b, 0).size() == c.size() - 1 && all_half && n_ess(b, 0) == 1, "lattice: N - 1 H0 bars all dying at a/2");
        check(finite(b, 1).empty() && n_ess(b, 1) == 1, "lattice: no finite H1 (all triangles fill at a/2); one loop around x");
    }

    // Random points: the finite H0 deaths are half the edges of the minimum spanning tree (Prim, brute force).
    {
        std::mt19937 rng(11);
        std::uniform_real_distribution<double> ux(0, 12), uy(0, 8);
        std::vector<std::array<double, 2>> c;
        for (int i = 0; i < 120; ++i) c.push_back({ux(rng), uy(rng)});
        const ObsFile o = make(12, 8, true, c);
        auto d = [&](std::size_t i, std::size_t j) {
            double dx = c[j][0] - c[i][0];
            dx -= 12 * std::round(dx / 12);
            return std::hypot(dx, c[j][1] - c[i][1]);
        };
        std::vector<double> key(c.size(), 1e300), mst;
        std::vector<char> in(c.size(), 0);
        key[0] = 0;
        for (std::size_t it = 0; it < c.size(); ++it) {
            std::size_t u = c.size();
            for (std::size_t i = 0; i < c.size(); ++i) if (!in[i] && (u == c.size() || key[i] < key[u])) u = i;
            in[u] = 1;
            if (it) mst.push_back(key[u]);
            for (std::size_t i = 0; i < c.size(); ++i) if (!in[i]) key[i] = std::min(key[i], d(u, i));
        }
        std::sort(mst.begin(), mst.end());
        b = rips_persistence(o, 2.9);   // 4 * 2.9 < 12, and above every MST edge / 2 for 120 points in 96 units of area
        std::vector<double> got;
        for (const auto& x : finite(b, 0)) got.push_back(2 * x.death);
        std::sort(got.begin(), got.end());
        bool same = got.size() == mst.size();
        for (std::size_t i = 0; same && i < got.size(); ++i) same = std::abs(got[i] - mst[i]) < 1e-12;
        check(same && n_ess(b, 0) == 1, "random: finite H0 deaths = half the minimum-spanning-tree edges");
        bool ok = true;
        for (const auto& x : finite(b, 1)) ok = ok && x.persistence() > 0 && x.birth > 0;
        check(ok, "random: H1 bars have positive length");
    }

    // Errors and parameters.
    auto throws = [&](const ObsFile& o, double r) {
        try { rips_persistence(o, r); } catch (const std::invalid_argument&) { return true; }
        return false;
    };
    check(throws(make(10, 6, true, {{1, 1}, {2, 2}}), 2.5), "rips: max_radius too large for the periodic box");
    check(throws(make(10, 6, true, {}), 1.0), "rips: empty arrangement");
    check(std::abs(rips_default_radius(make(10, 6, true, {{1, 1}, {2, 2}, {3, 3}, {4, 4}, {5, 5}, {6, 5}})) - 1.5 * std::sqrt(10.0)) < 1e-12,
          "rips: default radius = 1.5 x mean spacing");
    const Params p = parse_params("filtration: rips\nmax_radius: 1.25\n");
    check(p.persistence.filtration == Filtration::rips && p.continuum.max_radius == 1.25 && describe(p) == "filtration=rips max_radius=1.25",
          "params: filtration rips, max_radius");
    bool bad = false;
    try { parse_params("max_radius: -1\n"); } catch (const std::runtime_error&) { bad = true; }
    check(bad, "params: negative max_radius is an error");
}

int main() {
    try {
        test_empty_torus();
        test_single_obstacle();
        test_blocks();
        test_connectivity();
        test_euler_random();
        test_signed();
        test_distances();
        test_readers();
        test_quantile_and_params();
        test_continuum();
        test_rips();
    } catch (const std::exception& e) {
        std::cout << "  FAIL exception: " << e.what() << '\n';
        ++failures;
    }
    std::cout << (failures ? std::format("{} check(s) failed\n", failures) : std::string("all checks passed\n"));
    return failures ? 1 : 0;
}
