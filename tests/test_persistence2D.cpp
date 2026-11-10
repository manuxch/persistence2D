// Unit tests for persistence2D. No framework: each check prints its result and the exit code counts failures
// (run with ctest, or ./build/test_persistence2D).
#include "descriptors.hpp"
#include "distance.hpp"
#include "grid.hpp"
#include "output.hpp"
#include "params.hpp"
#include "persistence.hpp"

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
    } catch (const std::exception& e) {
        std::cout << "  FAIL exception: " << e.what() << '\n';
        ++failures;
    }
    std::cout << (failures ? std::format("{} check(s) failed\n", failures) : std::string("all checks passed\n"));
    return failures ? 1 : 0;
}
