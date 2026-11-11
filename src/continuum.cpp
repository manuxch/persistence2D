#include "continuum.hpp"

#include <gudhi/Bitmap_cubical_complex.h>
#include <gudhi/Bitmap_cubical_complex_periodic_boundary_conditions_base.h>
#include <gudhi/Persistent_cohomology.h>

#include <algorithm>
#include <cmath>
#include <format>
#include <fstream>
#include <limits>
#include <sstream>
#include <stdexcept>

namespace p2d {

namespace {

using Base = Gudhi::cubical_complex::Bitmap_cubical_complex_periodic_boundary_conditions_base<double>;
using Complex = Gudhi::cubical_complex::Bitmap_cubical_complex<Base>;
using Field = Gudhi::persistent_cohomology::Field_Zp;
using Cohomology = Gudhi::persistent_cohomology::Persistent_cohomology<Complex, Field>;

constexpr double kInf = std::numeric_limits<double>::infinity();
constexpr double kNaN = std::numeric_limits<double>::quiet_NaN();

std::string trim(const std::string& s) {
    const auto b = s.find_first_not_of(" \t\r");
    if (b == std::string::npos) return {};
    return s.substr(b, s.find_last_not_of(" \t\r") - b + 1);
}

// Cell bucket grid over the box, to find the nearest centre by expanding rings.
struct Buckets {
    double Lx, Ly, cell;
    int ncx, ncy;
    bool px;
    std::vector<std::vector<std::array<double, 2>>> items;

    Buckets(const ObsFile& o, double c) : Lx(o.Lx), Ly(o.Ly), px(o.periodic_x) {
        ncx = std::max(1, static_cast<int>(Lx / c));
        ncy = std::max(1, static_cast<int>(Ly / c));
        cell = c;
        items.resize(static_cast<std::size_t>(ncx) * static_cast<std::size_t>(ncy));
        for (const auto& p : o.centres) {
            double x = p[0];
            if (px) x -= Lx * std::floor(x / Lx);
            items[index(cx(x), cy(p[1]))].push_back({x, p[1]});
        }
    }
    int cx(double x) const { return std::clamp(static_cast<int>(x / Lx * ncx), 0, ncx - 1); }
    int cy(double y) const { return std::clamp(static_cast<int>(y / Ly * ncy), 0, ncy - 1); }
    std::size_t index(int i, int j) const { return static_cast<std::size_t>(i) + static_cast<std::size_t>(ncx) * static_cast<std::size_t>(j); }

    double nearest(double x, double y) const {
        const int i0 = cx(x), j0 = cy(y);
        const double wx = Lx / ncx, wy = Ly / ncy, wmin = std::min(wx, wy);
        double best2 = kInf;
        const int kmax = std::max(ncx, ncy);
        for (int k = 0; k <= kmax; ++k) {
            // Any point in ring k is at least (k - 1) cells away; stop once that exceeds the best distance.
            if (k >= 1 && std::sqrt(best2) <= (k - 1) * wmin) break;
            for (int dj = -k; dj <= k; ++dj) {
                const int j = j0 + dj;
                if (j < 0 || j >= ncy) continue;
                const bool edge_row = std::abs(dj) == k;
                for (int di = -k; di <= k; di += (edge_row ? 1 : std::max(1, 2 * k))) {
                    int i = i0 + di;
                    if (px) i = ((i % ncx) + ncx) % ncx;
                    else if (i < 0 || i >= ncx) continue;
                    // With periodic wrapping a ring wider than the box visits cells twice: harmless.
                    for (const auto& p : items[index(i, j)]) {
                        double dx = std::abs(p[0] - x);
                        if (px) dx = std::min(dx, Lx - dx);
                        const double dy = p[1] - y;
                        best2 = std::min(best2, dx * dx + dy * dy);
                    }
                }
            }
        }
        return std::sqrt(best2);
    }
};

struct RawPair {
    int dim;
    double fb, fd;
    bool essential;
    double xb, yb, xs, ys;   // cell centres in pixel units (vertex (i, j) is at (i, j))
};

// Sublevel persistence of f on the nx x ny vertices (vertex construction, 4-neighbourhood), periodic or not in
// each direction. f = +inf marks a vertex outside the filtration. Bars born at the sentinel are dropped and
// bars that die there are essential.
std::vector<RawPair> sublevel(std::size_t nx, std::size_t ny, std::vector<double> f, bool px, double min_pers) {
    double top = -kInf;
    for (double v : f)
        if (v != kInf) top = std::max(top, v);
    if (top == -kInf) return {};
    const double sentinel = top + 1.0;
    for (auto& v : f)
        if (v == kInf) v = sentinel;

    const std::vector<unsigned> dims{static_cast<unsigned>(nx), static_cast<unsigned>(ny)};
    Complex cpx(dims, f, std::vector<bool>{px, false}, false);
    Cohomology pcoh(cpx, true);
    pcoh.init_coefficients(2);
    pcoh.compute_persistent_cohomology(min_pers);

    // Cell layout of GUDHI: per axis, 2n cells if periodic, 2n - 1 for a vertex construction otherwise.
    const std::size_t ex = px ? 2 * nx : 2 * nx - 1;
    auto decode = [&](std::size_t cell, double& x, double& y) {
        x = 0.5 * static_cast<double>(cell % ex);
        y = 0.5 * static_cast<double>(cell / ex);
    };

    std::vector<RawPair> out;
    const auto null = cpx.null_simplex();
    for (const auto& pair : pcoh.get_persistent_pairs()) {
        const auto bc = std::get<0>(pair), dc = std::get<1>(pair);
        const double fb = cpx.filtration(bc);
        if (fb >= sentinel) continue;
        RawPair r{};
        r.dim = static_cast<int>(cpx.dimension(bc));
        r.fb = fb;
        decode(static_cast<std::size_t>(bc), r.xb, r.yb);
        if (dc == null) {
            r.essential = true;
            r.xs = r.ys = kNaN;
        } else {
            r.fd = cpx.filtration(dc);
            if (r.fd >= sentinel) {
                r.essential = true;
                r.xs = r.ys = kNaN;
            } else {
                decode(static_cast<std::size_t>(dc), r.xs, r.ys);
            }
        }
        out.push_back(r);
    }
    return out;
}

}  // namespace

bool is_obs_file(const std::string& path) {
    std::ifstream in(path);
    std::string line;
    return in && std::getline(in, line) && trim(line).rfind("# format: 1", 0) == 0;
}

ObsFile read_obs(const std::string& path) {
    std::ifstream in(path);
    if (!in) throw std::runtime_error(std::format("{}: cannot open", path));
    ObsFile o;
    bool have_lx = false, have_ly = false, have_r = false, have_b = false;
    std::string line;
    long declared = -1;
    for (int lineno = 1; std::getline(in, line); ++lineno) {
        const std::string t = trim(line);
        if (t.empty()) continue;
        if (t[0] == '#') {
            const auto colon = t.find(':');
            if (colon == std::string::npos) continue;
            const std::string key = trim(t.substr(1, colon - 1)), val = trim(t.substr(colon + 1));
            try {
                if (key == "Lx") { o.Lx = std::stod(val); have_lx = true; }
                else if (key == "Ly") { o.Ly = std::stod(val); have_ly = true; }
                else if (key == "R") { o.R = std::stod(val); have_r = true; }
                else if (key == "n_obstacles") declared = std::stol(val);
                else if (key == "sha256_body") o.sha256_body = val;
                else if (key == "boundary") {
                    if (val == "periodic_x") o.periodic_x = true;
                    else if (val == "walls_x") o.periodic_x = false;
                    else throw std::runtime_error(std::format("{}: unsupported boundary '{}'", path, val));
                    have_b = true;
                }
            } catch (const std::invalid_argument&) {
                throw std::runtime_error(std::format("{}:{}: bad value for '{}'", path, lineno, key));
            }
            continue;
        }
        std::istringstream ls(t);
        double x, y;
        if (!(ls >> x >> y)) throw std::runtime_error(std::format("{}:{}: expected 'x y'", path, lineno));
        o.centres.push_back({x, y});
    }
    if (!(have_lx && have_ly && have_r && have_b))
        throw std::runtime_error(std::format("{}: header lacks Lx, Ly, R or boundary", path));
    if (!(o.Lx > 0 && o.Ly > 0 && o.R > 0)) throw std::runtime_error(std::format("{}: Lx, Ly and R must be positive", path));
    if (declared >= 0 && static_cast<std::size_t>(declared) != o.centres.size())
        throw std::runtime_error(std::format("{}: header declares {} obstacles, file has {}", path, declared, o.centres.size()));
    return o;
}

std::vector<double> centre_distance(const ObsFile& obs, std::size_t nx, std::size_t ny) {
    if (obs.centres.empty()) throw std::invalid_argument("centre_distance: no obstacles");
    const double area = obs.Lx * obs.Ly;
    const double spacing = std::sqrt(area / static_cast<double>(obs.centres.size()));
    const Buckets b(obs, spacing);
    std::vector<double> d(nx * ny);
    for (std::size_t j = 0; j < ny; ++j) {
        const double y = (static_cast<double>(j) + 0.5) * obs.Ly / static_cast<double>(ny);
        for (std::size_t i = 0; i < nx; ++i) {
            const double x = (static_cast<double>(i) + 0.5) * obs.Lx / static_cast<double>(nx);
            d[i + nx * j] = b.nearest(x, y);
        }
    }
    return d;
}

ContinuumResult compute_continuum(const ObsFile& obs, const ContinuumOptions& opt) {
    if (!(opt.pixel > 0) || !(opt.grain > 0)) throw std::invalid_argument("pixel and grain must be positive");
    ContinuumResult res;
    if (opt.mode == ContinuumMode::height_exact) {
        res.clearance = opt.clearance >= 0 ? opt.clearance : obs.R + 0.5 * opt.grain;
        res.bars = height_exact(obs, res.clearance);
        return res;
    }
    res.nx = std::max<std::size_t>(2, static_cast<std::size_t>(std::llround(obs.Lx / opt.pixel)));
    res.ny = std::max<std::size_t>(2, static_cast<std::size_t>(std::llround(obs.Ly / opt.pixel)));
    const double hx = obs.Lx / static_cast<double>(res.nx), hy = obs.Ly / static_cast<double>(res.ny);
    const bool height = opt.mode == ContinuumMode::height;
    res.min_persistence = opt.min_persistence;
    res.clearance = opt.clearance >= 0 ? opt.clearance : (height ? obs.R + 0.5 * opt.grain : 0.0);

    const auto dist = centre_distance(obs, res.nx, res.ny);
    std::vector<double> f(res.nx * res.ny, kInf);
    for (std::size_t j = 0; j < res.ny; ++j)
        for (std::size_t i = 0; i < res.nx; ++i) {
            const std::size_t k = i + res.nx * j;
            const double x = (static_cast<double>(i) + 0.5) * hx, y = (static_cast<double>(j) + 0.5) * hy;
            if (dist[k] < res.clearance * (1.0 - 1e-12)) continue;
            if (height) {
                if (!obs.periodic_x && (x < 0.5 * opt.grain || obs.Lx - x < 0.5 * opt.grain)) continue;
                f[k] = y;
            } else {
                f[k] = -dist[k];
            }
            ++res.n_access;
        }

    const bool px = obs.periodic_x;
    for (const auto& r : sublevel(res.nx, res.ny, std::move(f), px, res.min_persistence)) {
        Bar b;
        b.dim = r.dim;
        b.essential = r.essential;
        b.birth = height ? r.fb : -r.fb;
        b.death = r.essential ? 0.0 : (height ? r.fd : -r.fd);
        b.wraps = r.essential && r.dim == 1 && px;
        auto to_x = [&](double c) {   // vertex coordinate -> length, reduced into [0, Lx) if periodic
            double x = (c + 0.5) * hx;
            return px ? x - obs.Lx * std::floor(x / obs.Lx) : x;
        };
        b.x = to_x(r.xb);
        b.y = (r.yb + 0.5) * hy;
        b.xs = std::isnan(r.xs) ? kNaN : to_x(r.xs);
        b.ys = std::isnan(r.ys) ? kNaN : (r.ys + 0.5) * hy;
        res.bars.push_back(b);
    }
    return res;
}

}  // namespace p2d
