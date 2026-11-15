#include "rips.hpp"

#include <gudhi/Persistent_cohomology.h>
#include <gudhi/Rips_complex.h>
#include <gudhi/Simplex_tree.h>

#include <cmath>
#include <format>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <tuple>

namespace p2d {

namespace {

using Complex = Gudhi::Simplex_tree<>;
using Field = Gudhi::persistent_cohomology::Field_Zp;
using Cohomology = Gudhi::persistent_cohomology::Persistent_cohomology<Complex, Field>;

constexpr double kNaN = std::numeric_limits<double>::quiet_NaN();

}  // namespace

PointCloud point_cloud(const ObsFile& obs) {
    PointCloud pc;
    pc.Lx = obs.Lx;
    pc.Ly = obs.Ly;
    pc.periodic_x = obs.periodic_x;
    pc.periodic_y = false;
    pc.points = obs.centres;
    return pc;
}

PointCloud point_cloud(const Grid& g, bool periodic, bool obstacles) {
    PointCloud pc;
    pc.Lx = static_cast<double>(g.nx);
    pc.Ly = static_cast<double>(g.ny);
    pc.periodic_x = pc.periodic_y = periodic;
    for (std::size_t y = 0; y < g.ny; ++y)
        for (std::size_t x = 0; x < g.nx; ++x)
            if ((g.obstacle[g.index(x, y)] != 0) == obstacles)
                pc.points.push_back({static_cast<double>(x), static_cast<double>(y)});
    return pc;
}

double rips_default_radius(const PointCloud& pc) {
    if (pc.points.empty()) throw std::invalid_argument("rips: no points");
    return 1.5 * std::sqrt(pc.Lx * pc.Ly / static_cast<double>(pc.points.size()));
}

double rips_default_radius(const ObsFile& obs) {
    if (obs.centres.empty()) throw std::invalid_argument("rips: no obstacles");
    return rips_default_radius(point_cloud(obs));
}

std::vector<Bar> rips_persistence(const ObsFile& obs, double max_radius) {
    if (obs.centres.empty()) throw std::invalid_argument("rips: no obstacles");
    return rips_persistence(point_cloud(obs), max_radius);
}

std::vector<Bar> rips_persistence(const PointCloud& pc, double max_radius) {
    if (pc.points.empty()) throw std::invalid_argument("rips: no points");
    if (!(max_radius > 0)) throw std::invalid_argument("rips: max_radius must be positive");
    for (const auto& [periodic, L, name] : {std::tuple{pc.periodic_x, pc.Lx, "width"}, std::tuple{pc.periodic_y, pc.Ly, "height"}})
        // 6 max_radius < L: with edges up to 2 max_radius >= L / 3, three edges can go round the periodic direction
        // and fill its loops (and 4 max_radius < L keeps the minimum image unique).
        if (periodic && !(6.0 * max_radius < L))
            throw std::invalid_argument(std::format("rips: max_radius {} too large for a periodic box of {} {} (need 6 max_radius < L)",
                                                    max_radius, name, L));
    const auto& c = pc.points;
    const bool px = pc.periodic_x, py = pc.periodic_y;
    const double Lx = pc.Lx, Ly = pc.Ly;

    // Coordinate differences reduced to the minimum image, so that |d| <= L / 2 in the periodic directions.
    auto dx_of = [&](std::size_t i, std::size_t j) {
        double dx = c[j][0] - c[i][0];
        if (px) dx -= Lx * std::round(dx / Lx);
        return dx;
    };
    auto dy_of = [&](std::size_t i, std::size_t j) {
        double dy = c[j][1] - c[i][1];
        if (py) dy -= Ly * std::round(dy / Ly);
        return dy;
    };
    auto dist = [&](std::size_t i, std::size_t j) { return std::hypot(dx_of(i, j), dy_of(i, j)); };

    std::vector<std::size_t> ids(c.size());
    std::iota(ids.begin(), ids.end(), std::size_t{0});
    Gudhi::rips_complex::Rips_complex<double> rips(ids, 2.0 * max_radius, dist);
    Complex st;
    rips.create_complex(st, 2);   // triangles are enough for H0 and H1
    st.initialize_filtration();
    Cohomology pcoh(st, true);   // true: also the top dimension, so a complex with no triangles still reports its loops
    pcoh.init_coefficients(2);
    pcoh.compute_persistent_cohomology(0.0);

    // Centroid of a simplex, unwrapped around its first vertex and reduced into the box in the periodic directions.
    auto centroid = [&](Complex::Simplex_handle sh) {
        std::size_t first = 0, n = 0;
        double sx = 0.0, sy = 0.0;
        for (auto v : st.simplex_vertex_range(sh)) {
            const auto i = static_cast<std::size_t>(v);
            if (n == 0) first = i;
            sx += dx_of(first, i);
            sy += dy_of(first, i);
            ++n;
        }
        double x = c[first][0] + sx / static_cast<double>(n);
        double y = c[first][1] + sy / static_cast<double>(n);
        if (px) x -= Lx * std::floor(x / Lx);
        if (py) y -= Ly * std::floor(y / Ly);
        return std::array<double, 2>{x, y};
    };

    std::vector<Bar> bars;
    for (const auto& pair : pcoh.get_persistent_pairs()) {
        const auto birth = std::get<0>(pair), death = std::get<1>(pair);
        Bar b;
        b.dim = st.dimension(birth);
        if (b.dim > 1) continue;
        b.birth = 0.5 * st.filtration(birth);
        b.essential = death == st.null_simplex();
        // A pair of equal length (up to rounding of the distances) is not a feature.
        if (!b.essential && st.filtration(death) - st.filtration(birth) <= 1e-12 * (1.0 + st.filtration(birth))) continue;
        const auto p = centroid(birth);
        b.x = p[0];
        b.y = p[1];
        if (b.essential) {
            b.death = 0.0;
        } else {
            b.death = 0.5 * st.filtration(death);
            const auto q = centroid(death);
            b.xs = q[0];
            b.ys = q[1];
        }
        bars.push_back(b);
    }
    return bars;
}

}  // namespace p2d
