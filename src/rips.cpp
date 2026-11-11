#include "rips.hpp"

#include <gudhi/Persistent_cohomology.h>
#include <gudhi/Rips_complex.h>
#include <gudhi/Simplex_tree.h>

#include <cmath>
#include <format>
#include <limits>
#include <numeric>
#include <stdexcept>

namespace p2d {

namespace {

using Complex = Gudhi::Simplex_tree<>;
using Field = Gudhi::persistent_cohomology::Field_Zp;
using Cohomology = Gudhi::persistent_cohomology::Persistent_cohomology<Complex, Field>;

constexpr double kNaN = std::numeric_limits<double>::quiet_NaN();

}  // namespace

double rips_default_radius(const ObsFile& obs) {
    if (obs.centres.empty()) throw std::invalid_argument("rips: no obstacles");
    return 1.5 * std::sqrt(obs.Lx * obs.Ly / static_cast<double>(obs.centres.size()));
}

std::vector<Bar> rips_persistence(const ObsFile& obs, double max_radius) {
    if (obs.centres.empty()) throw std::invalid_argument("rips: no obstacles");
    if (!(max_radius > 0)) throw std::invalid_argument("rips: max_radius must be positive");
    if (obs.periodic_x && !(4.0 * max_radius < obs.Lx))
        throw std::invalid_argument(std::format("rips: max_radius {} too large for a periodic box of width {} (need 4 max_radius < Lx)",
                                                max_radius, obs.Lx));
    const auto& c = obs.centres;
    const bool px = obs.periodic_x;
    const double Lx = obs.Lx;

    // dx reduced to the minimum image, so that |dx| <= Lx / 2.
    auto dx_of = [&](std::size_t i, std::size_t j) {
        double dx = c[j][0] - c[i][0];
        if (px) dx -= Lx * std::round(dx / Lx);
        return dx;
    };
    auto dist = [&](std::size_t i, std::size_t j) { return std::hypot(dx_of(i, j), c[j][1] - c[i][1]); };

    std::vector<std::size_t> ids(c.size());
    std::iota(ids.begin(), ids.end(), std::size_t{0});
    Gudhi::rips_complex::Rips_complex<double> rips(ids, 2.0 * max_radius, dist);
    Complex st;
    rips.create_complex(st, 2);   // triangles are enough for H0 and H1
    st.initialize_filtration();
    Cohomology pcoh(st, true);   // true: also the top dimension, so a complex with no triangles still reports its loops
    pcoh.init_coefficients(2);
    pcoh.compute_persistent_cohomology(0.0);

    // Centroid of a simplex, unwrapped around its first vertex and reduced into [0, Lx) if periodic.
    auto centroid = [&](Complex::Simplex_handle sh) {
        std::size_t first = 0, n = 0;
        double sx = 0.0, sy = 0.0;
        for (auto v : st.simplex_vertex_range(sh)) {
            const auto i = static_cast<std::size_t>(v);
            if (n == 0) first = i;
            sx += dx_of(first, i);
            sy += c[i][1];
            ++n;
        }
        double x = c[first][0] + sx / static_cast<double>(n);
        if (px) x -= Lx * std::floor(x / Lx);
        return std::array<double, 2>{x, sy / static_cast<double>(n)};
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
