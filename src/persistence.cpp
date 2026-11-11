#include "persistence.hpp"

#include <gudhi/Bitmap_cubical_complex.h>
#include <gudhi/Bitmap_cubical_complex_periodic_boundary_conditions_base.h>
#include <gudhi/Persistent_cohomology.h>

#include <stdexcept>

namespace p2d {

namespace {

using Base = Gudhi::cubical_complex::Bitmap_cubical_complex_periodic_boundary_conditions_base<double>;
using Complex = Gudhi::cubical_complex::Bitmap_cubical_complex<Base>;
using Field = Gudhi::persistent_cohomology::Field_Zp;
using Cohomology = Gudhi::persistent_cohomology::Persistent_cohomology<Complex, Field>;

// `free` mode: filtration value of the obstacles, above every free value (-w <= -1).
constexpr double kSentinel = 0.0;

}  // namespace

PersistenceResult compute_persistence(const Grid& g, const PersistenceOptions& opt) {
    if (opt.filtration == Filtration::height || opt.filtration == Filtration::centres || opt.filtration == Filtration::height_exact ||
        opt.filtration == Filtration::rips)
        throw std::invalid_argument("height, centres and rips are filtrations of .obs files");
    if (opt.connectivity != 4 && opt.connectivity != 8)
        throw std::invalid_argument("connectivity must be 4 or 8");
    if (g.nx < 2 || g.ny < 2) throw std::invalid_argument("the lattice must be at least 2 x 2");

    PersistenceResult res;
    res.signed_distance = opt.filtration == Filtration::signed_distance;
    res.w = distance_to_obstacles(g, opt.distance);
    if (res.signed_distance) {
        // Wall thickness: distance from each obstacle to the nearest free site (the boundary is not free space).
        Grid inverse = g;
        for (auto& c : inverse.obstacle) c = c ? 0 : 1;
        DistanceOptions inv = opt.distance;
        inv.boundary_obstacle = false;
        const auto d = distance_to_obstacles(inverse, inv);
        for (std::size_t i = 0; i < g.size(); ++i)
            if (g.obstacle[i]) res.w[i] = -d[i];
    }

    // GUDHI computes sublevel persistence: f = -w. In `free` mode the obstacles sit at the sentinel.
    std::vector<double> f(g.size());
    for (std::size_t i = 0; i < g.size(); ++i)
        f[i] = (g.obstacle[i] && !res.signed_distance) ? kSentinel : -res.w[i];

    const bool periodic = opt.distance.periodic;
    const std::vector<unsigned> dims{static_cast<unsigned>(g.nx), static_cast<unsigned>(g.ny)};  // x fastest
    const bool top_cells = opt.connectivity == 8;
    Complex cpx(dims, f, std::vector<bool>{periodic, periodic}, top_cells);

    // persistence_dim_max = true: without it GUDHI skips the top dimension (H2 of the whole torus).
    Cohomology pcoh(cpx, true);
    pcoh.init_coefficients(2);
    pcoh.compute_persistent_cohomology();

    const auto null = cpx.null_simplex();
    for (const auto& pair : pcoh.get_persistent_pairs()) {
        const auto birth_cell = std::get<0>(pair);
        const auto death_cell = std::get<1>(pair);
        const double fb = cpx.filtration(birth_cell);
        // `free` mode: a class born when the obstacles were added is not a class of the free space.
        if (!res.signed_distance && fb >= kSentinel) continue;
        Bar b;
        b.dim = cpx.dimension(birth_cell);
        b.birth = -fb;
        if (death_cell == null) {
            b.essential = true;
            b.wraps = b.dim == 1;  // survives into the whole domain: H1 of the torus (none if not periodic)
        } else {
            const double fd = cpx.filtration(death_cell);
            if (!res.signed_distance && fd >= kSentinel) b.essential = true;  // killed only by the obstacles
            else b.death = -fd;
        }
        res.bars.push_back(b);
    }
    return res;
}

int betti_at(const std::vector<Bar>& bars, int dim, double r) {
    int n = 0;
    for (const auto& b : bars)
        if (b.dim == dim && b.birth >= r && (b.essential || b.death < r)) ++n;
    return n;
}

int wraps_at(const std::vector<Bar>& bars, double r) {
    int n = 0;
    for (const auto& b : bars)
        if (b.wraps && b.birth >= r) ++n;
    return n;
}

}  // namespace p2d
