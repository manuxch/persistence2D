// Persistent homology of the free space of a binary map, filtered by corridor width.
#pragma once

#include "distance.hpp"
#include "grid.hpp"

#include <cmath>
#include <limits>
#include <vector>

namespace p2d {

// free:   superlevel filtration of the corridor width w on the free sites only, V_r = {free x : w(x) >= r},
//         r = w_max ... 1. It describes the free space alone.
// signed: the same filtration continued into the obstacles with the signed width s, equal to w on free sites
//         and to -d on obstacles, d = distance to the nearest free site (wall thickness). V_r = {x : s(x) >= r}
//         runs from r = w_max down to -d_max (the whole domain). For r >= 1 it coincides with `free`; below,
//         obstacles are absorbed from their surface inwards, so free clusters merge and loops close across
//         walls of thickness |r| (Robins et al., Water Resour. Res. 52, 315 (2016), signed distance transform).
enum class Filtration { free, signed_distance, height, centres };
// `height` and `centres` are the continuum filtrations of continuum.hpp (input: .obs files).

struct PersistenceOptions {
    DistanceOptions distance;
    Filtration filtration = Filtration::free;
    // 4: free sites are the vertices of the cubical complex, so two free sites are joined only through a shared
    //    edge (the 4-neighbourhood of lattice dynamics). 8: free sites are pixels (top cells), so diagonal
    //    neighbours touch at a corner and are connected.
    int connectivity = 4;
};

// One persistence interval of the superlevel filtration V_r. The class appears at width `birth` and
// disappears at the smaller width `death`. `essential` classes never die: in `free` mode, inside the free space;
// in `signed` mode, in the whole domain. On a periodic lattice, `wraps` marks the H1 classes that are
// non-trivial on the torus (they wind around it); there are at most two independent ones.
struct Bar {
    int dim = 0;
    double birth = 0.0;
    double death = 0.0;  // meaningless when essential
    bool essential = false;
    bool wraps = false;
    // Continuum filtrations only: position of the birth cell (x, y) and of the death cell (xs, ys; NaN if essential).
    double x = std::numeric_limits<double>::quiet_NaN(), y = x, xs = x, ys = x;
    double persistence() const { return std::abs(birth - death); }
};

struct PersistenceResult {
    // Filtration value of every site: w >= 1 on free sites; on obstacles 0 (`free`, not in the filtration) or
    // -d <= -1 (`signed`).
    std::vector<double> w;
    bool signed_distance = false;
    std::vector<Bar> bars;   // all intervals, H0 and H1 (and H2 if the whole torus is reached)

    bool in_filtration(const Grid& g, std::size_t i) const { return signed_distance || !g.obstacle[i]; }
};

// Persistence of V_r. In `free` mode the obstacles are appended at a sentinel level after the lowest free
// level, so that the last complex is the whole domain (the torus if periodic); its persistence truncated at
// the sentinel is that of the free space: bars born at the sentinel are dropped, bars alive there are
// essential. In both modes, the bars that never die in the whole torus are the classes that wind around it,
// and their birth is the widest corridor (or, if negative, the thinnest wall to cross) at which the
// domain wraps in that direction.
PersistenceResult compute_persistence(const Grid& g, const PersistenceOptions& opt);

// Betti numbers of V_r: number of bars alive at width r (birth >= r > death, or essential with birth >= r).
int betti_at(const std::vector<Bar>& bars, int dim, double r);
// Number of independent wrapping H1 classes in V_r (the wrap dimension of V_r).
int wraps_at(const std::vector<Bar>& bars, double r);

}  // namespace p2d
