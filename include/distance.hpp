// Distance from each site to the nearest obstacle (the corridor width w of the free space).
#pragma once

#include "grid.hpp"

#include <vector>

namespace p2d {

enum class Metric { l1, l2 };

struct DistanceOptions {
    Metric metric = Metric::l1;
    bool periodic = true;            // minimum-image distances on the torus
    bool boundary_obstacle = false;  // open boundary only: the sites just outside the lattice count as obstacles
};

// w(x) = min over obstacles y of d(x, y): 0 on obstacles, >= 1 on free sites (L1 values are integers).
// Without any obstacle (and no boundary obstacle) every free site gets w = 1, so that the filtration is a
// single level and the result is the topology of the whole domain.
std::vector<double> distance_to_obstacles(const Grid& g, const DistanceOptions& opt);

}  // namespace p2d
