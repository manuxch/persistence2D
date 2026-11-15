// Vietoris-Rips persistence of a point pattern: the obstacle centres of a .obs arrangement, or the obstacle (or
// free) sites of a lattice map.
#pragma once

#include "continuum.hpp"
#include "grid.hpp"

#include <array>
#include <vector>

namespace p2d {

// Points in an Lx x Ly box, periodic (minimum image) in x and/or y.
struct PointCloud {
    double Lx = 0.0, Ly = 0.0;
    bool periodic_x = false, periodic_y = false;
    std::vector<std::array<double, 2>> points;
};

PointCloud point_cloud(const ObsFile& obs);  // periodic in x if periodic_x, never in y
// Sites of a lattice map as points at their integer coordinates (x, y), in an nx x ny box periodic in both
// directions or in none. `obstacles` = true takes the obstacle sites, false the free sites.
PointCloud point_cloud(const Grid& g, bool periodic, bool obstacles);

// Rips complex of the centres, filtered by radius r = (edge length) / 2, so that its scale is that of the union of
// discs of radius r around the centres (and of the `centres` and alpha filtrations). A simplex enters at half the
// longest of its edges. Distances are Euclidean, with the minimum image in x if the arrangement is periodic_x
// (the metric of a cylinder; y is open). Only edges up to 2 * max_radius are kept, so the complex is the 2-skeleton
// clique complex of that graph; H0 and H1 are exact up to max_radius.
//
// Bars: H0 are born at 0 and die when two components join, at half the edge that joins them (x, y: the vertex that
// is born; xs, ys: midpoint of that edge). H1 are born at half the edge that closes a loop and die at half the
// longest edge of the triangle that fills it (x, y and xs, ys: centroids of the birth edge and of the killing
// triangle, reduced into the box in the periodic directions). `essential` bars are still alive at max_radius: one
// H0 per component of the graph, and H1 classes that the truncated complex has not filled (including the loops
// that wind around a periodic direction). Needs 4 max_radius < L in every periodic direction, so that the minimum
// image is unique.
std::vector<Bar> rips_persistence(const PointCloud& pc, double max_radius);
std::vector<Bar> rips_persistence(const ObsFile& obs, double max_radius);

// Default max_radius: 1.5 times the mean spacing sqrt(Lx Ly / N).
double rips_default_radius(const PointCloud& pc);
double rips_default_radius(const ObsFile& obs);

}  // namespace p2d
