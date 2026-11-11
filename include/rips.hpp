// Vietoris-Rips persistence of the obstacle centres of a .obs arrangement (point-pattern topology).
#pragma once

#include "continuum.hpp"

#include <vector>

namespace p2d {

// Rips complex of the centres, filtered by radius r = (edge length) / 2, so that its scale is that of the union of
// discs of radius r around the centres (and of the `centres` and alpha filtrations). A simplex enters at half the
// longest of its edges. Distances are Euclidean, with the minimum image in x if the arrangement is periodic_x
// (the metric of a cylinder; y is open). Only edges up to 2 * max_radius are kept, so the complex is the 2-skeleton
// clique complex of that graph; H0 and H1 are exact up to max_radius.
//
// Bars: H0 are born at 0 and die when two components join, at half the edge that joins them (x, y: the vertex that
// is born; xs, ys: midpoint of that edge). H1 are born at half the edge that closes a loop and die at half the
// longest edge of the triangle that fills it (x, y and xs, ys: centroids of the birth edge and of the killing
// triangle, reduced into [0, Lx) if periodic). `essential` bars are still alive at max_radius: one H0 per
// component of the graph, and H1 classes that the truncated complex has not filled (including, on a cylinder, the
// loop around x). Needs 2 * max_radius < Lx / 2 if periodic, so that the minimum image is unique.
std::vector<Bar> rips_persistence(const ObsFile& obs, double max_radius);

// Default max_radius: 1.5 times the mean spacing sqrt(Lx Ly / N).
double rips_default_radius(const ObsFile& obs);

}  // namespace p2d
