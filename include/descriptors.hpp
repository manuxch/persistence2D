// Scalar summaries of a persistence result.
#pragma once

#include "grid.hpp"
#include "persistence.hpp"

#include <vector>

namespace p2d {

struct DimSummary {
    int n_finite = 0;            // finite bars
    int n_essential = 0;         // essential bars (never die in the free space)
    double per_free = 0.0;       // n_finite / free sites
    double pers_median = 0.0;    // persistence of the finite bars: median, 0.9 quantile, maximum
    double pers_q90 = 0.0;
    double pers_max = 0.0;
    double total = 0.0;          // total persistence: sum of finite persistences
    double entropy = 0.0;        // persistence entropy -sum p log p, p = persistence / total (finite bars)
};

struct Summary {
    std::size_t n_free = 0;
    double phi = 0.0;            // obstacle fraction
    double w_max = 0.0;          // widest free site
    double w_mean = 0.0;         // mean corridor width over the free sites
    double d_max = 0.0;          // thickest wall: largest distance from an obstacle to free space (`signed`)
    double r_min = 0.0;          // smallest w on a free site (1 unless there are no obstacles)
    int beta0 = 0;               // Betti numbers of the whole free space (r = r_min)
    int beta1 = 0;
    int wrap_dim = 0;            // independent wrapping H1 classes of the free space (0, 1 or 2)
    // Width at which V_r first wraps around the torus in one direction (r_wrap1) and in both (r_star). In `signed`
    // mode they are always defined on a torus, and negative when walls of thickness -r must be crossed; in `free`
    // mode they are NaN if the free space does not wrap.
    double r_wrap1 = 0.0;
    double r_star = 0.0;
    double branching = 0.0;      // beta1 / n_free
    DimSummary h0;
    DimSummary h1;
};

// Linear-interpolation quantile (the default of numpy.quantile); q in [0, 1]. NaN for an empty sample.
double quantile(std::vector<double> v, double q);

Summary summarize(const Grid& g, const PersistenceResult& res);

// Distinct filtration levels (values of w on the sites of the filtration), decreasing.
std::vector<double> levels(const Grid& g, const PersistenceResult& res);

}  // namespace p2d
