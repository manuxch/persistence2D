#include "descriptors.hpp"

#include <algorithm>
#include <cmath>
#include <functional>
#include <limits>

namespace p2d {

namespace {

constexpr double kNaN = std::numeric_limits<double>::quiet_NaN();

DimSummary summarize_dim(const std::vector<Bar>& bars, int dim, std::size_t n_free) {
    DimSummary s;
    std::vector<double> p;
    for (const auto& b : bars) {
        if (b.dim != dim) continue;
        if (b.essential) { ++s.n_essential; continue; }
        p.push_back(b.persistence());
    }
    s.n_finite = static_cast<int>(p.size());
    s.per_free = n_free ? static_cast<double>(p.size()) / static_cast<double>(n_free) : kNaN;
    if (p.empty()) {
        s.pers_median = s.pers_q90 = s.pers_max = kNaN;
        return s;
    }
    s.pers_median = quantile(p, 0.5);
    s.pers_q90 = quantile(p, 0.9);
    s.pers_max = *std::max_element(p.begin(), p.end());
    for (double x : p) s.total += x;
    for (double x : p) {
        const double q = x / s.total;
        if (q > 0) s.entropy -= q * std::log(q);
    }
    return s;
}

}  // namespace

double quantile(std::vector<double> v, double q) {
    if (v.empty()) return kNaN;
    std::sort(v.begin(), v.end());
    const double pos = q * static_cast<double>(v.size() - 1);
    const auto lo = static_cast<std::size_t>(std::floor(pos));
    const std::size_t hi = std::min(lo + 1, v.size() - 1);
    return v[lo] + (pos - static_cast<double>(lo)) * (v[hi] - v[lo]);
}

std::vector<double> levels(const Grid& g, const PersistenceResult& res) {
    std::vector<double> lv;
    for (std::size_t i = 0; i < g.size(); ++i)
        if (res.in_filtration(g, i)) lv.push_back(res.w[i]);
    std::sort(lv.begin(), lv.end(), std::greater<>());
    lv.erase(std::unique(lv.begin(), lv.end()), lv.end());
    return lv;
}

Summary summarize(const Grid& g, const PersistenceResult& res) {
    Summary s;
    s.n_free = g.n_free();
    s.phi = static_cast<double>(g.n_obstacles()) / static_cast<double>(g.size());
    s.w_max = s.r_min = kNaN;
    s.w_mean = kNaN;
    if (s.n_free > 0) {
        double sum = 0.0, mx = 0.0, mn = std::numeric_limits<double>::infinity();
        for (std::size_t i = 0; i < g.size(); ++i) {
            if (g.obstacle[i]) continue;
            sum += res.w[i];
            mx = std::max(mx, res.w[i]);
            mn = std::min(mn, res.w[i]);
        }
        s.w_max = mx;
        s.w_mean = sum / static_cast<double>(s.n_free);
        s.r_min = mn;
        s.beta0 = betti_at(res.bars, 0, mn);
        s.beta1 = betti_at(res.bars, 1, mn);
        s.wrap_dim = wraps_at(res.bars, mn);
    }
    std::vector<double> wrap_births;
    for (const auto& b : res.bars)
        if (b.wraps) wrap_births.push_back(b.birth);
    std::sort(wrap_births.begin(), wrap_births.end(), std::greater<>());
    s.r_wrap1 = wrap_births.size() >= 1 ? wrap_births[0] : kNaN;
    s.r_star = wrap_births.size() >= 2 ? wrap_births[1] : kNaN;
    s.d_max = 0.0;
    if (res.signed_distance)
        for (std::size_t i = 0; i < g.size(); ++i)
            if (g.obstacle[i]) s.d_max = std::max(s.d_max, -res.w[i]);
    s.branching = s.n_free ? static_cast<double>(s.beta1) / static_cast<double>(s.n_free) : kNaN;
    s.h0 = summarize_dim(res.bars, 0, s.n_free);
    s.h1 = summarize_dim(res.bars, 1, s.n_free);
    return s;
}

}  // namespace p2d
