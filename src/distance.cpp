#include "distance.hpp"

#include <cmath>
#include <cstdint>
#include <limits>
#include <stdexcept>

namespace p2d {

namespace {

constexpr double kInf = std::numeric_limits<double>::infinity();

// Pads an open-boundary grid with one ring of obstacles (boundary_obstacle), so that both metrics can treat
// the boundary as ordinary obstacles.
Grid padded(const Grid& g) {
    Grid p;
    p.nx = g.nx + 2;
    p.ny = g.ny + 2;
    p.obstacle.assign(p.size(), 1);
    for (std::size_t y = 0; y < g.ny; ++y)
        for (std::size_t x = 0; x < g.nx; ++x) p.obstacle[p.index(x + 1, y + 1)] = g.obstacle[g.index(x, y)];
    return p;
}

std::vector<double> unpad(const Grid& g, const std::vector<double>& wp) {
    std::vector<double> w(g.size());
    const std::size_t px = g.nx + 2;
    for (std::size_t y = 0; y < g.ny; ++y)
        for (std::size_t x = 0; x < g.nx; ++x) w[g.index(x, y)] = wp[(x + 1) + px * (y + 1)];
    return w;
}

// Multi-source breadth-first search from the obstacles with the 4-neighbourhood: the graph distance of the
// lattice, which is the L1 distance (minimum image on the torus).
std::vector<double> l1_bfs(const Grid& g, bool periodic) {
    const std::size_t n = g.size();
    constexpr std::uint32_t unset = std::numeric_limits<std::uint32_t>::max();
    std::vector<std::uint32_t> d(n, unset);
    std::vector<std::size_t> queue;
    queue.reserve(n);
    for (std::size_t i = 0; i < n; ++i)
        if (g.obstacle[i]) { d[i] = 0; queue.push_back(i); }
    for (std::size_t head = 0; head < queue.size(); ++head) {
        const std::size_t i = queue[head];
        const std::size_t x = i % g.nx, y = i / g.nx;
        auto visit = [&](std::size_t j) {
            if (d[j] == unset) { d[j] = d[i] + 1; queue.push_back(j); }
        };
        if (x + 1 < g.nx) visit(i + 1); else if (periodic) visit(i + 1 - g.nx);
        if (x > 0) visit(i - 1); else if (periodic) visit(i + g.nx - 1);
        if (y + 1 < g.ny) visit(i + g.nx); else if (periodic) visit(x);
        if (y > 0) visit(i - g.nx); else if (periodic) visit(x + g.nx * (g.ny - 1));
    }
    std::vector<double> w(n);
    for (std::size_t i = 0; i < n; ++i) w[i] = d[i] == unset ? kInf : static_cast<double>(d[i]);
    return w;
}

// 1D squared Euclidean distance transform (Felzenszwalb & Huttenlocher, Theory Comput. 8, 415 (2012)):
// out[q] = min_p ((q - p)^2 + f[p]), with f = +inf where there is no source.
void dt1d(const std::vector<double>& f, std::vector<double>& out, std::vector<std::size_t>& v,
          std::vector<double>& z) {
    const std::size_t n = f.size();
    out.assign(n, kInf);
    v.resize(n);
    z.resize(n + 1);
    std::size_t k = 0;
    std::size_t first = 0;
    while (first < n && f[first] == kInf) ++first;
    if (first == n) return;
    v[0] = first;
    z[0] = -kInf;
    z[1] = kInf;
    for (std::size_t q = first + 1; q < n; ++q) {
        if (f[q] == kInf) continue;
        const double qd = static_cast<double>(q);
        double s;
        while (true) {
            const double vk = static_cast<double>(v[k]);
            s = ((f[q] + qd * qd) - (f[v[k]] + vk * vk)) / (2.0 * qd - 2.0 * vk);
            if (s <= z[k]) --k;  // z[0] = -inf stops the loop at k = 0
            else break;
        }
        ++k;
        v[k] = q;
        z[k] = s;
        z[k + 1] = kInf;
    }
    k = 0;
    for (std::size_t q = 0; q < n; ++q) {
        const double qd = static_cast<double>(q);
        while (z[k + 1] < qd) ++k;
        const double dv = qd - static_cast<double>(v[k]);
        out[q] = dv * dv + f[v[k]];
    }
}

// Squared Euclidean distance transform along one line of length n (stride between elements), periodic or
// open. The periodic case replicates the line three times, which is exact because the nearest image of any
// source is within n/2 < n of the central copy.
void edt_line(std::vector<double>& data, std::size_t start, std::size_t stride, std::size_t n, bool periodic,
              std::vector<double>& f, std::vector<double>& out, std::vector<std::size_t>& v, std::vector<double>& z) {
    const std::size_t m = periodic ? 3 * n : n;
    f.resize(m);
    for (std::size_t i = 0; i < m; ++i) f[i] = data[start + stride * (i % n)];
    dt1d(f, out, v, z);
    const std::size_t offset = periodic ? n : 0;
    for (std::size_t i = 0; i < n; ++i) data[start + stride * i] = out[offset + i];
}

std::vector<double> l2_edt(const Grid& g, bool periodic) {
    std::vector<double> d2(g.size());
    for (std::size_t i = 0; i < g.size(); ++i) d2[i] = g.obstacle[i] ? 0.0 : kInf;
    std::vector<double> f, out, z;
    std::vector<std::size_t> v;
    for (std::size_t x = 0; x < g.nx; ++x) edt_line(d2, x, g.nx, g.ny, periodic, f, out, v, z);       // columns
    for (std::size_t y = 0; y < g.ny; ++y) edt_line(d2, g.nx * y, 1, g.nx, periodic, f, out, v, z);   // rows
    for (auto& x : d2) x = std::sqrt(x);
    return d2;
}

}  // namespace

std::vector<double> distance_to_obstacles(const Grid& g, const DistanceOptions& opt) {
    if (g.size() == 0) throw std::invalid_argument("distance_to_obstacles: empty grid");
    const bool pad = !opt.periodic && opt.boundary_obstacle;
    const Grid padded_grid = pad ? padded(g) : Grid{};
    const Grid& src = pad ? padded_grid : g;
    std::vector<double> w = opt.metric == Metric::l1 ? l1_bfs(src, opt.periodic) : l2_edt(src, opt.periodic);
    if (pad) w = unpad(g, w);
    for (auto& x : w)
        if (x == kInf) x = 1.0;  // no obstacle at all: a single filtration level
    return w;
}

}  // namespace p2d
