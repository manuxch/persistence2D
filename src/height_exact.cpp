// Exact height persistence of the region accessible to a grain, without a pixel grid.
//
// A = {x : |x - c_i| >= rho for all i}, y up, gravity toward -y, x periodic. The sublevel sets {y <= h} of A change
// only at two kinds of boundary points:
//   - minima: a corner where two circles meet with the region above both tangents (a V trap, the rest point of
//     a grain on two obstacles);
//   - merges: the top of a disc, if no other disc covers it. Before it, the sublevel set has a piece on each side
//     of the disc; at its height the two pieces meet (or, if they were already one, a loop closes).
// The piece on each side is identified by letting a grain slide down: along the upper arc of the disc, onto the
// next circle at a corner, or in free fall from an equator, until it rests at a V trap or leaves through the bottom
// (the "ocean", the unbounded component below every obstacle). Saddles are then processed in increasing height with
// union-find and the elder rule, which gives the H0 pairs exactly. H1 (loops around obstacles, all essential) is not
// reported.
#include "continuum.hpp"

#include <algorithm>
#include <cmath>
#include <map>
#include <numeric>
#include <stdexcept>

namespace p2d {

namespace {

constexpr double kPi = 3.14159265358979323846;
constexpr double kEps = 1e-12;

struct Trap {
    double x, y;
    int i, j;
};

class Geometry {
public:
    Geometry(const ObsFile& o, double rho) : Lx_(o.Lx), rho_(rho) {
        const std::size_t n = o.centres.size();
        x_.resize(n);
        y_.resize(n);
        for (std::size_t i = 0; i < n; ++i) {
            x_[i] = o.centres[i][0] - Lx_ * std::floor(o.centres[i][0] / Lx_);
            y_[i] = o.centres[i][1];
        }
        order_.resize(n);
        std::iota(order_.begin(), order_.end(), 0);
        std::sort(order_.begin(), order_.end(), [&](int a, int b) { return x_[a] < x_[b]; });
        sorted_x_.resize(n);
        for (std::size_t k = 0; k < n; ++k) sorted_x_[k] = x_[order_[k]];
        // Neighbours within 2 rho, from the x-sorted order.
        nb_.assign(n, {});
        for (std::size_t i = 0; i < n; ++i)
            column(x_[i], 2 * rho_, [&](int k) {
                if (k == static_cast<int>(i)) return;
                const double dx = dx_(x_[k] - x_[i]), dy = y_[k] - y_[i];
                if (dx * dx + dy * dy < 4 * rho_ * rho_) nb_[i].push_back(k);
            });
    }

    std::size_t size() const { return x_.size(); }
    double dx_(double d) const { return d - Lx_ * std::round(d / Lx_); }
    double reduce(double x) const { return x - Lx_ * std::floor(x / Lx_); }
    double X(int i) const { return x_[i]; }
    double Y(int i) const { return y_[i]; }
    const std::vector<int>& near(int i) const { return nb_[i]; }
    double rho() const { return rho_; }

    // Calls f(k) for every disc whose x is within `half` of x (periodic).
    template <class F>
    void column(double x, double half, F&& f) const {
        const double w = std::min(half, 0.5 * Lx_);
        auto scan = [&](double a, double b) {
            auto lo = std::lower_bound(sorted_x_.begin(), sorted_x_.end(), a);
            for (auto it = lo; it != sorted_x_.end() && *it <= b; ++it) f(order_[it - sorted_x_.begin()]);
        };
        if (2 * w >= Lx_) { scan(-1, Lx_ + 1); return; }
        const double a = x - w, b = x + w;
        if (a >= 0 && b < Lx_) scan(a, b);
        else if (a < 0) { scan(0, b); scan(a + Lx_, Lx_); }
        else { scan(a, Lx_); scan(0, b - Lx_); }
    }

    bool covered(double px, double py, int skip1, int skip2) const {
        bool c = false;
        column(px, rho_, [&](int k) {
            if (k == skip1 || k == skip2) return;
            const double dx = dx_(px - x_[k]), dy = py - y_[k];
            if (dx * dx + dy * dy < rho_ * rho_ * (1 - kEps)) c = true;
        });
        return c;
    }

private:
    double Lx_, rho_;
    std::vector<double> x_, y_, sorted_x_;
    std::vector<int> order_;
    std::vector<std::vector<int>> nb_;
};

struct Counters {
    long underside = 0, missing_trap = 0, stuck = 0;
};

class Solver {
public:
    Solver(const ObsFile& o, double rho) : g_(o, rho), Lx_(o.Lx) {}

    std::vector<Bar> run(Counters& cnt) {
        enumerate_traps();
        const int n = static_cast<int>(g_.size());
        struct Saddle { double h, x; int a, b; };
        std::vector<Saddle> saddles;
        for (int c = 0; c < n; ++c) {
            const double top = g_.Y(c) + g_.rho();
            if (g_.covered(g_.X(c), top, c, -1)) continue;
            const int l = slide(c, 0.5 * kPi, -1, cnt), r = slide(c, 0.5 * kPi, +1, cnt);
            saddles.push_back({top, g_.X(c), l, r});
        }
        std::sort(saddles.begin(), saddles.end(), [](const Saddle& a, const Saddle& b) { return a.h < b.h; });

        const int ocean = static_cast<int>(traps_.size());
        std::vector<int> parent(ocean + 1);
        std::iota(parent.begin(), parent.end(), 0);
        std::vector<double> birth(ocean + 1);
        for (int t = 0; t < ocean; ++t) birth[t] = traps_[t].y;
        birth[ocean] = -std::numeric_limits<double>::infinity();
        auto find = [&](int a) {
            while (parent[a] != a) a = parent[a] = parent[parent[a]];
            return a;
        };
        std::vector<Bar> bars;
        for (const auto& s : saddles) {
            int a = find(s.a), b = find(s.b);
            if (a == b) continue;                       // a loop closes, not a merge
            if (birth[a] > birth[b]) std::swap(a, b);   // a is the elder
            const int dead = b;
            parent[b] = a;
            Bar bar;
            bar.dim = 0;
            bar.birth = birth[dead];
            bar.death = s.h;
            bar.x = traps_[dead].x;
            bar.y = traps_[dead].y;
            bar.xs = s.x;
            bar.ys = s.h;
            bars.push_back(bar);
        }
        for (int t = 0; t < ocean; ++t)
            if (find(t) == t && find(ocean) != t) {     // a closed void: its lowest trap never joins another basin
                Bar bar;
                bar.dim = 0;
                bar.essential = true;
                bar.birth = traps_[t].y;
                bar.x = traps_[t].x;
                bar.y = traps_[t].y;
                bars.push_back(bar);
            }
        return bars;
    }

private:
    Geometry g_;
    double Lx_;
    std::vector<Trap> traps_;
    std::map<std::pair<int, int>, int> trap_of_;

    // The corner of circles i and j that is the upper intersection point, if it is a V trap.
    bool corner_trap(int i, int j, double& gx, double& gy) const {
        const double dx = g_.dx_(g_.X(j) - g_.X(i)), dy = g_.Y(j) - g_.Y(i);
        const double D = std::hypot(dx, dy);
        const double rho = g_.rho();
        if (!(D > 0) || !(D < 2 * rho)) return false;
        double nx = -dy / D, ny = dx / D;
        if (ny < 0) nx = -nx, ny = -ny;
        if (!(ny > 0)) return false;
        const double h = std::sqrt(rho * rho - 0.25 * D * D);
        const double px = 0.5 * dx + h * nx, py = 0.5 * dy + h * ny;       // relative to disc i
        if (!(px > std::min(0.0, dx) && px < std::max(0.0, dx))) return false;
        for (int m : g_.near(i)) {
            if (m == j) continue;
            const double mx = g_.dx_(g_.X(m) - g_.X(i)) - px, my = g_.Y(m) - g_.Y(i) - py;
            if (mx * mx + my * my < rho * rho * (1 - kEps)) return false;
        }
        gx = g_.reduce(g_.X(i) + px);
        gy = g_.Y(i) + py;
        return true;
    }

    void enumerate_traps() {
        for (int i = 0; i < static_cast<int>(g_.size()); ++i)
            for (int j : g_.near(i)) {
                if (j < i) continue;
                double gx, gy;
                if (!corner_trap(i, j, gx, gy)) continue;
                trap_of_[{i, j}] = static_cast<int>(traps_.size());
                traps_.push_back({gx, gy, i, j});
            }
    }

    int trap_id(int i, int j, double gx, double gy, Counters& cnt) {
        const auto key = std::make_pair(std::min(i, j), std::max(i, j));
        const auto it = trap_of_.find(key);
        if (it != trap_of_.end()) return it->second;
        ++cnt.missing_trap;   // a corner reached by sliding that the enumeration did not list (numerical tie)
        trap_of_[key] = static_cast<int>(traps_.size());
        traps_.push_back({g_.reduce(gx), gy, key.first, key.second});
        return static_cast<int>(traps_.size()) - 1;
    }

    // A grain released at (px, py) falls vertically onto the highest upper surface below it, then slides.
    int fall(double px, double py, Counters& cnt) {
        int best = -1;
        double best_top = -std::numeric_limits<double>::infinity(), best_dx = 0;
        const double rho = g_.rho();
        g_.column(px, rho, [&](int k) {
            const double dx = g_.dx_(px - g_.X(k));
            if (std::abs(dx) >= rho) return;
            const double top = g_.Y(k) + std::sqrt(rho * rho - dx * dx);
            if (top < py - 1e-10 && top > best_top) { best = k; best_top = top; best_dx = dx; }
        });
        if (best < 0) return static_cast<int>(traps_.size());   // the ocean
        const double theta = std::atan2(best_top - g_.Y(best), best_dx);
        return slide(best, theta, best_dx > 0 ? +1 : -1, cnt);
    }

    // Slide down the upper arc of disc c from angle theta0 (the arc is theta in [0, pi], theta = pi/2 at the top),
    // to the right (s = +1, theta decreasing) or to the left (s = -1, theta increasing).
    int slide(int c, double theta0, int s, Counters& cnt) {
        const double rho = g_.rho();
        int prev = -1;
        for (int step = 0; step < 100000; ++step) {
            const double theta_end = s > 0 ? 0.0 : kPi;
            const double travel = std::abs(theta_end - theta0);
            double best = travel;
            int hit = -1;
            for (int k : g_.near(c)) {
                if (k == prev) continue;
                const double dx = g_.dx_(g_.X(k) - g_.X(c)), dy = g_.Y(k) - g_.Y(c);
                const double D = std::hypot(dx, dy), phi = std::atan2(dy, dx), alpha = std::acos(D / (2 * rho));
                // The arc of c inside disc k is (phi - alpha, phi + alpha); the grain enters at the end it reaches first.
                double delta = s > 0 ? theta0 - (phi + alpha) : (phi - alpha) - theta0;
                delta -= 2 * kPi * std::floor(delta / (2 * kPi));
                if (delta < best) { best = delta; hit = k; }
            }
            if (hit < 0) {   // no corner: leave the arc at the equator and fall
                return fall(g_.X(c) + s * rho, g_.Y(c), cnt);
            }
            const double th = theta0 - s * best;
            const double ncx = std::cos(th), ncy = std::sin(th);
            const double gx = g_.X(c) + rho * ncx, gy = g_.Y(c) + rho * ncy;
            const double kx = g_.dx_(g_.X(hit) - g_.X(c)) - rho * ncx, ky = g_.Y(hit) - g_.Y(c) - rho * ncy;
            const double nkx = -kx / rho, nky = -ky / rho;   // outward normal of disc `hit` at the corner
            double tx = -nky, ty = nkx;
            if (tx * ncx + ty * ncy < 0) tx = -tx, ty = -ty;   // along `hit`, away from disc c
            if (ty >= 0) return trap_id(c, hit, gx, gy, cnt);   // the boundary rises on both sides: a V trap
            if (nky > 0) {   // keep sliding on the upper half of `hit`
                prev = c;
                c = hit;
                theta0 = std::atan2(nky, nkx);
                s = tx > 0 ? +1 : -1;
                continue;
            }
            ++cnt.underside;   // the next surface faces down: the grain falls from the corner
            return fall(gx, gy, cnt);
        }
        ++cnt.stuck;
        return static_cast<int>(traps_.size());
    }
};

}  // namespace

std::vector<Bar> height_exact(const ObsFile& obs, double rho) {
    if (!obs.periodic_x) throw std::runtime_error("height_exact needs a periodic_x arrangement");
    if (!(rho > 0)) throw std::invalid_argument("height_exact: rho must be positive");
    if (obs.centres.empty()) return {};
    Counters cnt;
    Solver s(obs, rho);
    auto bars = s.run(cnt);
    if (cnt.stuck) throw std::runtime_error("height_exact: a sliding path did not terminate");
    return bars;
}

}  // namespace p2d
