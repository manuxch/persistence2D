#include "output.hpp"

#include <cmath>
#include <filesystem>
#include <format>
#include <fstream>
#include <limits>
#include <stdexcept>

namespace p2d {

namespace {

// Shortest representation that reads back exactly; NaN and infinity as numpy/pandas read them.
std::string num(double x) {
    if (std::isnan(x)) return "nan";
    if (std::isinf(x)) return x > 0 ? "inf" : "-inf";
    return std::format("{}", x);
}

// CSV field: quoted if it contains a comma, a quote or a newline.
std::string field(const std::string& s) {
    if (s.find_first_of(",\"\n") == std::string::npos) return s;
    std::string q = "\"";
    for (char c : s) q += c == '"' ? std::string("\"\"") : std::string(1, c);
    return q + '"';
}

std::ofstream open_csv(const std::string& path, const std::string& options, const char* what) {
    const auto parent = std::filesystem::path(path).parent_path();
    if (!parent.empty()) std::filesystem::create_directories(parent);
    std::ofstream out(path);
    if (!out) throw std::runtime_error(std::format("{}: cannot write", path));
    out << std::format("# persistence2D {}: {}\n# {}\n", P2D_VERSION, what, options);
    return out;
}

}  // namespace

MapResult analyse(const std::string& file, const Grid& g, const PersistenceOptions& opt, const std::string& format) {
    MapResult m;
    m.file = file;
    m.format = format;
    m.nx = g.nx;
    m.ny = g.ny;
    PersistenceResult res = compute_persistence(g, opt);
    m.levels = levels(g, res);
    m.summary = summarize(g, res);
    m.bars = std::move(res.bars);
    return m;
}

MapResult analyse_continuum(const std::string& file, const ObsFile& obs, const ContinuumOptions& opt) {
    MapResult m;
    m.file = file;
    m.format = opt.mode == ContinuumMode::height ? "obs-height" : "obs-centres";
    ContinuumResult r = compute_continuum(obs, opt);
    m.nx = r.nx;
    m.ny = r.ny;
    m.has_positions = true;
    constexpr double nan = std::numeric_limits<double>::quiet_NaN();
    Summary& s = m.summary;
    s.n_free = r.n_access;
    s.phi = 1.0 - static_cast<double>(r.n_access) / static_cast<double>(r.nx * r.ny);
    s.w_max = s.w_mean = s.r_min = s.r_wrap1 = s.r_star = s.branching = nan;
    s.h0 = summarize_dim(r.bars, 0, r.n_access);
    s.h1 = summarize_dim(r.bars, 1, r.n_access);
    m.bars = std::move(r.bars);
    return m;
}

void write_bars(const std::string& path, const std::vector<MapResult>& results, const std::string& options) {
    if (path.empty()) return;
    auto out = open_csv(path, options,
                        "persistence intervals (width: superlevel filtration of the corridor width w; height: sublevel sets of y on the accessible region; "
                        "x, y = birth cell, xs, ys = death cell, only for the continuum filtrations); "
                        "death = inf for essential classes; wraps = 1 for H1 classes winding around the torus");
    bool pos = false;
    for (const auto& m : results) pos = pos || m.has_positions;
    out << "file,dim,birth,death,persistence,essential,wraps" << (pos ? ",x,y,xs,ys" : "") << '\n';
    for (const auto& m : results) {
        if (!m.error.empty()) continue;
        const std::string f = field(m.file);
        for (const auto& b : m.bars) {
            const double inf = std::numeric_limits<double>::infinity();
            out << f << ',' << b.dim << ',' << num(b.birth) << ',' << num(b.essential ? inf : b.death) << ','
                << num(b.essential ? inf : b.persistence()) << ',' << int(b.essential) << ',' << int(b.wraps);
            if (pos) out << ',' << num(b.x) << ',' << num(b.y) << ',' << num(b.xs) << ',' << num(b.ys);
            out << '\n';
        }
    }
}

void write_betti(const std::string& path, const std::vector<MapResult>& results, const std::string& options) {
    if (path.empty()) return;
    auto out = open_csv(path, options,
                        "Betti numbers of V_r = {free x : w(x) >= r} at every distinct width r; "
                        "wraps = independent H1 classes winding around the torus");
    out << "file,r,beta0,beta1,wraps\n";
    for (const auto& m : results) {
        if (!m.error.empty()) continue;
        const std::string f = field(m.file);
        for (double r : m.levels)
            out << f << ',' << num(r) << ',' << betti_at(m.bars, 0, r) << ',' << betti_at(m.bars, 1, r) << ','
                << wraps_at(m.bars, r) << '\n';
    }
}

void write_summary(const std::string& path, const std::vector<MapResult>& results, const std::string& options) {
    if (path.empty()) return;
    auto out = open_csv(path, options, "one row per map; persistence statistics over the finite bars");
    out << "file,format,nx,ny,n_free,phi,w_max,w_mean,d_max,beta0,beta1,wrap_dim,r_wrap1,r_star,branching,"
           "h0_n,h0_essential,h0_per_free,h0_pers_median,h0_pers_q90,h0_pers_max,h0_total,h0_entropy,"
           "h1_n,h1_essential,h1_per_free,h1_pers_median,h1_pers_q90,h1_pers_max,h1_total,h1_entropy,error\n";
    for (const auto& m : results) {
        out << field(m.file) << ',' << m.format << ',' << m.nx << ',' << m.ny << ',';
        if (!m.error.empty()) {
            out << std::string(27, ',') << field(m.error) << '\n';
            continue;
        }
        const auto& s = m.summary;
        out << s.n_free << ',' << num(s.phi) << ',' << num(s.w_max) << ',' << num(s.w_mean) << ',' << num(s.d_max) << ',' << s.beta0 << ','
            << s.beta1 << ',' << s.wrap_dim << ',' << num(s.r_wrap1) << ',' << num(s.r_star) << ','
            << num(s.branching);
        for (const DimSummary* d : {&s.h0, &s.h1})
            out << ',' << d->n_finite << ',' << d->n_essential << ',' << num(d->per_free) << ','
                << num(d->pers_median) << ',' << num(d->pers_q90) << ',' << num(d->pers_max) << ','
                << num(d->total) << ',' << num(d->entropy);
        out << ",\n";
    }
}

}  // namespace p2d
