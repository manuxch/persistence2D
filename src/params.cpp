#include "params.hpp"

#include <format>
#include <fstream>
#include <set>
#include <sstream>
#include <stdexcept>
#include <thread>

#if defined(__linux__)
#include <sched.h>
#endif

namespace p2d {

namespace {

std::string trim(const std::string& s) {
    const auto b = s.find_first_not_of(" \t\r");
    if (b == std::string::npos) return {};
    const auto e = s.find_last_not_of(" \t\r");
    return s.substr(b, e - b + 1);
}

bool to_bool(const std::string& v, const std::string& where) {
    if (v == "true" || v == "1" || v == "yes") return true;
    if (v == "false" || v == "0" || v == "no") return false;
    throw std::runtime_error(std::format("{}: expected true or false, got '{}'", where, v));
}

int to_int(const std::string& v, const std::string& where) {
    try {
        std::size_t used = 0;
        const int x = std::stoi(v, &used);
        if (used == v.size()) return x;
    } catch (const std::exception&) {
    }
    throw std::runtime_error(std::format("{}: expected an integer, got '{}'", where, v));
}

}  // namespace

Params parse_params(const std::string& text, const std::string& origin) {
    Params p;
    std::istringstream in(text);
    std::string line;
    std::set<std::string> seen;
    for (int lineno = 1; std::getline(in, line); ++lineno) {
        const auto hash = line.find('#');
        const std::string content = trim(hash == std::string::npos ? line : line.substr(0, hash));
        if (content.empty()) continue;
        const std::string where = std::format("{}:{}", origin, lineno);
        const auto colon = content.find(':');
        if (colon == std::string::npos) throw std::runtime_error(std::format("{}: expected 'key: value'", where));
        const std::string key = trim(content.substr(0, colon));
        const std::string value = trim(content.substr(colon + 1));
        if (!seen.insert(key).second) throw std::runtime_error(std::format("{}: duplicate key '{}'", where, key));

        if (key == "input") {
            std::istringstream vs(value);
            for (std::string f; vs >> f;) p.inputs.push_back(f);
        } else if (key == "obstacle") {
            if (value == "dark") p.obstacle_color = ObstacleColor::dark;
            else if (value == "light") p.obstacle_color = ObstacleColor::light;
            else throw std::runtime_error(std::format("{}: obstacle must be dark or light", where));
        } else if (key == "metric") {
            if (value == "l1") p.persistence.distance.metric = Metric::l1;
            else if (value == "l2") p.persistence.distance.metric = Metric::l2;
            else throw std::runtime_error(std::format("{}: metric must be l1 or l2", where));
        } else if (key == "filtration") {
            if (value == "free") p.persistence.filtration = Filtration::free;
            else if (value == "signed") p.persistence.filtration = Filtration::signed_distance;
            else if (value == "height") p.persistence.filtration = Filtration::height;
            else if (value == "centres") p.persistence.filtration = Filtration::centres;
            else if (value == "height_exact") p.persistence.filtration = Filtration::height_exact;
            else throw std::runtime_error(std::format("{}: filtration must be free, signed, height, height_exact or centres", where));
        } else if (key == "pixel" || key == "grain" || key == "clearance" || key == "min_persistence") {
            double v;
            try {
                std::size_t used = 0;
                v = std::stod(value, &used);
                if (used != value.size()) throw std::invalid_argument("");
            } catch (const std::exception&) {
                throw std::runtime_error(std::format("{}: expected a number, got '{}'", where, value));
            }
            if (key == "clearance" || key == "min_persistence" ? v < 0 : !(v > 0))
                throw std::runtime_error(std::format("{}: {} must be {}", where, key, key == "clearance" || key == "min_persistence" ? "non-negative" : "positive"));
            (key == "pixel" ? p.continuum.pixel : key == "grain" ? p.continuum.grain : key == "clearance" ? p.continuum.clearance : p.continuum.min_persistence) = v;
        } else if (key == "periodic") {
            p.persistence.distance.periodic = to_bool(value, where);
        } else if (key == "boundary_obstacle") {
            p.persistence.distance.boundary_obstacle = to_bool(value, where);
        } else if (key == "connectivity") {
            p.persistence.connectivity = to_int(value, where);
            if (p.persistence.connectivity != 4 && p.persistence.connectivity != 8)
                throw std::runtime_error(std::format("{}: connectivity must be 4 or 8", where));
        } else if (key == "bars_file") {
            p.bars_file = value;
        } else if (key == "betti_file") {
            p.betti_file = value;
        } else if (key == "summary_file") {
            p.summary_file = value;
        } else if (key == "num_threads") {
            p.num_threads = to_int(value, where);
            if (p.num_threads < 0) throw std::runtime_error(std::format("{}: num_threads must be >= 0", where));
        } else {
            throw std::runtime_error(std::format("{}: unknown key '{}'", where, key));
        }
    }
    if (p.persistence.distance.periodic && p.persistence.distance.boundary_obstacle)
        throw std::runtime_error(std::format("{}: boundary_obstacle needs periodic: false", origin));
    return p;
}

Params read_params(const std::string& path) {
    std::ifstream in(path);
    if (!in) throw std::runtime_error(std::format("{}: cannot open", path));
    std::stringstream ss;
    ss << in.rdbuf();
    return parse_params(ss.str(), path);
}

std::string describe(const Params& p) {
    const auto& d = p.persistence.distance;
    const auto f = p.persistence.filtration;
    if (f == Filtration::height || f == Filtration::centres || f == Filtration::height_exact)
        return std::format("filtration={} pixel={} grain={} clearance={} min_persistence={}",
                           f == Filtration::height ? "height" : f == Filtration::height_exact ? "height_exact" : "centres", p.continuum.pixel, p.continuum.grain,
                           p.continuum.clearance < 0 ? std::string("auto") : std::format("{}", p.continuum.clearance),
                           p.continuum.min_persistence);
    return std::format("filtration={} metric={} periodic={} boundary_obstacle={} connectivity={} obstacle={}",
                       f == Filtration::free ? "free" : "signed",
                       d.metric == Metric::l1 ? "l1" : "l2", d.periodic, d.boundary_obstacle,
                       p.persistence.connectivity, p.obstacle_color == ObstacleColor::dark ? "dark" : "light");
}

int available_cpus() {
#if defined(__linux__)
    cpu_set_t set;
    CPU_ZERO(&set);
    if (sched_getaffinity(0, sizeof(set), &set) == 0) {
        const int n = CPU_COUNT(&set);
        if (n > 0) return n;
    }
#endif
    const auto n = static_cast<int>(std::thread::hardware_concurrency());
    return n > 0 ? n : 1;
}

}  // namespace p2d
