#include "grid.hpp"

#include <algorithm>
#include <cctype>
#include <format>
#include <fstream>
#include <iterator>
#include <sstream>
#include <stdexcept>

namespace p2d {

std::size_t Grid::n_obstacles() const {
    return static_cast<std::size_t>(std::count(obstacle.begin(), obstacle.end(), std::uint8_t{1}));
}

namespace {

std::string slurp(const std::string& path) {
    std::ifstream in(path, std::ios::binary);
    if (!in) throw std::runtime_error(std::format("{}: cannot open", path));
    return {std::istreambuf_iterator<char>(in), std::istreambuf_iterator<char>()};
}

[[noreturn]] void fail(const std::string& path, const std::string& what) {
    throw std::runtime_error(std::format("{}: {}", path, what));
}

// Cursor over the header of a PNM file: whitespace-separated tokens, '#' comments up to the end of line.
struct PnmCursor {
    const std::string& s;
    std::size_t pos = 2;  // after the magic number

    void skip() {
        while (pos < s.size()) {
            if (std::isspace(static_cast<unsigned char>(s[pos]))) ++pos;
            else if (s[pos] == '#') { while (pos < s.size() && s[pos] != '\n') ++pos; }
            else break;
        }
    }
    unsigned long number(const std::string& path) {
        skip();
        std::size_t start = pos;
        while (pos < s.size() && std::isdigit(static_cast<unsigned char>(s[pos]))) ++pos;
        if (start == pos) fail(path, "malformed PNM header or data");
        return std::stoul(s.substr(start, pos - start));
    }
};

Grid read_pnm(const std::string& path, const std::string& s, FileFormat fmt, ObstacleColor color) {
    PnmCursor c{s};
    Grid g;
    g.nx = c.number(path);
    g.ny = c.number(path);
    if (g.nx == 0 || g.ny == 0) fail(path, "empty image");
    const bool pbm = fmt == FileFormat::pbm_plain || fmt == FileFormat::pbm_raw;
    const unsigned long maxval = pbm ? 1 : c.number(path);
    if (maxval == 0 || maxval > 65535) fail(path, std::format("invalid maxval {}", maxval));
    g.obstacle.assign(g.size(), 0);

    auto set = [&](std::size_t i, bool dark) {
        g.obstacle[i] = (dark == (color == ObstacleColor::dark)) ? 1 : 0;
    };
    const std::size_t n = g.size();
    switch (fmt) {
        case FileFormat::pbm_plain:
            for (std::size_t i = 0; i < n; ++i) {
                c.skip();  // P1 digits may be packed without spaces
                if (c.pos >= s.size() || (s[c.pos] != '0' && s[c.pos] != '1')) fail(path, "truncated P1 data");
                set(i, s[c.pos++] == '1');
            }
            break;
        case FileFormat::pgm_plain:
            for (std::size_t i = 0; i < n; ++i) {
                const auto v = c.number(path);
                if (v > maxval) fail(path, std::format("value {} above maxval {}", v, maxval));
                set(i, 2 * v <= maxval);
            }
            break;
        case FileFormat::pbm_raw: {
            ++c.pos;  // exactly one whitespace byte after the header
            const std::size_t row_bytes = (g.nx + 7) / 8;
            if (s.size() < c.pos + row_bytes * g.ny) fail(path, "truncated P4 data");
            for (std::size_t y = 0; y < g.ny; ++y)
                for (std::size_t x = 0; x < g.nx; ++x) {
                    const auto byte = static_cast<unsigned char>(s[c.pos + y * row_bytes + x / 8]);
                    set(g.index(x, y), (byte >> (7 - x % 8)) & 1u);
                }
            break;
        }
        case FileFormat::pgm_raw: {
            ++c.pos;
            const std::size_t bpp = maxval < 256 ? 1 : 2;
            if (s.size() < c.pos + bpp * n) fail(path, "truncated P5 data");
            for (std::size_t i = 0; i < n; ++i) {
                unsigned long v = static_cast<unsigned char>(s[c.pos + bpp * i]);
                if (bpp == 2) v = (v << 8) | static_cast<unsigned char>(s[c.pos + 2 * i + 1]);
                set(i, 2 * v <= maxval);
            }
            break;
        }
        default:
            fail(path, "internal error: not a PNM format");
    }
    return g;
}

// Rows of map symbols; used by the text format and by the body of chasers .map files.
void parse_rows(const std::string& path, std::istream& in, Grid& g) {
    std::string line;
    std::vector<std::uint8_t> cells;
    while (std::getline(in, line)) {
        if (!line.empty() && line[0] == '%') continue;
        std::size_t width = 0;
        for (char ch : line) {
            if (std::isspace(static_cast<unsigned char>(ch))) continue;
            if (ch == '.' || ch == '0') cells.push_back(0);
            else if (ch == '#' || ch == '1') cells.push_back(1);
            else fail(path, std::format("unexpected character '{}' in row {}", ch, g.ny + 1));
            ++width;
        }
        if (width == 0) continue;
        if (g.nx == 0) g.nx = width;
        else if (width != g.nx) fail(path, std::format("row {} has {} sites, expected {}", g.ny + 1, width, g.nx));
        ++g.ny;
    }
    if (g.nx == 0) fail(path, "no map rows");
    g.obstacle = std::move(cells);
}

std::string trim(const std::string& s) {
    const auto b = s.find_first_not_of(" \t\r");
    if (b == std::string::npos) return {};
    const auto e = s.find_last_not_of(" \t\r");
    return s.substr(b, e - b + 1);
}

Grid read_chasers_map(const std::string& path, const std::string& s) {
    std::istringstream in(s);
    Grid g;
    std::string line;
    bool body = false;
    while (std::getline(in, line)) {
        if (trim(line) == "---") { body = true; break; }
        const auto hash = line.find('#');
        const std::string content = trim(hash == std::string::npos ? line : line.substr(0, hash));
        if (content.empty()) continue;
        const auto colon = content.find(':');
        if (colon == std::string::npos) fail(path, std::format("malformed header line '{}'", line));
        g.header[trim(content.substr(0, colon))] = trim(content.substr(colon + 1));
    }
    if (!body) fail(path, "no '---' line after the header");
    parse_rows(path, in, g);
    if (auto it = g.header.find("L"); it != g.header.end()) {
        const auto L = std::stoul(it->second);
        if (g.nx != L || g.ny != L) fail(path, std::format("header L = {} but the map is {} x {}", L, g.nx, g.ny));
    }
    if (auto it = g.header.find("n_obs"); it != g.header.end()) {
        if (std::stoul(it->second) != g.n_obstacles())
            fail(path, std::format("header n_obs = {} but the map has {}", it->second, g.n_obstacles()));
    }
    return g;
}

}  // namespace

FileFormat detect_format(const std::string& path) {
    std::ifstream in(path, std::ios::binary);
    if (!in) throw std::runtime_error(std::format("{}: cannot open", path));
    char head[16] = {};
    in.read(head, sizeof head);
    const std::string h(head, static_cast<std::size_t>(in.gcount()));
    if (h.size() >= 2 && h[0] == 'P') {
        switch (h[1]) {
            case '1': return FileFormat::pbm_plain;
            case '2': return FileFormat::pgm_plain;
            case '4': return FileFormat::pbm_raw;
            case '5': return FileFormat::pgm_raw;
            default: break;
        }
    }
    if (h.starts_with("# chasers-map")) return FileFormat::chasers_map;
    return FileFormat::text;
}

const char* format_name(FileFormat f) {
    switch (f) {
        case FileFormat::pbm_plain: return "pbm-P1";
        case FileFormat::pbm_raw: return "pbm-P4";
        case FileFormat::pgm_plain: return "pgm-P2";
        case FileFormat::pgm_raw: return "pgm-P5";
        case FileFormat::chasers_map: return "chasers-map";
        case FileFormat::text: return "text";
    }
    return "?";
}

Grid read_grid(const std::string& path, ObstacleColor obstacle_color) {
    const FileFormat fmt = detect_format(path);
    const std::string s = slurp(path);
    switch (fmt) {
        case FileFormat::chasers_map: return read_chasers_map(path, s);
        case FileFormat::text: {
            std::istringstream in(s);
            Grid g;
            parse_rows(path, in, g);
            return g;
        }
        default: return read_pnm(path, s, fmt, obstacle_color);
    }
}

void write_pgm(const std::string& path, const Grid& g, ObstacleColor obstacle_color) {
    std::ofstream out(path);
    if (!out) throw std::runtime_error(std::format("{}: cannot write", path));
    out << "P2\n" << g.nx << ' ' << g.ny << "\n1\n";
    const int obstacle_value = obstacle_color == ObstacleColor::dark ? 0 : 1;
    for (std::size_t y = 0; y < g.ny; ++y)
        for (std::size_t x = 0; x < g.nx; ++x)
            out << (g.obstacle[g.index(x, y)] ? obstacle_value : 1 - obstacle_value) << (x + 1 < g.nx ? ' ' : '\n');
}

void write_text(const std::string& path, const Grid& g) {
    std::ofstream out(path);
    if (!out) throw std::runtime_error(std::format("{}: cannot write", path));
    for (std::size_t y = 0; y < g.ny; ++y) {
        for (std::size_t x = 0; x < g.nx; ++x) out << (g.obstacle[g.index(x, y)] ? '#' : '.');
        out << '\n';
    }
}

}  // namespace p2d
