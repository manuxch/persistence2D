// Binary obstacle maps and their readers.
#pragma once

#include <cstddef>
#include <cstdint>
#include <map>
#include <string>
#include <vector>

namespace p2d {

// A binary map on an nx x ny lattice. Site (x, y) has flat index x + nx * y; row y = 0 is the first row of
// the file (top of the image).
struct Grid {
    std::size_t nx = 0;
    std::size_t ny = 0;
    std::vector<std::uint8_t> obstacle;          // 1 = obstacle, 0 = free
    std::map<std::string, std::string> header;   // key: value pairs of a chasers .map header, if any

    std::size_t size() const { return nx * ny; }
    std::size_t index(std::size_t x, std::size_t y) const { return x + nx * y; }
    std::size_t n_obstacles() const;
    std::size_t n_free() const { return size() - n_obstacles(); }
};

// Which pixels of a PBM/PGM image are obstacles: the dark ones (black, PGM value <= maxval/2, PBM bit 1) or
// the light ones. Text and .map files have fixed symbols and ignore it.
enum class ObstacleColor { dark, light };

enum class FileFormat { pbm_plain, pbm_raw, pgm_plain, pgm_raw, chasers_map, text };

// Detects the format from the first bytes of the file.
FileFormat detect_format(const std::string& path);
const char* format_name(FileFormat f);

// Reads any supported format. Throws std::runtime_error with the file name on malformed input.
//   PBM (P1, P4): bit 1 = black.   PGM (P2, P5): value v is dark if 2 v <= maxval.
//   chasers .map: header "key: value" lines up to "---", then rows of '.' (free) and '#' (obstacle).
//   text: rows of '.'/'0' (free) and '#'/'1' (obstacle); whitespace is ignored, '%' starts a comment line.
Grid read_grid(const std::string& path, ObstacleColor obstacle_color = ObstacleColor::dark);

// Writers used by the tests and for conversions.
void write_pgm(const std::string& path, const Grid& g, ObstacleColor obstacle_color = ObstacleColor::dark);
void write_text(const std::string& path, const Grid& g);

}  // namespace p2d
