// Parameter file: "key: value" lines, '#' starts a comment.
#pragma once

#include "grid.hpp"
#include "persistence.hpp"

#include <string>
#include <vector>

namespace p2d {

struct Params {
    std::vector<std::string> inputs;   // map files (the command line overrides `input`)
    ObstacleColor obstacle_color = ObstacleColor::dark;
    PersistenceOptions persistence;    // metric, periodic, boundary_obstacle, connectivity
    std::string bars_file = "bars.csv";
    std::string betti_file = "betti.csv";
    std::string summary_file = "summary.csv";
    int num_threads = 0;               // 0 = all CPUs available to the process
};

// Reads the parameter file; unknown keys and malformed values throw std::runtime_error.
Params read_params(const std::string& path);
// The same keys from a string (used by the tests).
Params parse_params(const std::string& text, const std::string& origin = "<string>");

// One-line description of the options, for the CSV headers and the log.
std::string describe(const Params& p);

int available_cpus();

}  // namespace p2d
