// CSV outputs. Every file starts with '#' comment lines (program version and options), then a header row.
#pragma once

#include "continuum.hpp"
#include "descriptors.hpp"
#include "grid.hpp"
#include "persistence.hpp"

#include <string>
#include <vector>

namespace p2d {

// Everything kept from one map (the distance field is not: it is large and can be recomputed).
struct MapResult {
    std::string file;
    std::string format;
    std::size_t nx = 0;
    std::size_t ny = 0;
    std::vector<Bar> bars;
    std::vector<double> levels;  // decreasing
    Summary summary;
    bool has_positions = false;  // bars carry the (x, y) of their cells (continuum filtrations)
    std::string error;           // non-empty if the map could not be processed
};

MapResult analyse(const std::string& file, const Grid& g, const PersistenceOptions& opt, const std::string& format);

// Continuum filtrations (height, centres) of a .obs arrangement.
MapResult analyse_continuum(const std::string& file, const ObsFile& obs, const ContinuumOptions& opt);

// Empty path = do not write that file.
void write_bars(const std::string& path, const std::vector<MapResult>& results, const std::string& options);
void write_betti(const std::string& path, const std::vector<MapResult>& results, const std::string& options);
void write_summary(const std::string& path, const std::vector<MapResult>& results, const std::string& options);

}  // namespace p2d
