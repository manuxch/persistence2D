// Continuum obstacle arrangements (hyperGalton2D .obs files) and two filtrations of the region a grain can
// reach, computed on a fine pixel grid with distances taken to the exact obstacle centres.
#pragma once

#include "persistence.hpp"

#include <array>
#include <string>
#include <vector>

namespace p2d {

// A hyperGalton2D obstacle file (format 1): `# key: value` header, then one `x y` centre per line.
struct ObsFile {
    double Lx = 0.0, Ly = 0.0;   // box; the arrangement is periodic in x if `periodic_x`
    double R = 0.0;              // obstacle radius
    bool periodic_x = true;      // boundary: periodic_x (true) or walls_x (false)
    std::string sha256_body;
    std::vector<std::array<double, 2>> centres;
};

bool is_obs_file(const std::string& path);      // first line is `# format: 1`
ObsFile read_obs(const std::string& path);      // throws std::runtime_error with the file name

// height : sublevel sets {y <= h} of the accessible region A = {x : dist(x, centres) >= clearance} (and, with
//          walls_x, at least grain/2 from the two walls). Gravity points to -y, so H0 births are the local
//          minima of A (the V traps: a grain at rest on two obstacles) and the persistence of each is the
//          height the grain must climb to reach a lower basin. clearance defaults to R + grain/2.
// centres: the superlevel filtration of the distance to the nearest centre, V_r = {x : dist >= r}; H0 and H1
//          of V_r are dual to H1 and H0 of the alpha filtration of the centres (bars are in radius units).
enum class ContinuumMode { height, centres, height_exact };

struct ContinuumOptions {
    ContinuumMode mode = ContinuumMode::height;
    double pixel = 0.05;         // target pixel side, in the length unit of the file (the grid fits the box)
    double grain = 1.0;          // grain diameter
    double clearance = -1.0;     // < 0: R + grain/2 (height) or 0 (centres)
    // Bars shorter than this are dropped (essential bars are kept). 0 keeps everything; a couple of pixels removes
    // the discretisation noise of the distance field (`centres`), but also the shallowest traps (`height`).
    double min_persistence = 0.0;
};

struct ContinuumResult {
    std::size_t nx = 0, ny = 0;
    std::size_t n_access = 0;    // pixels in the filtration
    double clearance = 0.0;
    double min_persistence = 0.0;   // the threshold used
    std::vector<Bar> bars;       // with positions (x, y) of the birth cell and (xs, ys) of the death cell
};

// Exact version of `height`, with no pixel grid: the V traps are the corners of the inflated discs and the
// barriers are the tops of discs, combined by union-find in order of height (see src/height_exact.cpp). rho is the
// clearance, R + grain/2. Needs a periodic_x arrangement. Returns the H0 bars (x, y: trap; xs, ys: top of the
// disc where the basin joins a lower one); the unbounded component below every obstacle has no bar, and closed
// voids give essential bars. H1 is not reported.
std::vector<Bar> height_exact(const ObsFile& obs, double rho);

ContinuumResult compute_continuum(const ObsFile& obs, const ContinuumOptions& opt);

// Distance from each pixel centre to the nearest obstacle centre (minimum image in x if periodic). Pixel (i, j)
// is centred at ((i + 1/2) Lx / nx, (j + 1/2) Ly / ny); x runs fastest. Exposed for the tests.
std::vector<double> centre_distance(const ObsFile& obs, std::size_t nx, std::size_t ny);

}  // namespace p2d
