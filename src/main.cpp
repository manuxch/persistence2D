// persistence2D: persistent homology of the free space of 2D binary maps, filtered by corridor width.
#include "output.hpp"
#include "params.hpp"

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cstring>
#include <exception>
#include <format>
#include <iostream>
#include <stdexcept>
#include <thread>
#include <vector>

namespace {

constexpr const char* kUsage = R"(usage: persistence2D PARAMS_FILE [MAP ...]

Persistent homology of the free space of 2D binary maps (PBM/PGM, chasers .map, or text), filtered by
corridor width: the distance w from each free site to the nearest obstacle. Writes three CSV files (see the
parameter file and README.md): the persistence bars, the Betti curves and one row of descriptors per map.

Maps given on the command line replace the `input` key of the parameter file.

  -h, --help       this help
  -V, --version    program and GUDHI versions
)";

}  // namespace

int main(int argc, char** argv) {
    if (argc < 2) {
        std::cerr << kUsage;
        return 2;
    }
    if (!std::strcmp(argv[1], "-h") || !std::strcmp(argv[1], "--help")) {
        std::cout << kUsage;
        return 0;
    }
    if (!std::strcmp(argv[1], "-V") || !std::strcmp(argv[1], "--version")) {
        std::cout << "persistence2D " << P2D_VERSION << '\n';
        return 0;
    }

    try {
        p2d::Params p = p2d::read_params(argv[1]);
        if (argc > 2) p.inputs.assign(argv + 2, argv + argc);
        if (p.inputs.empty()) throw std::runtime_error("no input maps (key `input` or command line)");

        const auto n = p.inputs.size();
        int threads = p.num_threads > 0 ? p.num_threads : p2d::available_cpus();
        threads = static_cast<int>(std::min<std::size_t>(static_cast<std::size_t>(threads), n));
        const std::string options = p2d::describe(p);
        std::cout << std::format("persistence2D {}: {} map(s), {} thread(s); {}\n", P2D_VERSION, n, threads, options);

        const auto t0 = std::chrono::steady_clock::now();
        std::vector<p2d::MapResult> results(n);
        std::atomic<std::size_t> next{0};
        auto work = [&] {
            for (std::size_t i; (i = next.fetch_add(1)) < n;) {
                const std::string& file = p.inputs[i];
                try {
                    const auto f = p.persistence.filtration;
                    if (f == p2d::Filtration::height || f == p2d::Filtration::centres || f == p2d::Filtration::height_exact || f == p2d::Filtration::rips) {
                        if (!p2d::is_obs_file(file)) throw std::runtime_error(file + ": filtration height/centres/rips needs a .obs file");
                        p2d::ContinuumOptions co = p.continuum;
                        co.mode = f == p2d::Filtration::height ? p2d::ContinuumMode::height
                                       : f == p2d::Filtration::height_exact ? p2d::ContinuumMode::height_exact
                                                                            : f == p2d::Filtration::rips ? p2d::ContinuumMode::rips
                                                                            : p2d::ContinuumMode::centres;
                        results[i] = p2d::analyse_continuum(file, p2d::read_obs(file), co);
                        continue;
                    }
                    if (p2d::is_obs_file(file)) throw std::runtime_error(file + ": an .obs file needs filtration height, height_exact, centres or rips");
                    const auto fmt = p2d::detect_format(file);
                    const p2d::Grid g = p2d::read_grid(file, p.obstacle_color);
                    results[i] = p2d::analyse(file, g, p.persistence, p2d::format_name(fmt));
                } catch (const std::exception& e) {
                    results[i].file = file;
                    results[i].error = e.what();
                }
            }
        };
        {
            std::vector<std::jthread> pool;
            for (int t = 1; t < threads; ++t) pool.emplace_back(work);
            work();
        }
        const double secs = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();

        p2d::write_bars(p.bars_file, results, options);
        p2d::write_betti(p.betti_file, results, options);
        p2d::write_summary(p.summary_file, results, options);

        int failed = 0;
        for (const auto& r : results)
            if (!r.error.empty()) {
                ++failed;
                std::cerr << "error: " << r.error << '\n';
            }
        std::cout << std::format("{} map(s) in {:.2f} s", n - failed, secs);
        for (const auto* f : {&p.bars_file, &p.betti_file, &p.summary_file})
            if (!f->empty()) std::cout << ", " << *f;
        std::cout << '\n';
        if (failed) std::cerr << std::format("{} map(s) failed (listed in the summary file)\n", failed);
        return failed ? 1 : 0;
    } catch (const std::exception& e) {
        std::cerr << "error: " << e.what() << '\n';
        return 2;
    }
}
