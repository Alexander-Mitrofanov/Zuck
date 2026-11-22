#ifndef PIPELINE_H
#define PIPELINE_H

#include "core/settings.h"
#include "core/graph_artifacts.h"
#include "sdbg/sdbg.h"
#include "core/sdbg_build.h"
#include "crispr_array/cycle_finder.h"
#include "crispr_array/post_processing.h"
#include <iostream>
#include <filesystem>

inline void step_build_graph(Settings& settings) {
    SDBGBuild sdbg_build(settings);
}

inline void step_load_graph(Settings& settings, SDBG& sdbg) {
    namespace fs = std::filesystem;
    std::string load_path;
    if (settings.graph_input.empty()) {
        load_path = settings.graph_folder + "/graph";
    } else {
        // If the user passed a directory, auto-append the prefix
        load_path = fs::is_directory(settings.graph_input)
            ? settings.graph_input + "/graph"
            : settings.graph_input;
    }
    std::cout << "  ▸ Loading graph: " << load_path << std::endl;
    sdbg.LoadFromFile(load_path.c_str());
    std::cout << "  ▸ Graph loaded  (" << sdbg.size() << " nodes, k=" << sdbg.k() << ")" << std::endl;
    settings.sdbg = &sdbg;
}

inline void step_find_cycles(Settings& settings) {
    CycleFinder cycle_finder(settings);
    std::cout << "  ▸ Raw cycle entries: " << cycle_finder.results.size() << " start nodes" << std::endl;
}

inline void step_post_process(Settings& settings) {
    PostProcessor processor(settings);
    processor.run_analysis();
}

inline void step_cleanup(const Settings& settings) {
    namespace fs = std::filesystem;
    std::error_code ec;
    if (!settings.retain_graph && settings.graph_input.empty()) {
        const fs::path graph_dir = settings.graph_folder;
        // Custom graph folders can contain user files. Remove only artifacts
        // belonging to this generated graph, and never recurse into directories.
        for (const auto& entry : fs::directory_iterator(graph_dir)) {
            const bool generated = Zuck::graph_artifact(entry.path().filename().string());
            if (generated && !entry.is_directory()) {
                fs::remove(entry.path(), ec);
                if (ec) std::cerr << "  Warning: could not remove graph file: " << ec.message() << std::endl;
            }
        }
        if (fs::is_empty(graph_dir)) fs::remove(graph_dir, ec);
    }
    if (!settings.retain_cycles) {
        fs::path cycles_file = fs::path(settings.output_folder) / "cycles.txt";
        fs::remove(cycles_file, ec);
        if (ec) std::cerr << "  Warning: could not remove cycles file: " << ec.message() << std::endl;
        const fs::path cycles_dir = settings.cycles_folder;
        if (fs::is_directory(cycles_dir) && fs::is_empty(cycles_dir)) fs::remove(cycles_dir, ec);
        if (ec) std::cerr << "  Warning: could not remove cycles folder: " << ec.message() << std::endl;
    }
    std::cout << "  ▸ Graph " << (settings.retain_graph || !settings.graph_input.empty() ? "retained" : "removed")
              << "; cycles " << (settings.retain_cycles ? "retained" : "removed") << std::endl;
}

#endif // PIPELINE_H
