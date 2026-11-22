#include "core/sdbg_build.h"

#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <vector>

#include "io/result_bundle.h"
#include "sequence/io/sequence_lib.h"
#include "sorting/read_to_sdbg.h"
#include "utils/utils.h"

namespace {
namespace fs = std::filesystem;

bool has_whitespace(const std::string& path) {
    return path.find_first_of(" \t\n\r\f\v") != std::string::npos;
}

// MEGAHIT's library descriptor splits paths on whitespace. Keep unique aliases
// alive while it reads the library, without copying potentially large FASTQs.
class ReadAliases {
public:
    explicit ReadAliases(std::vector<std::string> inputs) : paths(std::move(inputs)) {
        if (std::none_of(paths.begin(), paths.end(), has_whitespace)) return;
        fs::path temporary_root = fs::temp_directory_path();
        if (has_whitespace(temporary_root.string())) temporary_root = "/tmp";
        auto pattern = (temporary_root / "zuck-reads-XXXXXX").string();
        if (!mkdtemp(pattern.data()))
            throw std::runtime_error("Cannot create temporary read aliases");
        directory = pattern;
        try {
            for (std::size_t i = 0; i < paths.size(); ++i) {
                const auto alias = directory / ("read-" + std::to_string(i));
                fs::create_symlink(fs::absolute(paths[i]), alias);
                paths[i] = alias.string();
            }
        } catch (...) {
            cleanup();
            throw;
        }
    }
    ReadAliases(const ReadAliases&) = delete;
    ReadAliases& operator=(const ReadAliases&) = delete;
    ~ReadAliases() { cleanup(); }

    std::vector<std::string> paths;

private:
    fs::path directory;
    void cleanup() noexcept {
        if (directory.empty()) return;
        std::error_code ignored;
        fs::remove_all(directory, ignored);
    }
};
} // namespace

SDBGBuild::SDBGBuild(const Settings& settings) : settings(settings) {
    BuildLib();
    BuildSDBG();
}

std::string SDBGBuild::WriteLibFile(const std::vector<std::string>& inputs) const {
    namespace fs = std::filesystem;
    if (inputs.empty() || inputs.size() > 2)
        throw std::runtime_error("Graph construction requires one or two read files");
    fs::create_directories(settings.graph_folder);
    const auto path = fs::path(settings.graph_folder) / "data.lib";
    std::ofstream library(path);
    if (!library) throw std::runtime_error("Cannot create read library: " + path.string());
    library << "# Zuck read library\n" << (inputs.size() == 2 ? "pe" : "se");
    for (const auto& input : inputs) library << ' ' << input;
    library << '\n';
    result_bundle::close_checked(library, path);
    return path.string();
}

void SDBGBuild::BuildLib() const {
    const ReadAliases aliases(settings.effective_input_reads());
    const auto descriptor = WriteLibFile(aliases.paths);
    const auto output_prefix = (std::filesystem::path(settings.graph_folder) / "outfile_prefix").string();
    std::cout << "\nBuilding read library: " << descriptor << std::endl;
    AutoMaxRssRecorder recorder;
    SequenceLibCollection::Build(descriptor, output_prefix);
}

void SDBGBuild::BuildSDBG() const {
    // Settings uses GiB; MEGAHIT consumes bytes. Check before narrowing so a
    // large or non-finite request cannot wrap into an unrelated memory budget.
    const long double bytes = static_cast<long double>(settings.ram) * 1024 * 1024 * 1024;
    if (!std::isfinite(bytes) || bytes < 1 ||
        bytes > static_cast<long double>(std::numeric_limits<int64_t>::max()))
        throw std::runtime_error("Requested graph memory is outside the supported range");

    Read2SdbgOption options;
    options.k = 23;
    options.solid_threshold = 1;
    options.host_mem = static_cast<int64_t>(bytes);
    options.n_threads = settings.threads;
    options.read_lib_file = (std::filesystem::path(settings.graph_folder) / "outfile_prefix").string();
    options.output_prefix = (std::filesystem::path(settings.graph_folder) / "graph").string();
    std::cout << "  ▸ MEGAHIT memory budget: " << static_cast<int64_t>(bytes)
              << " bytes (" << settings.ram << " GiB)" << std::endl;

    AutoMaxRssRecorder recorder;
    SeqPkgWithSolidMarker package;
    {
        Read2SdbgS1 stage_one(options, &package);
        // With threshold one, initialization loads all edges without counting.
        stage_one.Initialize();
    }
    Read2SdbgS2 stage_two(options, &package);
    stage_two.Run();
}
