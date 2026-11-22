#include "core/cli.h"
#include "core/build_identity.h"
#include "core/settings.h"
#include "core/graph_artifacts.h"
#include "io/result_bundle.h"

#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <optional>
#include <set>

#ifdef _WIN32
#include <windows.h>
#elif defined(__APPLE__)
#include <sys/sysctl.h>
#endif

namespace {

constexpr char banner[] = "\nZuck " ZUCK_VERSION " — Metagenomic CRISPR Array Analysis\n\n";

std::optional<unsigned long long> read_memory_value(const fs::path& path) {
    std::ifstream input(path);
    std::string value;
    if (!(input >> value) || value == "max") return std::nullopt;
    try {
        const auto parsed = std::stoull(value);
        // cgroup v1 commonly uses a huge sentinel for an unlimited controller.
        if (parsed >= (1ULL << 60)) return std::nullopt;
        return parsed;
    } catch (...) {
        return std::nullopt;
    }
}

#ifdef __linux__
std::optional<unsigned long long> linux_mem_available() {
    std::ifstream meminfo("/proc/meminfo");
    std::string key;
    unsigned long long kib = 0;
    std::string unit;
    while (meminfo >> key >> kib >> unit) {
        if (key == "MemAvailable:") return kib * 1024ULL;
    }
    return std::nullopt;
}

std::optional<fs::path> process_cgroup_directory() {
    std::ifstream cgroups("/proc/self/cgroup");
    std::string line;
    while (std::getline(cgroups, line)) {
        const auto first = line.find(':');
        const auto second = line.find(':', first == std::string::npos
                                                ? first : first + 1);
        if (first == std::string::npos || second == std::string::npos) continue;
        const std::string controllers = line.substr(first + 1, second - first - 1);
        std::string relative = line.substr(second + 1);
        while (!relative.empty() && relative.front() == '/') relative.erase(0, 1);

        if (controllers.empty())
            return fs::path("/sys/fs/cgroup") / relative;
        if (controllers == "memory" ||
            controllers.find(",memory") != std::string::npos ||
            controllers.find("memory,") != std::string::npos)
            return fs::path("/sys/fs/cgroup/memory") / relative;
    }
    return std::nullopt;
}

std::optional<unsigned long long> cgroup_memory_available() {
    const auto directory = process_cgroup_directory();
    if (!directory) return std::nullopt;

    const bool v2 = fs::exists(*directory / "memory.max");
    const fs::path controller_root =
        v2 ? fs::path("/sys/fs/cgroup") : fs::path("/sys/fs/cgroup/memory");
    std::optional<unsigned long long> smallest_remaining;
    fs::path current = *directory;
    while (current.string().size() >= controller_root.string().size()) {
        const fs::path max_path = current /
            (v2 ? "memory.max" : "memory.limit_in_bytes");
        const fs::path used_path = current /
            (v2 ? "memory.current" : "memory.usage_in_bytes");
        const auto limit = read_memory_value(max_path);
        const auto used = read_memory_value(used_path);
        if (limit && used) {
            const unsigned long long remaining =
                *used >= *limit ? 0ULL : *limit - *used;
            if (!smallest_remaining || remaining < *smallest_remaining)
                smallest_remaining = remaining;
        }
        if (current == controller_root || current == current.parent_path()) break;
        current = current.parent_path();
    }
    return smallest_remaining;
}
#endif

}  // namespace

// Return the memory currently available to this process in GiB. A negative
// result means the platform did not provide a usable measurement; zero is a
// real, detected no-memory condition and must not trigger an unsafe fallback.
double get_available_system_ram() {
    std::optional<unsigned long long> available_bytes;

#ifdef __linux__
    available_bytes = linux_mem_available();
    const auto cgroup_available = cgroup_memory_available();
    if (cgroup_available &&
        (!available_bytes || *cgroup_available < *available_bytes))
        available_bytes = cgroup_available;
#elif defined(_WIN32)
    MEMORYSTATUSEX mem_info;
    mem_info.dwLength = sizeof(MEMORYSTATUSEX);
    if (GlobalMemoryStatusEx(&mem_info)) {
        available_bytes = mem_info.ullAvailPhys;
    }
#elif defined(__APPLE__)
    int mib[2] = {CTL_HW, HW_MEMSIZE};
    int64_t mem_size;
    size_t len = sizeof(mem_size);
    if (sysctl(mib, 2, &mem_size, &len, nullptr, 0) == 0) {
        available_bytes = static_cast<unsigned long long>(mem_size);
    }
#endif
    if (!available_bytes) return -1.0;
    return static_cast<double>(*available_bytes) /
           (1024.0 * 1024.0 * 1024.0);
}

void print_help() {
    cout << banner
                 << "Usage:\n"
                 << "  Zuck --input-files <file1> [file2] [options]\n"
                 << "  Zuck --graph <path> [options]\n"
                 << "  Zuck --cycles <cycles.txt> [options]\n"
                 << "\nInput (one required):\n"
                 << "  --input-files <file1> [file2]   One or two FASTA/FASTQ files (plain or gzipped)\n"
                 << "  --graph <path>                  Pre-built SDBG graph folder (skips construction)\n"
                 << "  --cycles <path>                 Retained cycles.txt (post-process only)\n"
                 << "\nOptions:\n"
                 << "  --output-folder <path>          Output directory  [default: Zuck_run_TIMESTAMP]\n"
                 << "  --ram <amount>                  Graph budget: B/K/M/G suffix  [default: 85% available]\n"
                 << "  --threads <num>                 Thread count  [default: CPU cores - 2]\n"
                 << "  --cycle-max-length <int>        Maximum cycle length  [default: 90]\n"
                 << "  --cycle-min-length <int>        Minimum cycle length  [default: 27]\n"
                 << "  --threshold-multiplicity <int>  Edge multiplicity cutoff  [default: 20]\n"
                 << "  --low-abundance <true|false>    Retain singleton graph edges [default: false]\n"
                 << "  --adaptive-coverage <true|false> Use effective graph edge multiplicity to adapt cycle discovery [default: false]\n"
                 << "                                   Explicit threshold-multiplicity and low-abundance settings override their adaptive values.\n"
                 << "  --read-first <true|false>       Discover repeat families in reads and export the combined spacer inventory [default: false]\n"
                 << "  --read-first-family-recut <true|false>  Reassign graph boundaries supported by distinct read repeat families [default: false]\n"
                 << "  --traversal-partitions <int>   Logical search partitions, independent of threads [default: 8]\n"
                 << "  --spacer-min-length <int>      Optional output length bound; 0 disables [default: 0]\n"
                 << "  --spacer-max-length <int>      Optional output length bound; 0 disables [default: 0]\n"
                 << "  --repeat-min-length <int>      Optional output length bound; 0 disables [default: 0]\n"
                 << "  --repeat-max-length <int>      Optional output length bound; 0 disables [default: 0]\n"
                 << "  --boundary-mode <legacy|two-sided|branch-aware>\n"
                 << "                                   Repeat boundary inference [default: two-sided]\n"
                 << "                                   branch-aware is experimental\n"
                 << "  --repeat-mixture-refinement <bool> Experimental nested branch refinement (default false)\n"
                 << "                                   Requires --boundary-mode branch-aware\n"
                 << "  --short-spacer-filter <true|false>\n"
                 << "                                   Reject arrays explained by one short spacer [default: false]\n"
                 << "                                   Applies to two-sided and branch-aware modes\n"
                 << "  --terminal-variant-rescue <true|false>\n"
                 << "                                   Retain equally supported terminal variants [default: true]\n"
                 << "  --read-linked-groups <true|false>\n"
                 << "                                   Report read-supported subgroups [default: true]\n"
                 << "  --grouping-reads <file1> [file2]\n"
                 << "                                   Reads for subgroup evidence in cycle/graph replay\n"
                 << "                                   Full runs use input reads automatically\n"
                 << "  --array-database <file.tsv>     Add confidence-weighted catalog-supported cycle groups\n"
                 << "  --database-read-validation <bool> Require exact linear FASTQ support (default true)\n"
                 << "  --database-novel-spacers <bool> Experimental observed novel-spacer rescue (default false)\n"
                 << "  --database-groups <bool>       Catalog-coherent content views (default true)\n"
                 << "  --database-include-possible <bool> Include Possible-tier evidence (default false)\n"
                 << "                                   Optional; requires two-sided boundaries\n"
                 << "  --confidence-scores <bool>     Annotate repeat, spacer and boundary support (default true)\n"
                 << "  --call-profile <standard|context|off>\n"
                 << "                                   Additional supported calls; primary evidence retained\n"
                 << "                                   standard: read context + 26-50 nt; context: any length\n"
                 << "                                   Default: standard; off omits supported-call files\n"
                 << "  --family-repeat-context <bool> Experimental next-repeat allele context (default false)\n"
                 << "  --read-recruitment <bool>     Experimental graph-seeded read recovery (default false)\n"
                 << "  --local-boundary-graph <bool> Experimental read-local boundary repair (default false)\n"
                 << "  --motif-cliff-boundary <bool> Experimental replacement-only motif repair (default false)\n"
                 << "  --singleton-precision-prune <bool> Experimental supported-view demotion (default false)\n"
                 << "  --circular-unit-dedup <bool> One supported partition per circular DNA unit (default false)\n"
                 << "  --directed-primary-groups <bool> Experimental disjoint directed grouping (default false)\n"
                 << "  --candidate-evidence <bool>    Score all exact pre-collapse de novo pairs (default false)\n"
                 << "  --boundary-radius <0..3>      Independent shifts at both junctions (default 3)\n"
                 << "  --boundary-policy <mode>      off|catalogue|catalogue-all|catalogue-best|anchored|anchored-best\n"
                 << "                                Add read-supported cuts using database anchors (default catalogue)\n"
                 << "  --experimental-alternative-views\n"
                 << "                                   Emit bounded alternative-repeat sidecar [off]\n"
                 << "  --settings <path>               Key=value settings file (CLI overrides)\n"
                 << "  --retain-graph                  Keep the large graph after a full run\n"
                 << "  --cleanup-cycles                Remove cycles.txt after a full run\n"
                 << "  --version                       Print the release version\n"
                 << "  --help, -h                      Show this message\n\n";
}

Settings parse_arguments(int argc, char* argv[]) {
    vector<string> input_paths;
    bool input_option_provided = false;
    Settings settings;

    // Get timestamp for fallback folder naming
    string timestamp = settings.get_timestamp();

    const auto require_value = [&](int& index, const std::string& option) -> std::string {
        if (++index >= argc || std::string(argv[index]).rfind("--", 0) == 0)
            throw runtime_error("Error: Missing value for " + option);
        return argv[index];
    };

    // Load the file first so explicit CLI values always override it.
    std::string settings_file;
    for (int i = 1; i < argc; ++i) {
        if (std::string(argv[i]) != "--settings") continue;
        if (!settings_file.empty()) throw runtime_error("Error: --settings may be specified only once");
        settings_file = require_value(i, "--settings");
    }
    if (!settings_file.empty() && !settings.LoadFromFile(settings_file))
        throw runtime_error("Error: could not load settings from " + settings_file);

    bool output_folder_provided = false;
    for (int i = 1; i < argc; ++i) {
        string arg = argv[i];

        if (arg == "--help" || arg == "-h") {
            print_help();
            exit(0);
        }
        if (arg == "--version") {
            cout << "Zuck " ZUCK_VERSION "\n";
            exit(0);
        }
       
        if (arg == "--input-files" || arg == "--input_files" || arg == "-i") {
            input_option_provided = true;
            const auto previous_count = input_paths.size();
            while (++i < argc && argv[i][0] != '-') input_paths.push_back(argv[i]);
            --i;
            if (input_paths.size() == previous_count)
                throw runtime_error("Error: --input-files requires one or two files");
        } else if (arg == "--graph") {
            settings.graph_input = require_value(i, arg);
        } else if (arg == "--cycles") {
            settings.cycles_input = require_value(i, arg);
        } else if (arg == "--retain-graph") {
            settings.retain_graph = true;
        } else if (arg == "--retain-cycles") {
            settings.retain_cycles = true;
        } else if (arg == "--cleanup-cycles") {
            settings.retain_cycles = false;
        } else if (arg == "--experimental-alternative-views") {
            settings.experimental_alternative_views = true;
        } else if (arg == "--boundary-mode") {
            settings.set_boundary_mode(require_value(i, arg));
        } else if (arg == "--boundary-radius") {
            settings.set_boundary_radius(require_value(i, arg));
        } else if (arg == "--boundary-policy") {
            settings.set_boundary_policy(require_value(i, arg));
        } else if (arg == "--call-profile") {
            settings.set_call_profile(require_value(i, arg));
        } else if (arg == "--short-spacer-filter") {
            settings.set_short_spacer_filter(require_value(i, arg));
        } else if (arg == "--terminal-variant-rescue") {
            settings.set_terminal_variant_rescue(require_value(i, arg));
        } else if (arg == "--read-linked-groups") {
            settings.set_read_linked_groups(require_value(i, arg));
        } else if (arg == "--grouping-reads") {
            settings.grouping_reads.clear();
            while (++i < argc && argv[i][0] != '-')
                settings.grouping_reads.push_back(argv[i]);
            --i;
            if (settings.grouping_reads.empty() || settings.grouping_reads.size() > 2)
                throw runtime_error("Error: --grouping-reads requires one or two files");
        } else if (arg == "--array-database") {
            settings.array_database = require_value(i, arg);
        } else if (arg == "--database-read-validation" || arg == "--database-novel-spacers" || arg == "--database-groups" || arg == "--database-include-possible" || arg == "--confidence-scores" || arg == "--candidate-evidence" || arg == "--family-repeat-context" || arg == "--repeat-mixture-refinement" || arg == "--read-recruitment" || arg == "--local-boundary-graph" || arg == "--motif-cliff-boundary" || arg == "--singleton-precision-prune" || arg == "--circular-unit-dedup" || arg == "--directed-primary-groups" || arg == "--cross-variant-scaffolds") {
            const bool value = Settings::database_boolean(require_value(i, arg), arg);
            if (arg == "--database-read-validation") settings.database_read_validation = value;
            else if (arg == "--database-novel-spacers") settings.database_novel_spacers = value;
            else if (arg == "--database-groups") settings.database_groups = value;
            else if (arg == "--confidence-scores") settings.confidence_scores = value;
            else if (arg == "--candidate-evidence") settings.candidate_evidence = value;
            else if (arg == "--family-repeat-context") settings.family_repeat_context = value;
            else if (arg == "--repeat-mixture-refinement") settings.repeat_mixture_refinement = value;
            else if (arg == "--read-recruitment") settings.read_recruitment = value;
            else if (arg == "--local-boundary-graph") settings.local_boundary_graph = value;
            else if (arg == "--motif-cliff-boundary") settings.motif_cliff_boundary = value;
            else if (arg == "--singleton-precision-prune") settings.singleton_precision_prune = value;
            else if (arg == "--circular-unit-dedup") settings.circular_unit_dedup = value;
            else if (arg == "--directed-primary-groups") settings.directed_primary_groups = value;
            else if (arg == "--cross-variant-scaffolds") settings.cross_variant_scaffolds = value;
            else settings.database_include_possible = value;
        } else if (arg == "--benchmark") {
            throw runtime_error("Error: --benchmark is not supported by the production pipeline; use the external scoring workflow");
        } else if (arg == "--ram") {
            settings.ram = Zuck::options::memory_gib(require_value(i, arg));
        } else if (arg == "--threads") {
            settings.threads = Zuck::options::unsigned_integer(
                require_value(i, arg), arg, std::numeric_limits<std::size_t>::max());
        } else if (arg == "--output-folder" || arg == "--output_folder") {
            settings.output_folder = require_value(i, arg);
            output_folder_provided = true;
        } else if (arg == "--cycle-max-length") {
            settings.cycle_finder_settings.cycle_max_length =
                Settings::nonnegative_integer(require_value(i, arg), arg);
        } else if (arg == "--cycle-min-length") {
            settings.cycle_finder_settings.cycle_min_length =
                Settings::nonnegative_integer(require_value(i, arg), arg);
        } else if (arg == "--adaptive-coverage") {
            settings.adaptive_coverage = Zuck::options::boolean(require_value(i, arg), arg);
        } else if (arg == "--read-first") {
            settings.read_first = Zuck::options::boolean(require_value(i, arg), arg);
        } else if (arg == "--read-first-family-recut") {
            settings.read_first_family_recut = Zuck::options::boolean(require_value(i, arg), arg);
        } else if (arg == "--threshold-multiplicity") {
            settings.cycle_finder_settings.threshold_multiplicity =
                Zuck::options::unsigned_integer(require_value(i, arg), arg);
            settings.cycle_finder_settings.threshold_explicit = true;
        } else if (arg == "--low-abundance") {
            settings.cycle_finder_settings.low_abundance =
                Zuck::options::boolean(require_value(i, arg), arg);
            settings.cycle_finder_settings.low_abundance_explicit = true;
        } else if (arg == "--traversal-partitions") {
            settings.cycle_finder_settings.traversal_partitions =
                Settings::nonnegative_integer(require_value(i, arg), arg);
        } else if (arg == "--spacer-min-length" || arg == "--spacer-max-length" ||
                   arg == "--repeat-min-length" || arg == "--repeat-max-length") {
            settings.set_length(arg.substr(2), require_value(i, arg));
        } else if (arg == "--settings") {
            if (++i >= argc)
                throw runtime_error("Error: Missing value for --settings");
            // The file was loaded in the pre-scan; skip its value here.
        } else {
            throw runtime_error("Error: Unknown option: " + arg);
        }
    }
    if (!input_option_provided) input_paths = settings.effective_input_reads();
    settings.set_input_reads(input_paths);

    // Validate resolved values equally whether they came from CLI or a file.
    if (!settings.cycles_input.empty() && !fs::is_regular_file(settings.cycles_input))
        throw runtime_error("Error: Cannot read cycles file: " + settings.cycles_input);
    if (!settings.graph_input.empty()) {
        const auto prefix = fs::is_directory(settings.graph_input)
            ? fs::path(settings.graph_input) / "graph" : fs::path(settings.graph_input);
        if (!fs::is_regular_file(prefix.string() + ".sdbg_info"))
            throw runtime_error("Error: No SDBG graph found at: " + settings.graph_input);
    }

    const int input_mode_count = (!input_paths.empty() ? 1 : 0)
        + (!settings.graph_input.empty() ? 1 : 0)
        + (!settings.cycles_input.empty() ? 1 : 0);
    if (input_mode_count != 1) {
        throw runtime_error(
            "Error: Choose exactly one input mode: --input-files, --graph, or --cycles");
    }
    if (settings.read_first && (settings.call_profile == "off" || settings.effective_grouping_reads().empty()))
        throw runtime_error("Error: --read-first requires reads and an enabled call profile");
    if (settings.read_first_family_recut && !settings.read_first)
        throw runtime_error("Error: --read-first-family-recut requires --read-first true");
    if (settings.candidate_evidence && !settings.two_sided_boundaries)
        throw runtime_error("Error: --candidate-evidence requires two-sided or branch-aware boundaries");
    if (settings.repeat_mixture_refinement && !settings.branch_aware_boundaries)
        throw runtime_error("Error: --repeat-mixture-refinement requires --boundary-mode branch-aware");
    if (settings.read_recruitment && (!settings.two_sided_boundaries || settings.effective_grouping_reads().empty()))
        throw runtime_error("Error: --read-recruitment requires two-sided or branch-aware boundaries and --input-files or --grouping-reads");
    if (settings.local_boundary_graph && (!settings.two_sided_boundaries || settings.effective_grouping_reads().empty()))
        throw runtime_error("Error: --local-boundary-graph requires two-sided or branch-aware boundaries and --input-files or --grouping-reads");
    if (settings.motif_cliff_boundary && !settings.two_sided_boundaries)
        throw runtime_error("Error: --motif-cliff-boundary requires two-sided or branch-aware boundaries");
    if (settings.singleton_precision_prune && settings.call_profile == "off")
        throw runtime_error("Error: --singleton-precision-prune requires an enabled call profile");
    if (settings.circular_unit_dedup && settings.call_profile == "off")
        throw runtime_error("Error: --circular-unit-dedup requires an enabled call profile");
    if (settings.directed_primary_groups && (settings.call_profile == "off" ||
        !settings.two_sided_boundaries || settings.effective_grouping_reads().empty()))
        throw runtime_error("Error: --directed-primary-groups requires an enabled call profile, two-sided boundaries, and --input-files or --grouping-reads");
    if (settings.cross_variant_scaffolds &&
        (!settings.directed_primary_groups || !settings.read_linked_groups))
        throw runtime_error("Error: --cross-variant-scaffolds requires --directed-primary-groups true and --read-linked-groups true");
    if (settings.read_recruitment || settings.local_boundary_graph || settings.directed_primary_groups || settings.cross_variant_scaffolds || settings.confidence_scores || settings.call_profile != "off" || settings.candidate_evidence || (settings.read_linked_groups && settings.two_sided_boundaries) ||
        (!settings.array_database.empty() && settings.database_read_validation)) {
        const auto grouping_inputs = settings.effective_grouping_reads();
        if (grouping_inputs.size() > 2)
            throw runtime_error("Error: read-linked grouping accepts at most two input files");
        for (const auto& path : grouping_inputs)
            if (!fs::is_regular_file(path))
                throw runtime_error("Error: Grouping reads file not found: " + path);
    }
    if (!settings.array_database.empty()) {
        if (!settings.two_sided_boundaries)
            throw runtime_error("Error: --array-database requires two-sided or branch-aware boundaries");
        if (!fs::is_regular_file(settings.array_database))
            throw runtime_error("Error: Array database file not found: " + settings.array_database);
        if (settings.database_read_validation && settings.effective_grouping_reads().empty())
            throw runtime_error("Error: database read validation requires --input-files or --grouping-reads; use --database-read-validation false only for a cycle-only diagnostic");
        if (settings.database_novel_spacers && !settings.database_read_validation)
            throw runtime_error("Error: database novel spacers require read validation");
    }
    // Set default output folder if not provided (only if not set by settings file)
    if (!output_folder_provided && settings.output_folder.empty()) {
        settings.output_folder = "Zuck_run_" + timestamp;
    }

    // Set subfolders and output file (but allow them to be set in settings file)
    if (settings.graph_folder.empty())
        settings.graph_folder = (fs::path(settings.output_folder) / "graph").string();
    if (settings.cycles_folder.empty())
        settings.cycles_folder = (fs::path(settings.output_folder) / "cycles").string();
    if (settings.output_file.empty())
        settings.output_file = (fs::path(settings.output_folder) / "CRISPR_Arrays.txt").string();

    // Protect all source files, including the settings file, from result
    // replacement. Resolve symbolic and hard-link aliases before any writes.
    auto protected_files = settings.effective_input_reads();
    protected_files.insert(protected_files.end(), settings.grouping_reads.begin(), settings.grouping_reads.end());
    if (!settings.array_database.empty()) protected_files.push_back(settings.array_database);
    if (!settings.cycles_input.empty()) protected_files.push_back(settings.cycles_input);
    if (!settings_file.empty()) protected_files.push_back(settings_file);
    if (!settings.graph_input.empty()) {
        const auto prefix = fs::is_directory(settings.graph_input)
            ? fs::path(settings.graph_input) / "graph" : fs::path(settings.graph_input);
        const auto directory = prefix.has_parent_path() ? prefix.parent_path() : fs::path(".");
        const auto filename = prefix.filename().string();
        for (const auto& entry : fs::directory_iterator(directory)) {
            const auto name = entry.path().filename().string();
            if (name == filename + ".sdbg_info" || name.rfind(filename + ".sdbg.", 0) == 0)
                protected_files.push_back(entry.path().string());
        }
    }
    const auto output_directory = fs::weakly_canonical(settings.output_folder);
    for (const auto& input : protected_files) {
        const auto resolved = fs::weakly_canonical(input);
        if (resolved.parent_path() == output_directory &&
            result_bundle::managed_output(resolved.filename().string()))
            throw runtime_error("Error: Input conflicts with a managed result file: " + input);
        if (fs::is_directory(output_directory)) {
            for (const auto& entry : fs::directory_iterator(output_directory))
                if (result_bundle::managed_output(entry.path().filename().string()) &&
                    entry.is_regular_file() && fs::exists(resolved) && fs::equivalent(entry.path(), resolved))
                    throw runtime_error("Error: Input aliases a managed result file: " + input);
        }
        if (!settings.cycles_input.empty()) continue;
        const auto generated_cycles = fs::path(settings.output_folder) / "cycles.txt";
        if (resolved == fs::weakly_canonical(generated_cycles) ||
            (fs::exists(resolved) && fs::exists(generated_cycles) && fs::equivalent(resolved, generated_cycles)))
            throw runtime_error("Error: Input conflicts with generated cycles.txt: " + input);
        const auto relative = resolved.lexically_relative(fs::weakly_canonical(settings.graph_folder));
        if (settings.graph_input.empty() && !relative.empty() && *relative.begin() != "..")
            throw runtime_error("Error: Input is inside the generated graph directory: " + input);
    }

    if (settings.cycles_input.empty() && settings.graph_input.empty() && fs::is_directory(settings.graph_folder)) {
        for (const auto& entry : fs::directory_iterator(settings.graph_folder)) {
            if (!Zuck::graph_artifact(entry.path().filename().string())) continue;
            if (entry.is_symlink() || !entry.is_regular_file())
                throw runtime_error("Error: Unsafe generated graph destination: " + entry.path().string());
            for (const auto& input : protected_files)
                if (fs::exists(input) && fs::equivalent(entry.path(), input))
                    throw runtime_error("Error: Input aliases a generated graph file: " + input);
        }
    }

    // Apply defaults before printing
    if (settings.threads == 0) {
        size_t hc = thread::hardware_concurrency();
        settings.threads = (hc > 2) ? hc - 2 : 1;
    }
    if (settings.ram == 0.0) {
        const double available_ram = get_available_system_ram();
        if (available_ram < 0.0) {
            cerr << "Warning: could not read system RAM, defaulting to 4G" << endl;
            settings.ram = 4.0;
        } else {
            settings.ram = available_ram * 0.85;
            if (settings.ram < 1.0)
                throw runtime_error(
                    "Error: less than 1 GiB is available to this process; "
                    "free memory or provide a larger container limit");
        }
    }

    const double available_ram = get_available_system_ram();
    if (available_ram >= 0.0 && settings.ram > available_ram)
        throw runtime_error("Error: RAM budget exceeds currently available memory");
    if (!settings.print_settings().empty())
        throw runtime_error("Invalid settings; existing output files were preserved.");

    // Print settings summary
    cout << banner << "  ▸ Output      " << settings.output_folder << "\n";
    if (!settings.cycles_input.empty())
        cout << "  ▸ Cycles      " << settings.cycles_input << "  (post-process only)\n";
    else if (!settings.graph_input.empty())
        cout << "  ▸ Graph       " << settings.graph_input << "  (pre-built)\n";
    else
        cout << "  ▸ Input       " << settings.input_files << "\n";
    cout << fixed << setprecision(1)
         << "  ▸ RAM         " << settings.ram << " GiB"
         << "    Threads  " << settings.threads << "\n"
         << "  ▸ Cycles      len " << settings.cycle_finder_settings.cycle_min_length
         << "–" << settings.cycle_finder_settings.cycle_max_length
         << "   threshold ×" << settings.cycle_finder_settings.threshold_multiplicity
         << "   low-abundance " << (settings.cycle_finder_settings.low_abundance ? "on" : "off") << "\n"
         << "  ▸ Alternative views  "
         << (settings.experimental_alternative_views ? "experimental/on" : "off")
         << "\n  ▸ Boundaries  "
         << (settings.branch_aware_boundaries ? "branch-aware (experimental)" :
             settings.two_sided_boundaries ? "two-sided" : "legacy")
         << "\n  ▸ Short-spacer filter  "
         << (!settings.two_sided_boundaries ? "not applied (legacy)" :
             settings.short_spacer_filter ? "on" : "off")
         << "\n  ▸ Terminal-variant rescue  "
         << (!settings.two_sided_boundaries ? "not applied (legacy)" :
             settings.terminal_variant_rescue ? "on" : "off")
         << "\n  ▸ Read-linked groups  "
         << (!settings.two_sided_boundaries ? "not applied (legacy)" :
             !settings.read_linked_groups ? "off" :
             settings.effective_grouping_reads().empty() ? "no reads supplied" : "on")
         << "\n  ▸ Local boundary graph  "
         << (!settings.two_sided_boundaries ? "not applied (legacy)" :
             !settings.local_boundary_graph ? "off" :
             settings.effective_grouping_reads().empty() ? "no reads supplied" : "experimental/on")
         << "\n  ▸ Motif-cliff boundary  "
         << (!settings.two_sided_boundaries ? "not applied (legacy)" :
             settings.motif_cliff_boundary ? "experimental/on" : "off")
         << "\n  ▸ Singleton precision prune  "
         << (settings.singleton_precision_prune ? "experimental/on" : "off")
         << "\n  ▸ Circular-unit dedup  "
         << (settings.circular_unit_dedup ? "experimental/on" : "off")
         << "\n  ▸ Directed primary groups  "
         << (settings.directed_primary_groups ? "experimental/on" : "off")
         << "\n  ▸ Cross-variant scaffolds  "
         << (settings.cross_variant_scaffolds ? "experimental/on" : "off")
         << "\n\n";

    // Validate input files before creating any directories (skipped when --graph is used)
    if (settings.graph_input.empty() && settings.cycles_input.empty()) {
        if (input_paths.size() < 1 || input_paths.size() > 2) {
            throw runtime_error("Error: You must provide one or two input files.");
        }

        static const set<string> valid_exts = {".fa", ".fasta", ".fq", ".fastq", ".gz"};
        for (const auto& file : input_paths) {
            if (!fs::is_regular_file(file)) {
                throw runtime_error("Error: Input file " + file + " does not exist.");
            }
            fs::path p(file);
            if (valid_exts.find(p.extension().string()) == valid_exts.end()) {
                cerr << "Warning: " << file << " has an unexpected extension, expected FASTA/FASTQ (.fa/.fasta/.fq/.fastq/.gz)" << endl;
            }
        }
    }

    // Create directories (after input validation)
    try {
        fs::create_directories(settings.output_folder);
        if (settings.cycles_input.empty()) {
            fs::create_directories(settings.graph_folder);
            fs::create_directories(settings.cycles_folder);
        }
    } catch (const fs::filesystem_error& e) {
        throw runtime_error("Error: Could not create directories: " + string(e.what()));
    }

    return settings;
}

