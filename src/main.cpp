#include <iostream>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <chrono>
#include <cstring>
#include <cctype>
#include <set>
#include <unordered_map>
#include <limits>
#include <optional>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef __linux__
#include <sys/sysinfo.h>
#elif defined(_WIN32)
#include <windows.h>
#elif defined(__APPLE__)
#include <sys/sysctl.h>
#endif
#ifdef DEBUG
#include "io/io_ops.h"
#endif

#include "core/main_run_and_debug.h"
#include "crispr_array/cycle_finder.h"
#include "crispr_array/filters.h"
#include "crispr_array/post_processing.h"
#include "core/pipeline.h"
#include "sdbg/sdbg.h"
#include "core/sdbg_build.h"
#include "core/settings.h"
#include "protospacer_detection/phage_curator.h"
#include "protospacer_detection/isolate_protospacers.h"

using namespace std;
namespace fs = std::filesystem;
#ifdef DEBUG
#pragma message("DEBUG is defined")
#else
#pragma message("DEBUG is NOT defined")
#endif
void print_usage(const char* program_name) {
    cout << "\nmy_mcaat — Metagenomic CRISPR Array Analysis\n\n";
}

bool check_for_error(Settings& settings) {
    string erroneous_message = settings.print_settings();
    if (erroneous_message.empty()) return false;
    cout << "  Errors found in: " << erroneous_message << "\n";
    return true;
}

namespace {

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

Settings parse_arguments(int argc, char* argv[]) {
    vector<string> input_files_default;
    Settings settings;

    // Get timestamp for fallback folder naming
    string timestamp = settings.get_timestamp();

    // Pre-scan for --settings so file values can be used as defaults and CLI overrides them
    for (int j = 1; j < argc; ++j) {
        if (string(argv[j]) == "--settings" && j + 1 < argc) {
            string settings_file = argv[j + 1];
            if (!settings.LoadFromFile(settings_file))
                throw runtime_error("Error: could not load settings from " + settings_file);
            break;
        }
    }

     bool output_folder_provided = false;
 bool required_files_provided = false;
     bool input_files_from_settings = false;
    for (int i = 1; i < argc; ++i) {
        string arg = argv[i];

        if (arg == "--help" || arg == "-h") {
            cout << "\nmy_mcaat — Metagenomic CRISPR Array Analysis\n\n"
                 << "Usage:\n"
                 << "  my_mcaat --input-files <file1> [file2] [options]\n"
                 << "  my_mcaat --graph <path> [options]\n"
                 << "  my_mcaat --cycles <cycles.txt> [options]\n"
                 << "\nInput (one required):\n"
                 << "  --input-files <file1> [file2]   One or two FASTA/FASTQ files (plain or gzipped)\n"
                 << "  --graph <path>                  Pre-built SDBG graph folder (skips construction)\n"
                 << "  --cycles <path>                 Retained cycles.txt (post-process only)\n"
                 << "\nOptions:\n"
                 << "  --output-folder <path>          Output directory  [default: my_mcaat_run_TIMESTAMP]\n"
                 << "  --ram <amount>                  Memory cap: B/K/M/G suffix  [default: 85% available]\n"
                 << "  --threads <num>                 Thread count  [default: CPU cores - 2]\n"
                 << "  --cycle-max-length <int>        Maximum cycle length  [default: 90]\n"
                 << "  --cycle-min-length <int>        Minimum cycle length  [default: 27]\n"
                 << "  --threshold-multiplicity <int>  Edge multiplicity cutoff  [default: 20]\n"
                 << "  --low-abundance <true|false>    Low-abundance mode  [default: true]\n"
                 << "  --boundary-mode <legacy|two-sided>\n"
                 << "                                   Repeat boundary inference [default: legacy]\n"
                 << "                                   two-sided is experimental\n"
                 << "  --experimental-alternative-views\n"
                 << "                                   Emit bounded alternative-repeat sidecar [off]\n"
                 << "  --benchmark <file>              Expected sequences for evaluation\n"
                 << "  --settings <path>               Key=value settings file (CLI overrides)\n"
                 << "  --retain-graph                  Keep the large graph after a full run\n"
                 << "  --cleanup-cycles                Remove cycles.txt after a full run\n"
                 << "  --help, -h                      Show this message\n\n";
            exit(0);
        }
       
        if (arg == "--input-files" || arg == "--input_files" || arg == "-i") {
            while (++i < argc && argv[i][0] != '-') {
                input_files_default.push_back(argv[i]);
                
            }
            --i;
            required_files_provided = true;
        } else if (arg == "--graph") {
            if (++i < argc) {
                settings.graph_input = argv[i];
                // Accept a directory (auto-appends /graph in step_load_graph)
                // or a file prefix — validate by checking the .sdbg_info file exists
                std::string prefix = fs::is_directory(settings.graph_input)
                    ? settings.graph_input + "/graph"
                    : settings.graph_input;
                if (!fs::exists(prefix + ".sdbg_info")) {
                    throw runtime_error("Error: No SDBG graph found at: " + settings.graph_input +
                        "\n  Expected file: " + prefix + ".sdbg_info");
                }
                required_files_provided = true;
            } else {
                throw runtime_error("Error: Missing value for --graph");
            }
        } else if (arg == "--cycles") {
            if (++i < argc) {
                settings.cycles_input = argv[i];
                if (!fs::is_regular_file(settings.cycles_input)) {
                    throw runtime_error("Error: Cycles file not found: " + settings.cycles_input);
                }
                required_files_provided = true;
            } else {
                throw runtime_error("Error: Missing value for --cycles");
            }
        } else if (arg == "--retain-graph") {
            settings.retain_graph = true;
        } else if (arg == "--retain-cycles") {
            settings.retain_cycles = true;
        } else if (arg == "--cleanup-cycles") {
            settings.retain_cycles = false;
        } else if (arg == "--experimental-alternative-views") {
            settings.experimental_alternative_views = true;
        } else if (arg == "--boundary-mode") {
            if (++i >= argc)
                throw runtime_error("Error: Missing value for --boundary-mode");
            settings.set_boundary_mode(argv[i]);
        } else if (arg == "--benchmark") {
            if (++i < argc) {
                settings.benchmark_file = argv[i];
            } else {
                throw runtime_error("Error: Missing value for --benchmark");
            }
        } else if (arg == "--ram") {
            if (++i < argc) {
                string ram_input = argv[i];
                try {
                    double value = 0.0;
                    char unit = 'G';
                    size_t pos = ram_input.find_first_not_of("0123456789.");
                    if (pos != string::npos) {
                        value = stod(ram_input.substr(0, pos));
                        unit = toupper(ram_input[pos]);
                    } else {
                        value = stod(ram_input);
                    }

                    switch (unit) {
                        case 'B': settings.ram = value / (1024.0 * 1024.0 * 1024.0); break;
                        case 'K': settings.ram = value / (1024.0 * 1024.0); break;
                        case 'M': settings.ram = value / 1024.0; break;
                        case 'G': settings.ram = value; break;
                        default:
                            throw runtime_error("Error: Invalid RAM unit. Use B, K, M, or G.");
                    }

                    double available_ram = get_available_system_ram();
                    if (settings.ram < 1.0) {
                        throw runtime_error("Error: RAM value " + to_string(settings.ram) +
                                            " GB is too low (must be at least 1 GB)");
                    }
                    if (available_ram >= 0.0 && settings.ram > available_ram) {
                        throw runtime_error("Error: RAM value " + to_string(settings.ram) +
                                            " GiB exceeds currently available memory of " +
                                            to_string(available_ram) + " GiB");
                    }
                } catch (const invalid_argument&) {
                    throw runtime_error("Error: Invalid RAM value provided: " + ram_input);
                }
            } else {
                throw runtime_error("Error: Missing value for --ram");
            }
        } else if (arg == "--threads") {
            if (++i < argc) {
                try {
                    settings.threads = stoul(argv[i]);
                } catch (...) {
                    throw runtime_error("Error: Invalid value for --threads: " + string(argv[i]));
                }
            } else {
                throw runtime_error("Error: Missing value for --threads");
            }
        } else if (arg == "--output-folder" || arg == "--output_folder") {
            if (++i < argc) {
                settings.output_folder = string(argv[i]);
                output_folder_provided = true;
            } else {
                throw runtime_error("Error: Missing value for --output-folder");
            }
        }
        else if (arg == "--cycle-max-length") {
            if (++i < argc) {
                try {
                    settings.cycle_finder_settings.cycle_max_length = stoi(argv[i]);
                } catch (...) {
                    throw runtime_error("Error: Invalid value for --cycle-max-length: " + string(argv[i]));
                }
            } else {
                throw runtime_error("Error: Missing value for --cycle-max-length");
            }
        } else if (arg == "--cycle-min-length") {
            if (++i < argc) {
                try {
                    settings.cycle_finder_settings.cycle_min_length = stoi(argv[i]);
                } catch (...) {
                    throw runtime_error("Error: Invalid value for --cycle-min-length: " + string(argv[i]));
                }
            } else {
                throw runtime_error("Error: Missing value for --cycle-min-length");
            }
        } else if (arg == "--threshold-multiplicity") {
            if (++i < argc) {
                try {
                    settings.cycle_finder_settings.threshold_multiplicity = stoull(argv[i]);
                } catch (...) {
                    throw runtime_error("Error: Invalid value for --threshold-multiplicity: " + string(argv[i]));
                }
            } else {
                throw runtime_error("Error: Missing value for --threshold-multiplicity");
            }
        } else if (arg == "--low-abundance") {
            if (++i < argc) {
                string value = argv[i];
                std::transform(value.begin(), value.end(), value.begin(), ::tolower);
                if (value == "1" || value == "true" || value == "yes")
                    settings.cycle_finder_settings.low_abundance = true;
                else if (value == "0" || value == "false" || value == "no")
                    settings.cycle_finder_settings.low_abundance = false;
                else
                    throw runtime_error(
                        "Error: Invalid value for --low-abundance: " + value +
                        " (expected true or false)");
            } else {
                throw runtime_error("Error: Missing value for --low-abundance");
            }
        } else if (arg == "--settings") {
            if (++i >= argc)
                throw runtime_error("Error: Missing value for --settings");
            // The file was loaded in the pre-scan; skip its value here.
        } else {
            throw runtime_error("Error: Unknown option: " + arg);
        }
    }
    // If settings file provided input_files, allow it as equivalent to CLI --input-files
    if (input_files_default.empty() && !settings.input_files.empty()) {
        std::istringstream iss(settings.input_files);
        string token;
        while (iss >> token) {
            input_files_default.push_back(token);
        }
        required_files_provided = true;
        input_files_from_settings = true;
    }

    const int input_mode_count = (!input_files_default.empty() ? 1 : 0)
        + (!settings.graph_input.empty() ? 1 : 0)
        + (!settings.cycles_input.empty() ? 1 : 0);
    if (input_mode_count != 1) {
        throw runtime_error(
            "Error: Choose exactly one input mode: --input-files, --graph, or --cycles");
    }
    // Set default output folder if not provided (only if not set by settings file)
    if (!output_folder_provided && settings.output_folder.empty()) {
        settings.output_folder = "my_mcaat_run_" + timestamp;
    }

    // Set subfolders and output file (but allow them to be set in settings file)
    if (settings.graph_folder.empty())
        settings.graph_folder = (fs::path(settings.output_folder) / "graph").string();
    if (settings.cycles_folder.empty())
        settings.cycles_folder = (fs::path(settings.output_folder) / "cycles").string();
    if (settings.output_file.empty())
        settings.output_file = (fs::path(settings.output_folder) / "CRISPR_Arrays.txt").string();

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

    // Print settings summary
    cout << "  ▸ Output      " << settings.output_folder << "\n";
    if (!settings.cycles_input.empty())
        cout << "  ▸ Cycles      " << settings.cycles_input << "  (post-process only)\n";
    else if (!settings.graph_input.empty())
        cout << "  ▸ Graph       " << settings.graph_input << "  (pre-built)\n";
    else
        cout << "  ▸ Input       " << settings.input_files << "\n";
    cout << fixed << setprecision(1)
         << "  ▸ RAM         " << settings.ram << " GB"
         << "    Threads  " << settings.threads << "\n"
         << "  ▸ Cycles      len " << settings.cycle_finder_settings.cycle_min_length
         << "–" << settings.cycle_finder_settings.cycle_max_length
         << "   threshold ×" << settings.cycle_finder_settings.threshold_multiplicity
         << "   low-abundance " << (settings.cycle_finder_settings.low_abundance ? "on" : "off") << "\n"
         << "  ▸ Alternative views  "
         << (settings.experimental_alternative_views ? "experimental/on" : "off")
         << "\n  ▸ Boundaries  "
         << (settings.two_sided_boundaries ? "two-sided (experimental)" : "legacy")
         << "\n\n";

    // Validate input files before creating any directories (skipped when --graph is used)
    if (settings.graph_input.empty() && settings.cycles_input.empty()) {
        if (input_files_default.size() < 1 || input_files_default.size() > 2) {
            throw runtime_error("Error: You must provide one or two input files.");
        }

        static const set<string> valid_exts = {".fa", ".fasta", ".fq", ".fastq", ".gz"};
        for (const auto& file : input_files_default) {
            if (!fs::exists(file)) {
                throw runtime_error("Error: Input file " + file + " does not exist.");
            }
            fs::path p(file);
            if (valid_exts.find(p.extension().string()) == valid_exts.end()) {
                cerr << "Warning: " << file << " has an unexpected extension, expected FASTA/FASTQ (.fa/.fasta/.fq/.fastq/.gz)" << endl;
            }
            if (required_files_provided && !input_files_from_settings) {
                if (!settings.input_files.empty()) settings.input_files += " ";
                settings.input_files += file;
            }
        }
    } else {
        // graph_input already shown in settings summary above
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


string fetchNodeLabel(Settings& settings, uint64_t node) {
    std::string label;            
    uint8_t seq[settings.sdbg->k()];
    settings.sdbg->GetLabel(node, seq);
    for (int i = settings.sdbg->k() - 1; i >= 0; --i) label.append(1, "ACGT"[seq[i] - 1]);
    reverse(label.begin(), label.end());
    return label;
}

uint64_t findNodeFromKmer(Settings& settings, const std::string& kmer) {
    int k = settings.sdbg->k();
    if (static_cast<int>(kmer.size()) != k) {
        std::cerr << "Warning: kmer size " << kmer.size() << " does not match k=" << k << std::endl;
        return UINT64_MAX; // or handle error
    }
    std::vector<uint8_t> seq(k);
    for (int i = 0; i < k; ++i) {
        seq[i] = "ACGT"s.find(kmer[i]) + 1;
    }
    int64_t result = settings.sdbg->IndexBinarySearch(seq.data());
    return (result == -1) ? UINT64_MAX : static_cast<uint64_t>(result);
}

std::map<uint64_t, std::vector<std::vector<uint64_t>>> createRepeatToSpacerNodes(Settings& settings, const std::map<std::string, std::vector<std::string>>& systems_from_analyzer) {
    std::map<uint64_t, std::vector<std::vector<uint64_t>>> repeat_to_spacer_nodes;
    int k = settings.sdbg->k();
    for (const auto& [repeat, spacers] : systems_from_analyzer) {
        if (static_cast<int>(repeat.size()) < k) continue;
        std::string first_kmer = repeat.substr(0, k);
        uint64_t key_node = findNodeFromKmer(settings, first_kmer);
        if (key_node == UINT64_MAX) continue; // invalid
        std::vector<std::vector<uint64_t>> spacer_node_vectors;
        for (const auto& spacer : spacers) {
            std::vector<uint64_t> nodes;
            int L = spacer.size();
            for (int i = 0; i <= L - k; ++i) {
                std::string kmer = spacer.substr(i, k);
                uint64_t node = findNodeFromKmer(settings, kmer);
                if (node != UINT64_MAX) {
                    nodes.push_back(node);
                }
            }
            if (!nodes.empty()) {
                spacer_node_vectors.push_back(std::move(nodes));
            }
        }
        if (!spacer_node_vectors.empty()) {
            repeat_to_spacer_nodes[key_node] = std::move(spacer_node_vectors);
        }
    }
    return repeat_to_spacer_nodes;
}

#ifdef DEBUG
int main(int argc, char** argv) {
    // %% PARSE ARGUMENTS %%
    Settings settings = parse_arguments(argc, argv);
    string name_of_genome = "test";
    if (check_for_error(settings)){
        cerr << "Invalid settings; existing output files were preserved.\n";
        return 1;
    }
    // %% PARSE ARGUMENTS %%

    // %% BUILD GRAPH %%
    SDBGBuild sdbg_build(settings);
    // %% BUILD GRAPH %%
    
   
    // cycle finder max/min length are read from settings.cycle_finder_settings
    SDBG sdbg;
    vector<string> folders = {"/vol/d/development/git/mcaat_master/mcaat/_build/mcaat_run_2025-08-27_09-53-11/graph/graph","/vol/d/development/git/mcaat_master/mcaat/_build/mcaat_run_2025-08-28_13-26-39/graph/graph",
    "/vol/d/development/git/mcaat_master/mcaat/build/mcaat_run_2025-10-24_12-27-47/graph/graph","/vol/d/data/real/graph"};
    string graph_folder_old = settings.graph_folder;///vol/d/development/git/mcaat_master/mcaat/_build/mcaat_run_2025-08-28_13-23-46
    settings.graph_folder=settings.graph_folder+"/graph";
    char * cstr = new char [settings.graph_folder.length()+1];
    std::strcpy (cstr, settings.graph_folder.c_str());
    cout << "Graph folder: " << cstr << endl;
    sdbg.LoadFromFile(cstr);
    cout << "Loaded the graph" << endl;
    settings.sdbg = &sdbg;

    // %% LOAD GRAPH %%
    
    delete[] cstr;
    string ending_seq = "CAGAGATAGAAATTATTTTTATTATACGTTTTTTTGT";
    vector<int64_t> end_nodes;
    //using findNodeFromKmer find all the node ids in the ending_seq
    for (size_t i = 0; i <= ending_seq.size() - settings.sdbg->k(); ++i) {
        string kmer = ending_seq.substr(i, sdbg.k());
        uint64_t node = findNodeFromKmer(settings, kmer);
        if (node != UINT64_MAX) {
            end_nodes.push_back(node);
        }
    }
    //print all the end nodes
    cout << "End nodes: ";
    for (const auto& node : end_nodes) {
        cout << node << " ";
    }
    cout << endl;
    // make some distribution of multiplicities of all nodes in the graph
    std::map<int, int> multiplicity_distribution;
    for (uint64_t node = 0; node < sdbg.size(); ++node) {
        int mult = sdbg.EdgeMultiplicity(node);
        multiplicity_distribution[mult]++;
    }
    cout << "Node Multiplicity Distribution:" << endl;
    //write distribution to file
    ofstream mult_file("node_multiplicities.txt");
    for (const auto& [mult, count] : multiplicity_distribution) {
        mult_file << "Multiplicity " << mult << ": " << count << " nodes" << endl;
    }

    mult_file.close();
    // for (uint64_t node = 0; node < sdbg.size(); ++node) {
    //     int mult = sdbg.EdgeMultiplicity(node);
    //     mult_file << node << "\t" << mult << "\n";
    // }
    // mult_file.close();
    // %% FBCE ALGORITHM %%
    cout << "FBCE FROM DEBUG START:" << endl;
    auto start_time = chrono::high_resolution_clock::now();
    CycleFinder cycle_finder(settings);
    // number_of_spacers_total not used; remove to avoid compiler warning
    auto cycles = cycle_finder.results;
    cout << "Number of nodes in results: " << cycles.size() << endl;
    // %% FBCE ALGORITHM %%
    
    int number_of_spacers = 0;
    
    // %% FILTERS %%
    cout << "FILTERS START:" << endl;
    Filters filters(sdbg, cycles);
    auto  SYSTEMS = filters.ListArrays(number_of_spacers);
    cout<< "Number of spacers: " << number_of_spacers << " before cleaning"<<endl;
    // %% FILTERS %%

    //%% POST PROCESSING %%
    cout << "POST PROCESSING START:" << endl;
    CRISPRAnalyzer analyzer(SYSTEMS, settings.output_file);
    analyzer.run_analysis();
    cout << "Saved in: " << settings.output_file << endl;
    auto systems_from_analyzer = analyzer.getSystems();
    auto repeat_to_spacer_nodes = createRepeatToSpacerNodes(settings, systems_from_analyzer);
    cout << "Created repeat_to_spacer_nodes map with " << repeat_to_spacer_nodes.size() << " entries." << endl;
    //%% POST PROCESSING %%

    //%% PROTOSPACER ISOLATION %%
    IsolateProtospacers isolator(sdbg, repeat_to_spacer_nodes);
    pair<std::map<uint64_t,std::set<uint64_t>>,std::map<uint64_t,std::set<uint64_t>>> protospacer_nodes = isolator.getProtospacerNodes();
    auto grouped_paths_protospacers = isolator.DepthLimitedPathsFromInToOut(protospacer_nodes.first, protospacer_nodes.second, 50,1);
       isolator.WritePathsToFile(grouped_paths_protospacers, "grouped_paths_protospacers.txt");

    //auto grouped_paths_protospacers = isolator.ReadPathsFromFile("grouped_paths_protospacers.txt");
    
    //%% PROTOSPACER ISOLATION %%
    PhageCurator phage_curator(sdbg, grouped_paths_protospacers, cycles);

    // Find and output quality paths using DLS with multiplicity filters, writing live to file
    //phage_curator.FindQualityPathsDLSFromGroupedPaths(3000, 3010, "QualityPaths.fasta");  // min_length 3000, max_length arbitrary large
    vector<int> beam_widths = {50};
    // Find and output quality paths using beam search with multiplicity filters, writing live to file
    for (int beam_width : beam_widths) {
        string filename = "QualityPaths_BeamWidth" + to_string(beam_width) + ".fasta";
        phage_curator.FindQualityPathsBeamSearchFromGroupedPaths(3000, 3010, filename, beam_width);  // min_length 3000, max_length arbitrary large
    }
    //print multiplicities of nodes into a file

    //io_ops::write_nodes_gfa("output.gfa", sdbg);
    //io_ops::write_nodes_gfa("output.gfa", sdbg);
    

    //%% POST PROCESSING %%
    // %% DELETE THE GRAPH FOLDER %%
    //fs::remove_all(graph_folder_old);
    // %% DELETE THE GRAPH FOLDER %%          
    
}

#else
int main(int argc, char** argv) {
    try {
        cout << "\nmy_mcaat — Metagenomic CRISPR Array Analysis\n\n";

        Settings settings = parse_arguments(argc, argv);
        if (check_for_error(settings)) {
            cerr << "Invalid settings; existing output files were preserved.\n";
            return 1;
        }

        auto t_start = chrono::high_resolution_clock::now();
        auto step_elapsed = [](chrono::high_resolution_clock::time_point t0) {
            return chrono::duration<double>(chrono::high_resolution_clock::now() - t0).count();
        };

        if (!settings.cycles_input.empty()) {
            auto t0 = chrono::high_resolution_clock::now();
            cout << "[1/1] Post-processing retained cycles...\n";
            step_post_process(settings);
            cout << "      done in " << fixed << setprecision(1) << step_elapsed(t0) << " s\n\n"
                 << "Output: " << settings.output_folder << "\n\n";
            return 0;
        }

        SDBG sdbg;
        auto t0 = chrono::high_resolution_clock::now();
        if (settings.graph_input.empty()) {
            cout << "[1/4] Building de Bruijn graph...\n";
            step_build_graph(settings);
        } else {
            cout << "[1/4] Using pre-built graph\n";
        }
        step_load_graph(settings, sdbg);
        cout << "      done in " << fixed << setprecision(1) << step_elapsed(t0) << " s\n\n";

        t0 = chrono::high_resolution_clock::now();
        cout << "[2/4] Finding CRISPR cycles...\n";
        step_find_cycles(settings);
        cout << "      done in " << fixed << setprecision(1) << step_elapsed(t0) << " s\n\n";

        t0 = chrono::high_resolution_clock::now();
        cout << "[3/4] Post-processing arrays...\n";
        step_post_process(settings);
        cout << "      done in " << fixed << setprecision(1) << step_elapsed(t0) << " s\n\n";

        cout << "[4/4] Cleaning up...\n";
        step_cleanup(settings);

        double total = step_elapsed(t_start);
        cout << "\nDone in " << fixed << setprecision(1) << total << " s\n"
             << "Output: " << settings.output_folder << "\n\n";

        return 0;
    } catch (const std::exception& error) {
        cerr << error.what() << "\n";
        return 2;
    }
}
#endif
