#ifndef SETTINGS_H
#define SETTINGS_H

#include <iostream>
#include <string>
#include <map>
#include <vector>
#include <thread>
#include <filesystem>
#include <chrono>
#include <iomanip>
#include <sstream>
#include <algorithm>
#include <stdexcept>
#include <cctype>
#include <cstdint>
#include <limits>
#include <cmath>
#include "core/option_values.h"
#include "sdbg/sdbg.h"
#include <fstream>

using namespace std;
namespace fs = std::filesystem;

struct Settings {
    std::string input_files; // Human-readable input summary (legacy callers may set this).
    std::vector<std::string> input_paths; // Preserve exact paths, including spaces.
    std::string graph_input;  // Path to a pre-built SDBG graph folder (skips build step)
    std::string cycles_input; // Retained legacy cycles.txt (post-process only)
    bool retain_graph = false; // Graphs are large; remove them unless explicitly retained
    bool retain_cycles = true; // Keep compact cycle evidence for fast replay
    bool adaptive_coverage = false;
    bool read_first = false;
    bool read_first_family_recut = false;
    bool read_first_family_admission = false;
    bool read_first_shallow_rescue = false;
    struct CoverageSupport {
        bool measured = false;
        bool low_support = false;
        uint64_t real_edges = 0;
        uint64_t repeated_edges = 0;
        unsigned q25 = 0;
        unsigned q50 = 0;
        unsigned q75 = 0;
    } coverage_support;
    // Experimental, bounded views over already accepted primary rows.
    bool experimental_alternative_views = false;
    // Correct circular representation and infer both repeat boundaries.
    // Best measured recovery/concordance on the frozen development benchmark.
    bool two_sided_boundaries = true;
    // Experimental divergent-family separation before two-sided extraction.
    bool branch_aware_boundaries = false;
    // Experimental nested repeat-family mixtures; requires branch-aware mode.
    bool repeat_mixture_refinement = false;
    // Reject groups explained by one spacer below the 23-mer collapse length.
    bool short_spacer_filter = false;
    // Preserve bounded, equally supported terminal spacer variants.
    bool terminal_variant_rescue = true;
    // Use raw read junctions to report linked subgroups without replacing parents.
    bool read_linked_groups = true;
    std::vector<std::string> grouping_reads;
    // Optional, confidence-bearing CRISPRidentify array evidence export.
    std::string array_database;
    bool database_read_validation = true;
    bool database_novel_spacers = false;
    bool database_groups = true;
    bool database_include_possible = false;
    // Annotate observed elements; confidence never silently deletes evidence.
    bool confidence_scores = true;
    // Additional selected-call view; inclusive primary rows remain unchanged.
    std::string call_profile = "standard";
    // Experimental support for adjacent observed repeat alleles within a family.
    bool family_repeat_context = false;
    // Graph-seeded read discovery and independently witnessed junction support.
    bool read_recruitment = false;
    // Experimental read-local repeat/spacer recutting before grouping.
    bool local_boundary_graph = false;
    // Sequence-only terminal motif-cliff repair over exact oriented haplotypes.
    bool motif_cliff_boundary = false;
    // Demote threshold-edge types found only in selected singleton families.
    bool singleton_precision_prune = false;
    // Emit one evidence-ranked spacer partition per exact circular DNA unit.
    bool circular_unit_dedup = false;
    // Complete row-disjoint partition from directed exact-haplotype read paths.
    bool directed_primary_groups = false;
    // Merge directed components only through native nonbranched read-link trees.
    bool cross_variant_scaffolds = false;
    // Optional evidence for every exact de novo pair, before allele collapse.
    bool candidate_evidence = false;
    int boundary_radius = 3;
    std::string boundary_policy = "catalogue";
    double ram = 0.0; // Maximum RAM budget in GiB
    size_t threads = 0; // Number of threads to use
    std::string output_folder; // Output directory path. Empty until chosen by CLI or settings file.
    std::string graph_folder; // Folder for graph data
    std::string cycles_folder; // Folder for cycle data
    std::string output_file; // Path to the main output file
    std::string benchmark_file;
    // Sdbg

    struct CycleFinderSettings {
        uint64_t threshold_multiplicity = 20; // Minimum multiplicity threshold
        bool low_abundance = false; // Preserve historical effective singleton removal
        int traversal_partitions = 8; // Logical search histories, independent of workers
        int cycle_max_length = 90; // Maximum repeat+spacer cycle length
        int cycle_min_length = 27; // Minimum length of a cycle
        bool threshold_explicit = false;
        bool low_abundance_explicit = false;
    } cycle_finder_settings;
    struct DNASequenceSettings {
        // Optional output constraints. Zero leaves that bound unconstrained;
        // the existing discovery/boundary limits still apply.
        int spacer_min_length = 0;
        int spacer_max_length = 0;
        int repeat_min_length = 0;
        int repeat_max_length = 0;

        bool active() const {
            return spacer_min_length || spacer_max_length || repeat_min_length || repeat_max_length;
        }
        bool repeat_allowed(std::size_t n) const {
            return (!repeat_min_length || n >= static_cast<std::size_t>(repeat_min_length)) &&
                   (!repeat_max_length || n <= static_cast<std::size_t>(repeat_max_length));
        }
        bool pair_allowed(const std::string& repeat, const std::string& spacer) const {
            return repeat_allowed(repeat.size()) &&
                   (!spacer_min_length || spacer.size() >= static_cast<std::size_t>(spacer_min_length)) &&
                   (!spacer_max_length || spacer.size() <= static_cast<std::size_t>(spacer_max_length));
        }
    } dna_sequence_settings;

    SDBG* sdbg = nullptr; // Pointer to the SDBG graph

    void set_boundary_mode(const std::string& mode) {
        if (mode != "legacy" && mode != "two-sided" && mode != "branch-aware")
            throw std::runtime_error(
            "Invalid boundary mode: " + mode +
            " (expected legacy, two-sided, or branch-aware)");
        two_sided_boundaries = mode != "legacy";
        branch_aware_boundaries = mode == "branch-aware";
    }

    void set_boundary_radius(const std::string& value) {
        boundary_radius = nonnegative_integer(value, "boundary-radius");
        if (boundary_radius > 3) throw std::runtime_error("boundary-radius must be in [0,3]");
    }
    void set_call_profile(const std::string& value) {
        if (value != "standard" && value != "context" && value != "off")
            throw std::runtime_error("call-profile must be standard, context, or off");
        call_profile = value;
    }
    void set_boundary_policy(const std::string& value) {
        if (value != "off" && value != "catalogue" && value != "catalogue-all" && value != "catalogue-best" &&
            value != "anchored" && value != "anchored-best")
            throw std::runtime_error("boundary-policy must be off, catalogue, catalogue-all, catalogue-best, anchored, or anchored-best");
        boundary_policy = value;
    }

    void set_short_spacer_filter(const std::string& value) {
        short_spacer_filter = Zuck::options::boolean(value, "short-spacer-filter");
    }

    void set_terminal_variant_rescue(const std::string& value) {
        terminal_variant_rescue = Zuck::options::boolean(value, "terminal-variant-rescue");
    }

    void set_read_linked_groups(const std::string& value) {
        read_linked_groups = Zuck::options::boolean(value, "read-linked-groups");
    }

    static bool database_boolean(const std::string& value, const std::string& name) {
        return Zuck::options::boolean(value, name);
    }

    static int nonnegative_integer(const std::string& value, const std::string& name) {
        return Zuck::options::nonnegative_integer(value, name);
    }

    void set_input_reads(std::vector<std::string> paths) {
        input_paths = std::move(paths);
        std::ostringstream summary;
        for (const auto& path : input_paths) {
            if (summary.tellp() > 0) summary << ' ';
            summary << std::quoted(path);
        }
        input_files = summary.str();
    }

    std::vector<std::string> effective_input_reads() const {
        return input_paths.empty() ? Zuck::options::paths(input_files) : input_paths;
    }

    void set_length(const std::string& name, const std::string& value) {
        const int n = nonnegative_integer(value, name);
        if (name == "spacer-min-length") dna_sequence_settings.spacer_min_length = n;
        else if (name == "spacer-max-length") dna_sequence_settings.spacer_max_length = n;
        else if (name == "repeat-min-length") dna_sequence_settings.repeat_min_length = n;
        else if (name == "repeat-max-length") dna_sequence_settings.repeat_max_length = n;
        else throw std::runtime_error("Unknown length setting: " + name);
    }

    std::vector<std::string> effective_grouping_reads() const {
        return grouping_reads.empty() ? effective_input_reads() : grouping_reads;
    }

    // Generate timestamp in YYYY-MM-DD_HH-MM-SS format
    string get_timestamp() {
        auto now = chrono::system_clock::now();
        auto in_time_t = chrono::system_clock::to_time_t(now);
        stringstream ss;
        ss << put_time(localtime(&in_time_t), "%Y-%m-%d_%H-%M-%S");
        return ss.str();
    }

    // Get maximum available threads
    size_t max_threads() const {
        return std::max(1u, thread::hardware_concurrency());
    }

    // Validate settings and return status map
    map<string, pair<bool, string>> validate_settings() const {
        map<string, pair<bool, string>> settings_message_map;
        
        // Check input files — skipped when a pre-built graph is provided
        const auto reads = effective_input_reads();
        const bool input_valid = !cycles_input.empty() || !graph_input.empty() ||
            (!reads.empty() && reads.size() <= 2 &&
             std::all_of(reads.begin(), reads.end(), [](const auto& path) {
                 return fs::is_regular_file(path);
             }));
        settings_message_map["Input Files"] = make_pair(input_valid,
            input_valid ? (!cycles_input.empty() ? "(cycles replay mode)" :
                           (graph_input.empty() ? input_files : "(graph mode)"))
                        : "File(s) not found: " + input_files);

        // RAM is a finite GiB budget; zero is resolved to an automatic budget before validation.
        string ram_str = to_string(ram).substr(0, to_string(ram).find(".") + 3);
        bool ram_valid = std::isfinite(ram) && ram >= 1.0;
        settings_message_map["RAM"] = make_pair(ram_valid, 
            ram_valid ? ram_str + " GB" : 
            "Value " + ram_str + " GB is invalid (must be at least 1 GB)");

        // Check threads
        size_t max_t = max_threads();
        bool threads_valid = threads > 0 && threads <= max_t;
        settings_message_map["Threads"] = make_pair(threads_valid, 
            threads_valid ? to_string(threads) + " thread(s)" : 
            "Value " + to_string(threads) + " is invalid (must be between 1 and " + to_string(max_t) + ")");

        const auto& dna = dna_sequence_settings;
        const bool lengths_valid =
            (!dna.spacer_max_length || dna.spacer_min_length <= dna.spacer_max_length) &&
            (!dna.repeat_max_length || dna.repeat_min_length <= dna.repeat_max_length);
        settings_message_map["Sequence lengths"] = {lengths_valid,
            lengths_valid ? "Optional output constraints" : "Minimum length exceeds maximum length"};
        const auto& search = cycle_finder_settings;
        settings_message_map["Cycle search"] = {
            search.cycle_max_length > search.cycle_min_length && search.cycle_min_length >= 0 &&
            search.traversal_partitions > 0,
            "Require maximum > minimum >= 0 and traversal-partitions > 0"};

        // Check folder paths
        bool output_valid = !output_folder.empty();
        settings_message_map["Output Folder"] = make_pair(output_valid, 
            output_valid ? output_folder : "Invalid output folder");

        return settings_message_map;
    }

    // Print settings summary and return erroneous properties
    string print_settings() {
        string erroneous_properties;
        auto settings_map = validate_settings();
        for (const auto& [key, value] : settings_map) {
            if (!value.first) {
                erroneous_properties += key + " ";
                cout << "  [✗] " << key << ": " << value.second << "\n";
            }
        }
        return erroneous_properties;
    }

    // Read settings values from a simple key=value file. Lines starting with # or // are ignored.
    // Keys match CLI flag names (dashes); underscore variants are also accepted.
    // Examples:
    // input-files=/path/a.fa /path/b.fa
    // graph=/path/to/graph
    // cycles=/path/to/cycles.txt
    // ram=4G
    // threads=4
    // cycle-max-length=77
    // cycle-min-length=27
    // threshold-multiplicity=20
    // low-abundance=true
    bool LoadFromFile(const std::string& path) {
        std::ifstream file(path);
        if (!file.is_open()) {
            std::cerr << "Could not open settings file: " << path << std::endl;
            return false;
        }
        std::string line;
        auto trim = [](std::string s) {
            // trim in place
            const char* ws = " \t\n\r\f\v";
            s.erase(0, s.find_first_not_of(ws));
            s.erase(s.find_last_not_of(ws) + 1);
            return s;
        };
        while (std::getline(file, line)) {
            char quote = 0;
            for (std::size_t i = 0; i < line.size(); ++i) {
                const char c = line[i];
                if (quote) { if (c == quote) quote = 0; }
                else if (c == '\'' || c == '"') quote = c;
                else if ((i == 0 || std::isspace(static_cast<unsigned char>(line[i - 1]))) &&
                         (c == '#' || line.compare(i, 2, "//") == 0)) {
                    line.resize(i);
                    break;
                }
            }
            string s = trim(line);
            if (s.empty()) continue;
            size_t eq = s.find('=');
            if (eq == std::string::npos) continue;
            string key = trim(s.substr(0, eq));
            string val = trim(s.substr(eq + 1));

            // Normalize key: replace underscores with dashes so both formats are accepted
            std::replace(key.begin(), key.end(), '_', '-');

            // Interpret known keys
            if (key == "input-files") {
                set_input_reads(Zuck::options::paths(val));
            } else if (key == "graph") {
                graph_input = val;
            } else if (key == "cycles") {
                cycles_input = val;
            } else if (key == "ram") {
                ram = Zuck::options::memory_gib(val);
            } else if (key == "threads") {
                threads = Zuck::options::unsigned_integer(val, key, std::numeric_limits<std::size_t>::max());
            } else if (key == "output-folder") { this->output_folder = val; }
            else if (key == "graph-folder") { this->graph_folder = val; }
            else if (key == "cycles-folder") { this->cycles_folder = val; }
            else if (key == "output-file") { this->output_file = val; }
            else if (key == "cycle-max-length") { this->cycle_finder_settings.cycle_max_length = nonnegative_integer(val, key); }
            else if (key == "cycle-min-length") { this->cycle_finder_settings.cycle_min_length = nonnegative_integer(val, key); }
            else if (key == "adaptive-coverage") { adaptive_coverage = database_boolean(val, key); }
            else if (key == "read-first") { read_first = database_boolean(val, key); }
            else if (key == "read-first-family-recut") { read_first_family_recut = database_boolean(val, key); }
            else if (key == "read-first-family-admission") { read_first_family_admission = database_boolean(val, key); }
            else if (key == "read-first-shallow-rescue") { read_first_shallow_rescue = database_boolean(val, key); }
            else if (key == "threshold-multiplicity") {
                cycle_finder_settings.threshold_multiplicity = Zuck::options::unsigned_integer(val, key);
                cycle_finder_settings.threshold_explicit = true;
            }
            else if (key == "low-abundance") {
                cycle_finder_settings.low_abundance = database_boolean(val, key);
                cycle_finder_settings.low_abundance_explicit = true;
            }
            else if (key == "traversal-partitions") {
                cycle_finder_settings.traversal_partitions = nonnegative_integer(val, key);
            }
            else if (key == "retain-graph") { retain_graph = database_boolean(val, key); }
            else if (key == "retain-cycles") { retain_cycles = database_boolean(val, key); }
            else if (key == "experimental-alternative-views") {
                experimental_alternative_views = database_boolean(val, key);
            }
            else if (key == "boundary-mode") { set_boundary_mode(val); }
            else if (key == "repeat-mixture-refinement") { repeat_mixture_refinement = database_boolean(val, key); }
            else if (key == "boundary-radius") { set_boundary_radius(val); }
            else if (key == "boundary-policy") { set_boundary_policy(val); }
            else if (key == "short-spacer-filter") { set_short_spacer_filter(val); }
            else if (key == "terminal-variant-rescue") { set_terminal_variant_rescue(val); }
            else if (key == "read-linked-groups") { set_read_linked_groups(val); }
            else if (key == "array-database") { array_database = val; }
            else if (key == "database-read-validation") { database_read_validation = database_boolean(val, key); }
            else if (key == "database-novel-spacers") { database_novel_spacers = database_boolean(val, key); }
            else if (key == "database-groups") { database_groups = database_boolean(val, key); }
            else if (key == "database-include-possible") { database_include_possible = database_boolean(val, key); }
            else if (key == "confidence-scores") { confidence_scores = database_boolean(val, key); }
            else if (key == "call-profile") { set_call_profile(val); }
            else if (key == "family-repeat-context") { family_repeat_context = database_boolean(val, key); }
            else if (key == "read-recruitment") { read_recruitment = database_boolean(val, key); }
            else if (key == "local-boundary-graph") { local_boundary_graph = database_boolean(val, key); }
            else if (key == "motif-cliff-boundary") { motif_cliff_boundary = database_boolean(val, key); }
            else if (key == "singleton-precision-prune") { singleton_precision_prune = database_boolean(val, key); }
            else if (key == "circular-unit-dedup") { circular_unit_dedup = database_boolean(val, key); }
            else if (key == "directed-primary-groups") { directed_primary_groups = database_boolean(val, key); }
            else if (key == "cross-variant-scaffolds") { cross_variant_scaffolds = database_boolean(val, key); }
            else if (key == "candidate-evidence") { candidate_evidence = database_boolean(val, key); }
            else if (key == "grouping-reads") {
                grouping_reads = Zuck::options::paths(val);
                if (grouping_reads.empty() || grouping_reads.size() > 2)
                    throw std::runtime_error("grouping-reads requires one or two files");
            }
            else if (key == "spacer-min-length" || key == "spacer-max-length" ||
                     key == "repeat-min-length" || key == "repeat-max-length") { set_length(key, val); }
            // unknown keys are ignored for forward-compatibility
        }
        if (file.bad()) throw std::runtime_error("Error reading settings file: " + path);
        return true;
    }
};

#endif // SETTINGS_H
