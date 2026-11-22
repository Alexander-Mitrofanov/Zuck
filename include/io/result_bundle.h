#ifndef ZUCK_RESULT_BUNDLE_H
#define ZUCK_RESULT_BUNDLE_H

#include <algorithm>
#include <array>
#include <chrono>
#include <cctype>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>
#include "xxhash/xxhash.h"

namespace result_bundle {
namespace fs = std::filesystem;

inline std::string json_string(const std::string& value) {
    std::ostringstream out;
    out << '"';
    for (unsigned char c : value) {
        if (c == '"' || c == '\\') out << '\\' << c;
        else if (c < 0x20)
            out << "\\u" << std::hex << std::setw(4) << std::setfill('0')
                << static_cast<unsigned int>(c);
        else out << c;
    }
    out << '"';
    return out.str();
}

// XXH64(seed=0) is a fast file identity check, not a cryptographic digest.
inline std::string file_checksum(const fs::path& path) {
    std::ifstream input(path, std::ios::binary);
    if (!input || !fs::is_regular_file(path))
        throw std::runtime_error("Cannot checksum result/input file: " + path.string());
    std::unique_ptr<XXH64_state_t, decltype(&XXH64_freeState)>
        state(XXH64_createState(), &XXH64_freeState);
    if (!state) throw std::bad_alloc();
    XXH64_reset(state.get(), 0);
    std::array<char, 65536> buffer;
    while (input.read(buffer.data(), buffer.size()) || input.gcount())
        XXH64_update(state.get(), buffer.data(), static_cast<std::size_t>(input.gcount()));
    if (!input.eof())
        throw std::runtime_error("Error checksumming result/input file: " + path.string());
    std::ostringstream out;
    out << std::hex << std::setw(16) << std::setfill('0') << XXH64_digest(state.get());
    return out.str();
}

inline void close_checked(std::ofstream& out, const fs::path& path) {
    out.close();
    if (!out) throw std::runtime_error("Error writing result file: " + path.string());
}

inline bool numbered_output(const std::string& name) {
    const std::string prefix = "CRISPR_Arrays_";
    if (name.size() <= prefix.size() + 4 || name.compare(0, prefix.size(), prefix) ||
        name.compare(name.size() - 4, 4, ".txt")) return false;
    return std::all_of(name.begin() + prefix.size(), name.end() - 4,
                       [](unsigned char c) { return std::isdigit(c); });
}

inline bool managed_output(const std::string& name) {
    return numbered_output(name) || name == "arrays.tsv" || name == "entries.tsv" ||
           name == "read_linked_groups.tsv" || name == "read_linked_edges.tsv" ||
           name == "read_linked_scaffolds.tsv" ||
           name == "read_linked_scaffold_edges.tsv" ||
           name == "read_linked_scaffold_components.tsv" ||
           name == "read_linked_scaffold_component_edges.tsv" ||
           name == "database_evidence.tsv" || name == "database_groups.tsv" ||
           name == "database_group_members.tsv" ||
           name == "element_confidence.tsv" ||
           name == "call_decisions.tsv" || name == "supported_entries.tsv" ||
           name == "supported_arrays.tsv" || name == "supported_spacers.fasta" ||
           name == "read_first_entries.tsv" || name == "read_first_spacers.fasta" ||
           name == "read_first_families.tsv" || name == "read_first_observations.tsv" ||
           name == "read_first_graph_boundary.tsv" || name == "read_first_variant_evidence.tsv" ||
           name == "read_first_family_recut.tsv" ||
           name == "repeat_haplotype_groups.tsv" ||
           name == "singleton_precision_pruned.tsv" ||
           name == "circular_unit_dedup.tsv" ||
           name == "directed_haplotype_paths.tsv" ||
           name == "directed_haplotype_edges.tsv" ||
           name == "cross_variant_scaffolds.tsv" ||
           name == "candidate_alleles.tsv" ||
           name == "read_recruitment.tsv" || name == "read_recruitment_candidates.tsv" ||
           name == "read_recruitment_original_rows.tsv" || name == "read_recruitment_context.tsv" ||
           name == "boundary_candidates.tsv" || name == "boundary_evidence.tsv" ||
           name == "local_boundary_graph.tsv" || name == "motif_cliff_boundary.tsv" ||
           name == "run_manifest.json";
}

// Stage every result before touching previous outputs. Roll back rename errors;
// unrelated files and cycle evidence are never included in the replacement.
// This is an exception-safe replacement, not an atomic multi-file transaction:
// concurrent writers/readers and process/power failure require external locking.
class Transaction {
    fs::path output_;
    fs::path scratch_;
    bool keep_recovery_ = false;
public:
    Transaction(const fs::path& output, const fs::path& input,
                const std::vector<std::string>& additional_inputs = {}) : output_(output) {
        std::vector<fs::path> inputs{input};
        for (const auto& path : additional_inputs) inputs.emplace_back(path);
        for (std::size_t i = 0; i < inputs.size(); ++i) {
            const auto resolved_input = fs::weakly_canonical(inputs[i]);
            if (resolved_input.parent_path() == fs::weakly_canonical(output_) &&
                managed_output(resolved_input.filename().string()))
                throw std::runtime_error(std::string(i ? "Read input" : "Cycles input") +
                    " conflicts with a managed result file: " + inputs[i].string());
        }
        for (const auto& entry : fs::directory_iterator(output_)) {
            if (!managed_output(entry.path().filename().string())) continue;
            if (fs::is_symlink(entry.symlink_status()) || !entry.is_regular_file())
                throw std::runtime_error("Cannot replace arrays file/result file: " + entry.path().string());
            for (std::size_t i = 0; i < inputs.size(); ++i)
                if (fs::equivalent(entry.path(), inputs[i]))
                    throw std::runtime_error(std::string(i ? "Read input" : "Cycles input") +
                        " aliases a managed result file: " + inputs[i].string());
        }
        const auto nonce = std::chrono::steady_clock::now().time_since_epoch().count();
        for (unsigned int attempt = 0; attempt < 100; ++attempt) {
            const auto candidate = output_ / (".Zuck-results-" +
                std::to_string(nonce) + "-" + std::to_string(attempt));
            if (!fs::create_directory(candidate)) continue;
            scratch_ = candidate;
            try {
                fs::create_directory(scratch_ / "new");
                fs::create_directory(scratch_ / "old");
            } catch (...) {
                std::error_code ignored;
                fs::remove_all(scratch_, ignored);
                throw;
            }
            return;
        }
        throw std::runtime_error("Cannot create result staging directory: " + output_.string());
    }
    Transaction(const Transaction&) = delete;
    Transaction& operator=(const Transaction&) = delete;
    ~Transaction() {
        if (!keep_recovery_) {
            std::error_code ignored;
            fs::remove_all(scratch_, ignored);
        }
    }
    fs::path directory() const { return scratch_ / "new"; }

    void commit() {
        std::vector<fs::path> old_names, new_names, backed_up, installed;
        for (const auto& entry : fs::directory_iterator(output_)) {
            if (!managed_output(entry.path().filename().string())) continue;
            if (fs::is_symlink(entry.symlink_status()) || !entry.is_regular_file())
                throw std::runtime_error("Cannot replace arrays file/result file: " + entry.path().string());
            old_names.push_back(entry.path().filename());
        }
        for (const auto& entry : fs::directory_iterator(directory()))
            new_names.push_back(entry.path().filename());
        std::sort(old_names.begin(), old_names.end());
        // Install the manifest last: its inventory/checksums describe the completed bundle.
        std::sort(new_names.begin(), new_names.end(), [](const auto& left, const auto& right) {
            if (left == "run_manifest.json") return false;
            if (right == "run_manifest.json") return true;
            return left < right;
        });
        backed_up.reserve(old_names.size());
        installed.reserve(new_names.size());
        try {
            for (const auto& name : old_names) {
                backed_up.push_back(name);
                try { fs::rename(output_ / name, scratch_ / "old" / name); }
                catch (...) { backed_up.pop_back(); throw; }
            }
            for (const auto& name : new_names) {
                installed.push_back(name);
                try { fs::rename(directory() / name, output_ / name); }
                catch (...) { installed.pop_back(); throw; }
            }
        } catch (const std::exception& error) {
            for (const auto& name : installed) {
                std::error_code ec;
                fs::remove(output_ / name, ec);
                keep_recovery_ |= static_cast<bool>(ec);
            }
            for (const auto& name : backed_up) {
                std::error_code ec;
                fs::rename(scratch_ / "old" / name, output_ / name, ec);
                keep_recovery_ |= static_cast<bool>(ec);
            }
            throw std::runtime_error(std::string("Result replacement failed: ") + error.what() +
                (keep_recovery_ ? "; recovery files retained at " + scratch_.string()
                                : "; previous results restored"));
        }
    }
};
} // namespace result_bundle
#endif
