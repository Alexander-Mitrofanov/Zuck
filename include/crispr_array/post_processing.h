#ifndef POST_PROCESSING_H
#define POST_PROCESSING_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <map>
#include <set>
#include <unordered_map>
#include <unordered_set>
#include <algorithm>
#include <cctype>
#include <numeric>
#include <filesystem>
#include <functional>
#include <chrono>
#include <iomanip>
#include <ctime>
#include <optional>
#include "spoa/spoa.hpp"
#include "core/settings.h"
#include "crispr_array/alternative_repeat_view.h"
#include "crispr_array/post_processing_helpers.h"
#include "crispr_array/short_spacer_filter.h"
#include "crispr_array/terminal_variant_rescue.h"
#include "crispr_array/read_linked_groups.h"
#include "crispr_array/array_database.h"
#include "crispr_array/database_groups.h"
#include "crispr_array/element_confidence.h"
#include "crispr_array/supported_calls.h"
#include "crispr_array/read_recruitment_output.h"
#include "crispr_array/read_first_output.h"
#include "crispr_array/local_boundary_graph.h"
#include "crispr_array/motif_cliff_boundary.h"
#include "crispr_array/directed_haplotype_paths.h"
#include "crispr_array/cross_variant_scaffolds.h"
#include "crispr_array/candidate_evidence.h"
#include "crispr_array/boundary_candidates.h"
#include "crispr_array/strand_equivalent_groups.h"
#include "io/run_manifest.h"

namespace fs = std::filesystem;

/**
 * PostProcessor: Aggregate cycle evidence and output CRISPR arrays.
 *
 * Pipeline:
 * 1-4: Glue, group by 23-mer, extend, and extract spacers. The optional
 *      two-sided mode reconstructs both repeat portions from circular cores.
 * 5:   Aggregate exact spacer alleles within compatible repeat families
 * 6:   Cluster similar repeats by edit distance, SPOA consensus (front repeat)
 * 7:   Per array: detect tail repeat in spacer suffixes, validate it matches
 *       a prefix of front consensus, SPOA consensus on tails, strip + prepend
 * 8:   Sanity filter on repeat/spacer dimensions
 * 9:   Output
 */
class PostProcessor {
private:
    Settings& settings;

    static constexpr int GROUP_KMER_SIZE = 23;
    static constexpr int TAIL_KMER_SIZE = 22;
    static constexpr int MAX_LINES_PER_FILE = 10000;

    // SPOA merge parameters
    static constexpr int MAX_EDIT_DIST = 4;

    // Tail detection parameters
    static constexpr int MAX_TAIL_SCAN = 40;
    static constexpr double TAIL_AGREE_THRESH = 0.6;
    static constexpr int TAIL_SUSTAIN_WINDOW = 3;
    static constexpr int MIN_TAIL_LEN = 5;
    // Boundary correction must never create a spacer shorter than the same
    // 20-nt floor used during extraction.
    static constexpr int MIN_SPACER_AFTER_TRIM = 20;

    // Tail validation: max edit distance between tail_consensus and front_consensus prefix
    static constexpr int MAX_TAIL_PREFIX_EDIT_DIST = 2;

    // Sanity filters
    static constexpr int MAX_REPEAT_LEN = 55;
    static constexpr int MIN_REPEAT_LEN = 20;
    static constexpr int MIN_SPACER_LEN = 20;
    static constexpr double MIN_MEDIAN_SPACER_REPEAT_RATIO = 0.5;
    struct ConsensusArray {
        std::string consensus;
        std::vector<std::pair<std::string, std::string>> entries;
    };
    std::vector<ConsensusArray> consensus_arrays;

    std::map<std::string, std::vector<std::string>> crispr_arrays;

    struct SpacerData {
        std::string spacer;
        std::string repeat;
        std::string cycle;
        int boundary_front = 0;
        int boundary_tail = 0;
    };

    static int edit_distance(const std::string& a, const std::string& b, int limit) {
        return Zuck::post_processing::bounded_edit_distance(a, b, limit);
    }

    static std::string revcomp(const std::string& s) {
        static const auto ct = [](){
            std::array<char,256> t{};
            t['A']='T'; t['T']='A'; t['C']='G'; t['G']='C';
            t['a']='T'; t['t']='A'; t['c']='G'; t['g']='C';
            t['N']='N'; t['n']='N';
            return t;
        }();
        std::string rc(s.size(), 'N');
        for (size_t i = 0; i < s.size(); ++i)
            rc[s.size()-1-i] = ct[(unsigned char)s[i]];
        return rc;
    }

    // Canonical form of a sequence: min(seq, RC(seq)).
    // Ensures that a sequence and its reverse complement map to the same key,
    // so the graph finding both strands does not double-count arrays.
    static std::string canonical_seq(const std::string& s) {
        std::string rc = revcomp(s);
        return (s <= rc) ? s : rc;
    }

    std::string spoa_consensus(const std::vector<std::pair<std::string, int>>& seq_weights) {
        if (seq_weights.empty()) return "";
        if (seq_weights.size() == 1) return seq_weights[0].first;

        auto alignment_engine = spoa::AlignmentEngine::Create(
            spoa::AlignmentType::kNW, 5, -4, -8, -6);

        spoa::Graph graph;

        for (const auto& [seq, weight] : seq_weights) {
            for (int w = 0; w < weight; ++w) {
                auto alignment = alignment_engine->Align(seq, graph);
                graph.AddAlignment(alignment, seq);
            }
        }

        return graph.GenerateConsensus();
    }

    /**
     * Detect tail repeat length by reverse majority voting on spacer suffixes.
     */
    int detect_tail_length(const std::vector<std::string>& spacers) {
        if (spacers.size() < 2) return 0;

        int min_len = (int)spacers[0].size();
        for (const auto& s : spacers)
            min_len = std::min(min_len, (int)s.size());

        int scan_limit = std::min(MAX_TAIL_SCAN, min_len - MIN_SPACER_AFTER_TRIM);
        if (scan_limit < MIN_TAIL_LEN) return 0;

        int tail_len = 0;
        int low_streak = 0;

        for (int rpos = 0; rpos < scan_limit; ++rpos) {
            int counts[4] = {0, 0, 0, 0};
            int total = 0;
            for (const auto& s : spacers) {
                int pos = (int)s.size() - 1 - rpos;
                switch (s[pos]) {
                    case 'A': counts[0]++; total++; break;
                    case 'C': counts[1]++; total++; break;
                    case 'G': counts[2]++; total++; break;
                    case 'T': counts[3]++; total++; break;
                }
            }
            if (total == 0) break;

            int best = *std::max_element(counts, counts + 4);
            double agreement = (double)best / total;

            if (agreement >= TAIL_AGREE_THRESH) {
                low_streak = 0;
                tail_len = rpos + 1;
            } else {
                low_streak++;
                if (low_streak >= TAIL_SUSTAIN_WINDOW) break;
            }
        }

        return (tail_len >= MIN_TAIL_LEN) ? tail_len : 0;
    }

    /**
     * Validate that tail_consensus approximately matches a prefix of front_consensus.
     * This is the key biological constraint: the tail is the beginning of the next
     * repeat occurrence in the cycle, so it must look like the front repeat.
     */
    bool validate_tail_matches_front(const std::string& tail, const std::string& front) {
        if (tail.empty() || front.empty()) return false;

        int check_len = std::min((int)tail.size(), (int)front.size());
        std::string tail_portion = tail.substr(tail.size() - check_len);
        std::string front_prefix = front.substr(0, check_len);

        int dist = edit_distance(tail_portion, front_prefix, MAX_TAIL_PREFIX_EDIT_DIST + 1);
        return dist <= MAX_TAIL_PREFIX_EDIT_DIST;
    }

    /**
     * Compute median of a vector of ints.
     */
    double median(std::vector<int>& vals) {
        if (vals.empty()) return 0;
        std::sort(vals.begin(), vals.end());
        int n = (int)vals.size();
        if (n % 2 == 0) return (vals[n/2 - 1] + vals[n/2]) / 2.0;
        return vals[n/2];
    }

    /**
     * Frozen release rotation. The legacy primary pipeline deliberately
     * rotates the complete glued record, including its 22-base closure, and
     * keeps the first maximum-frequency 23-mer. Supplemental views must not
     * change this primary boundary decision.
     */
    std::string freq_rotation(
        const std::string& sequence,
        const std::unordered_map<std::string, int>& frequencies) const {
        const int length = static_cast<int>(sequence.size());
        if (length < GROUP_KMER_SIZE) return sequence;

        int best_position = 0;
        int best_frequency = 0;
        for (int position = 0;
             position + GROUP_KMER_SIZE <= length; ++position) {
            const auto found = frequencies.find(
                sequence.substr(position, GROUP_KMER_SIZE));
            const int frequency =
                found == frequencies.end() ? 0 : found->second;
            if (frequency > best_frequency) {
                best_frequency = frequency;
                best_position = position;
            }
        }
        if (best_position == 0) return sequence;
        return sequence.substr(best_position) +
               sequence.substr(0, best_position);
    }

    // A cycle's closure repeats its first 22 bases only for serialization.
    // Rotate the circular core and break equal-frequency anchor-sequence ties
    // lexically rather than selecting a different anchor from the record start.
    std::string circular_freq_rotation(
        const std::string& core,
        const std::unordered_map<std::string, int>& frequencies) const {
        if (core.size() < GROUP_KMER_SIZE) return core;
        const std::string doubled = core + core;
        std::size_t best_position = 0;
        int best_frequency = 0;
        std::string best_kmer;
        for (std::size_t position = 0; position < core.size(); ++position) {
            const std::string kmer = doubled.substr(position, GROUP_KMER_SIZE);
            const auto found = frequencies.find(kmer);
            const int frequency = found == frequencies.end() ? 0 : found->second;
            if (frequency > best_frequency ||
                (frequency == best_frequency &&
                 (best_kmer.empty() || kmer < best_kmer))) {
                best_frequency = frequency;
                best_position = position;
                best_kmer = kmer;
            }
        }
        return core.substr(best_position) + core.substr(0, best_position);
    }

    std::vector<SpacerData> extract_two_sided_group_single(
        const std::string& anchor,
        const std::vector<std::string>& group_cycles) {
        const std::set<std::string> distinct(group_cycles.begin(), group_cycles.end());
        if (distinct.size() < 2) return {};

        int shortest_core = static_cast<int>(distinct.begin()->size());
        for (const auto& core : distinct)
            shortest_core = std::min(shortest_core, static_cast<int>(core.size()));

        // Keep an exact shared-boundary fallback for rows that do not support
        // the majority repeat. Both sides must leave a complete spacer.
        int front_length = GROUP_KMER_SIZE;
        while (front_length < MAX_REPEAT_LEN &&
               front_length < shortest_core - MIN_SPACER_LEN) {
            std::set<char> bases;
            for (const auto& core : distinct) bases.insert(core[front_length]);
            if (bases.size() != 1) break;
            ++front_length;
        }
        int tail_length = 0;
        while (front_length + tail_length < MAX_REPEAT_LEN &&
               front_length + tail_length < shortest_core - MIN_SPACER_LEN) {
            std::set<char> bases;
            for (const auto& core : distinct)
                bases.insert(core[core.size() - 1 - tail_length]);
            if (bases.size() != 1) break;
            ++tail_length;
        }

        // Count distinct rotated cores, rather than repeated traversal records,
        // for each boundary vote. The repeat prefix before the anchor and its
        // suffix after the anchor are different portions of the same repeat;
        // there is no requirement that their sequences match each other.
        auto majority_base = [&](int position, bool from_end) {
            std::map<char, std::size_t> counts;
            for (const auto& core : distinct)
                ++counts[from_end ? core[core.size() - 1 - position] : core[position]];
            char best = 'N';
            std::size_t support = 0;
            for (const auto& [base, count] : counts) {
                if (count > support) {
                    best = base;
                    support = count;
                }
            }
            return static_cast<double>(support) / distinct.size() >= TAIL_AGREE_THRESH
                ? best : '\0';
        };

        std::string majority_front = anchor;
        std::string majority_tail;
        while (static_cast<int>(majority_front.size()) < MAX_REPEAT_LEN &&
               static_cast<int>(majority_front.size()) < shortest_core - MIN_SPACER_LEN) {
            const char next = majority_base(static_cast<int>(majority_front.size()), false);
            if (!next) break;
            majority_front.push_back(next);
        }
        while (static_cast<int>(majority_front.size() + majority_tail.size()) < MAX_REPEAT_LEN &&
               static_cast<int>(majority_front.size() + majority_tail.size()) <
                   shortest_core - MIN_SPACER_LEN) {
            const char next = majority_base(static_cast<int>(majority_tail.size()), true);
            if (!next) break;
            majority_tail.insert(majority_tail.begin(), next);
        }
        const std::string majority_repeat = majority_tail + majority_front;

        std::vector<SpacerData> extracted;
        for (const auto& core : group_cycles) {
            int row_front = front_length;
            int row_tail = tail_length;
            const std::string proposed_repeat =
                core.substr(core.size() - majority_tail.size()) +
                core.substr(0, majority_front.size());
            if (edit_distance(proposed_repeat, majority_repeat, MAX_TAIL_PREFIX_EDIT_DIST) <=
                MAX_TAIL_PREFIX_EDIT_DIST) {
                row_front = static_cast<int>(majority_front.size());
                row_tail = static_cast<int>(majority_tail.size());
            }
            if (static_cast<int>(core.size()) < row_front + row_tail + MIN_SPACER_LEN)
                continue;
            extracted.push_back({
                core.substr(row_front, core.size() - row_front - row_tail),
                core.substr(core.size() - row_tail) + core.substr(0, row_front),
                core, row_front, row_tail});
        }
        return extracted;
    }

    // Opt-in mixture inference: a shared anchor need not imply one repeat.
    // Linked flank signatures propose branches; per-side provenance, distinct
    // core support, and the existing repeat-family radius constrain adoption.
    // This is not a calibrated significance test or evidence of unique loci.
    std::vector<SpacerData> extract_two_sided_group(
        const std::string& anchor,
        const std::vector<std::string>& group_cycles,
        std::size_t refinement_depth = 0) {
        auto baseline = extract_two_sided_group_single(anchor, group_cycles);
        if (!settings.branch_aware_boundaries) return baseline;
        // Two partition levels suffice to examine a nested family mixture.
        // The bound is independent of the number of input cycles.
        if (settings.repeat_mixture_refinement && refinement_depth >= 2) return baseline;
        const std::set<std::string> distinct(group_cycles.begin(), group_cycles.end());
        if (distinct.size() < 6 || baseline.empty()) return baseline;
        int shortest = static_cast<int>(distinct.begin()->size());
        for (const auto& core : distinct)
            shortest = std::min(shortest, static_cast<int>(core.size()));
        int front = GROUP_KMER_SIZE;
        int tail = 0;
        while (front < MAX_REPEAT_LEN && front < shortest - MIN_SPACER_LEN) {
            std::set<char> bases;
            for (const auto& core : distinct) bases.insert(core[front]);
            if (bases.size() != 1) break;
            ++front;
        }
        while (front + tail < MAX_REPEAT_LEN &&
               front + tail < shortest - MIN_SPACER_LEN) {
            std::set<char> bases;
            for (const auto& core : distinct)
                bases.insert(core[core.size() - 1 - tail]);
            if (bases.size() != 1) break;
            ++tail;
        }
        const int minimum_signature = settings.repeat_mixture_refinement ? 5 : 6;
        if (front + tail + minimum_signature > MAX_REPEAT_LEN ||
            front + tail + minimum_signature > shortest - MIN_SPACER_LEN) return baseline;
        std::map<std::string, SpacerData> original;
        for (const auto& observation : baseline) original[observation.cycle] = observation;
        std::map<std::string, SpacerData> best;
        std::size_t best_gain = 0;
        // Keep historical shapes and tie order unchanged unless refinement is
        // explicitly enabled. Five bases can resolve a short remaining repeat
        // allele block without conditioning on a sixth base in the spacer.
        std::vector<std::pair<int, int>> shapes{{6, 0}, {0, 6}, {3, 3}};
        if (settings.repeat_mixture_refinement) {
            shapes.emplace_back(5, 0);
            shapes.emplace_back(0, 5);
        }
        for (const auto& [f, t] : shapes) {
            if (front + tail + f + t > MAX_REPEAT_LEN ||
                front + tail + f + t > shortest - MIN_SPACER_LEN) continue;
            std::map<std::string, std::vector<std::string>> partitions;
            for (const auto& core : distinct)
                partitions[core.substr(front, f) +
                           core.substr(core.size() - tail - t, t)].push_back(core);
            std::map<std::string, SpacerData> replacement;
            std::size_t covered = 0;
            std::size_t branches = 0;
            std::size_t gain = 0;
            bool invalid = false;
            std::vector<std::string> branch_repeats;
            for (const auto& partition : partitions) {
                const auto& members = partition.second;
                if (members.size() < 3) continue;
                ++branches;
                if (branches > 4) { invalid = true; break; }
                auto candidate = settings.repeat_mixture_refinement
                    ? extract_two_sided_group(anchor, members, refinement_depth + 1)
                    : extract_two_sided_group_single(anchor, members);
                std::map<std::string, std::size_t> repeat_counts;
                for (const auto& row : candidate) {
                    // The signature must lie within both inferred repeat
                    // portions, not merely within their combined length.
                    if (row.boundary_front < front + f || row.boundary_tail < tail + t) {
                        invalid = true;
                        break;
                    }
                    replacement[row.cycle] = row;
                    ++repeat_counts[row.repeat];
                    const auto old = original.find(row.cycle);
                    if (old != original.end() &&
                        row.repeat.size() > old->second.repeat.size() + MAX_EDIT_DIST)
                        gain += row.repeat.size() - old->second.repeat.size();
                }
                if (invalid) break;
                const auto modal = std::max_element(repeat_counts.begin(), repeat_counts.end(),
                    [](const auto& a, const auto& b) { return a.second < b.second; });
                if (modal == repeat_counts.end()) { invalid = true; break; }
                // The historical rule rejects all near-identical branches.
                // Refinement also admits a block of >=3 substitutions, which
                // can be a repeat allele rather than one/two random spacer
                // endpoint bases. Conserved spacer termini can remain
                // observationally ambiguous; this is an experimental model.
                const int minimum_separation = settings.repeat_mixture_refinement ? 2 : MAX_EDIT_DIST;
                for (const auto& repeat : branch_repeats)
                    if (edit_distance(repeat, modal->first, minimum_separation) <= minimum_separation) {
                        invalid = true;
                        break;
                    }
                if (invalid) break;
                branch_repeats.push_back(modal->first);
                covered += members.size();
            }
            if (invalid || branches < 2 || covered * 5 < distinct.size() * 4 || gain <= best_gain)
                continue;
            best_gain = gain;
            best = std::move(replacement);
        }
        for (auto& row : baseline) {
            const auto found = best.find(row.cycle);
            // Recover substantially truncated repeats without shifting an
            // already supported boundary by a few chance-consensus bases.
            if (found != best.end() &&
                found->second.repeat.size() > row.repeat.size() + MAX_EDIT_DIST)
                row = found->second;
        }
        return baseline;
    }

    void write_structured_results(
        const fs::path& directory, const std::string& input_path,
        const std::string& input_checksum,
        const std::map<std::string, std::size_t>& counts,
        const std::vector<std::size_t>& supplemental_parents, double seconds) const {
        const auto arrays_path = directory / "arrays.tsv";
        const auto entries_path = directory / "entries.tsv";
        std::ofstream arrays(arrays_path), entries(entries_path);
        if (!arrays || !entries)
            throw std::runtime_error("Cannot open structured result tables: " + directory.string());
        arrays << "array_id\tkind\tparent_array_id\trepeat_consensus\tentry_count"
                  "\tdistinct_spacer_count\tdistinct_canonical_spacer_count\trepeat_allele_count\n";
        entries << "array_id\tentry_index\trepeat_variant\tspacer\n";
        const auto primary_count = counts.at("primary_arrays");
        for (std::size_t index = 0; index < consensus_arrays.size(); ++index) {
            const auto& array = consensus_arrays[index];
            const auto id = "Array_" + std::to_string(index + 1);
            const bool primary = index < primary_count;
            std::set<std::string> spacers, canonical_spacers, repeats;
            for (std::size_t row = 0; row < array.entries.size(); ++row) {
                const auto& [repeat, spacer] = array.entries[row];
                entries << id << '\t' << row + 1 << '\t' << repeat << '\t' << spacer << '\n';
                spacers.insert(spacer);
                canonical_spacers.insert(canonical_seq(spacer));
                repeats.insert(repeat);
            }
            const bool database_group = primary && counts.count("database_de_novo_groups") &&
                index >= counts.at("database_de_novo_groups");
            const bool boundary_group = primary && counts.count("boundary_previous_groups") &&
                index >= counts.at("boundary_previous_groups");
            arrays << id << '\t' << (boundary_group ? "boundary_supported" : database_group ? "database_supported" :
                primary ? "primary" : "experimental_alternative_view") << '\t';
            if (!primary) arrays << "Array_" << supplemental_parents.at(index - primary_count) + 1;
            arrays << '\t' << array.consensus << '\t' << array.entries.size() << '\t'
                   << spacers.size() << '\t' << canonical_spacers.size() << '\t' << repeats.size() << '\n';
        }
        result_bundle::close_checked(arrays, arrays_path);
        result_bundle::close_checked(entries, entries_path);
        result_bundle::write_manifest(directory, input_path, input_checksum, settings, counts, {
            {"group_kmer_size", GROUP_KMER_SIZE}, {"closure_length", TAIL_KMER_SIZE},
            {"repeat_cluster_max_edit_distance", MAX_EDIT_DIST},
            {"tail_scan_limit", MAX_TAIL_SCAN}, {"tail_agreement", TAIL_AGREE_THRESH},
            {"tail_sustain_window", TAIL_SUSTAIN_WINDOW}, {"tail_minimum_length", MIN_TAIL_LEN},
            {"spacer_minimum_after_trim", MIN_SPACER_AFTER_TRIM},
            {"tail_prefix_max_edit_distance", MAX_TAIL_PREFIX_EDIT_DIST},
            {"repeat_minimum_length", MIN_REPEAT_LEN}, {"repeat_maximum_length", MAX_REPEAT_LEN},
            {"spacer_minimum_length", MIN_SPACER_LEN},
            {"terminal_rescue_minimum_observations", Zuck::TERMINAL_RESCUE_MIN_OBSERVATIONS},
            {"terminal_rescue_maximum_component_size", Zuck::TERMINAL_RESCUE_MAX_COMPONENT_SIZE},
            {"terminal_rescue_length_difference", Zuck::TERMINAL_RESCUE_LENGTH_DIFFERENCE},
            {"short_spacer_anchor_length", Zuck::SHORT_SPACER_ANCHOR_LENGTH},
            {"short_spacer_max_terminal_extension", Zuck::SHORT_SPACER_MAX_TERMINAL_EXTENSION},
            {"read_linked_minimum_anchor", Zuck::read_linked::Config{}.minimum_anchor},
            {"read_linked_minimum_quality", Zuck::read_linked::Config{}.minimum_quality},
            {"read_linked_maximum_expected_errors", Zuck::read_linked::Config{}.maximum_expected_errors},
            {"read_linked_minimum_sequence_units", Zuck::read_linked::Config{}.minimum_sequence_units},
            {"read_linked_minimum_start_offsets", Zuck::read_linked::Config{}.minimum_start_offsets},
            {"read_linked_minimum_component_spacers", Zuck::read_linked::Config{}.minimum_component_spacers},
            {"read_linked_minimum_repeat_length", Zuck::read_linked::Config{}.minimum_repeat},
            {"read_linked_maximum_repeat_length", Zuck::read_linked::Config{}.maximum_repeat},
            {"database_minimum_bona_fide_spacers", Zuck::array_database::MIN_BONA_FIDE_SPACERS},
            {"database_minimum_possible_spacers", Zuck::array_database::MIN_POSSIBLE_SPACERS},
            {"database_minimum_bona_fide_confidence", Zuck::array_database::MIN_BONA_FIDE_CONFIDENCE},
            {"database_minimum_possible_confidence", Zuck::array_database::MIN_POSSIBLE_CONFIDENCE},
            {"database_maximum_core_interpretations", Zuck::array_database::MAX_CORE_INTERPRETATIONS},
            {"database_group_minimum_bona_fide_spacers", Zuck::database_groups::MIN_BONA_FIDE_SPACERS},
            {"database_group_minimum_possible_spacers", Zuck::database_groups::MIN_POSSIBLE_SPACERS},
            {"database_group_maximum_views_per_parent", Zuck::database_groups::MAX_VIEWS_PER_PARENT},
            {"database_read_next_repeat_anchor", Zuck::database_reads::NEXT_REPEAT_ANCHOR},
            {"database_read_minimum_quality", Zuck::database_reads::MINIMUM_QUALITY},
            {"database_read_maximum_expected_errors", Zuck::database_reads::MAXIMUM_EXPECTED_ERRORS},
            {"database_read_minimum_sequence_units", Zuck::database_reads::MINIMUM_SEQUENCE_UNITS},
            {"database_read_minimum_start_offsets", Zuck::database_reads::MINIMUM_START_OFFSETS},
            {"minimum_median_spacer_repeat_ratio", MIN_MEDIAN_SPACER_REPEAT_RATIO}
        }, seconds);
    }

    // Parse replay records, enforce closure, then select deterministic anchors.
    std::vector<std::string> load_rotated_cycles(const std::string& input_path) {
        std::vector<std::string> cycles;
        {
            // First pass: read and glue k-mers into raw sequences
            std::vector<std::string> raw_cycles;
            {
                std::ifstream in(input_path);
                if (!in.is_open()) {
                    throw std::runtime_error("Cannot open cycles file: " + input_path);
                }
                std::string line;
                std::size_t line_number = 0;
                while (std::getline(in, line)) {
                    ++line_number;
                    if (!line.empty()) {
                        line.erase(0, line.find_first_not_of(" \t\r\n"));
                        line.erase(line.find_last_not_of(" \t\r\n") + 1);
                        if (!line.empty() && line.front() == '#') continue;
                        if (!line.empty()) {
                            std::string glued;
                            try {
                                glued = Zuck::post_processing::glue_kmers(line);
                            } catch (const std::invalid_argument&) {
                                throw std::runtime_error(
                                    "Invalid compact/legacy cycle record at " +
                                    input_path + ":" + std::to_string(line_number));
                            }
                            if (settings.two_sided_boundaries && !glued.empty()) {
                                if (glued.size() < 2 * TAIL_KMER_SIZE ||
                                    glued.compare(0, TAIL_KMER_SIZE, glued,
                                                  glued.size() - TAIL_KMER_SIZE,
                                                  TAIL_KMER_SIZE) != 0)
                                    throw std::runtime_error(
                                        "Two-sided boundaries require an exact cycle closure at " +
                                        input_path + ":" + std::to_string(line_number));
                                glued.resize(glued.size() - TAIL_KMER_SIZE);
                            }
                            if (!glued.empty()) raw_cycles.push_back(std::move(glued));
                        }
                    }
                }
                if (in.bad() || (in.fail() && !in.eof()))
                    throw std::runtime_error("Error reading cycles file: " + input_path);
            }

            // Legacy mode counts linear 23-mers on the complete glued record,
            // including its closure. Two-sided mode counts each circular-core
            // position and applies a deterministic anchor tie break.
            std::unordered_map<std::string,int> kmer_freq;
            kmer_freq.reserve(raw_cycles.size() * 4);
            for (const auto& seq : raw_cycles) {
                const std::string counted = settings.two_sided_boundaries ? seq + seq : seq;
                const int positions = settings.two_sided_boundaries
                    ? static_cast<int>(seq.size())
                    : static_cast<int>(seq.size()) - GROUP_KMER_SIZE + 1;
                for (int position = 0; position < positions; ++position)
                    ++kmer_freq[counted.substr(position, GROUP_KMER_SIZE)];
            }

            cycles = std::move(raw_cycles);
            for (auto& seq : cycles)
                seq = settings.two_sided_boundaries
                    ? circular_freq_rotation(seq, kmer_freq)
                    : freq_rotation(seq, kmer_freq);
        }
        return cycles;
    }

    // Group rotated cores by anchor and infer observed repeat/spacer pairs.
    std::vector<SpacerData> extract_spacer_observations(
        std::vector<std::string> cycles, std::size_t& anchor_group_count) {
        // ==========================================================
        // Group by the selected 23-mer anchor
        // ==========================================================
        std::unordered_map<std::string, std::vector<std::string>> groups;
        for (auto& seq : cycles) {
            if ((int)seq.size() >= GROUP_KMER_SIZE)
                groups[seq.substr(0, GROUP_KMER_SIZE)].push_back(std::move(seq));
        }
        std::vector<std::string>().swap(cycles);
        anchor_group_count = groups.size();
        std::cout << "  ▸ Grouped into " << groups.size() << " repeat groups" << std::endl;

        // ==========================================================
        // Infer boundaries and extract spacers
        // ==========================================================
        std::vector<SpacerData> spacer_data;

        for (const auto& [kmer, group_cycles] : groups) {
            if (settings.two_sided_boundaries) {
                for (auto& observation : extract_two_sided_group(kmer, group_cycles))
                    spacer_data.push_back(std::move(observation));
                continue;
            }
            int t = 0;
            while (true) {
                std::unordered_set<char> chars;
                bool valid = true;
                for (const auto& c : group_cycles) {
                    int pos = GROUP_KMER_SIZE + t;
                    if (pos >= (int)c.size() - TAIL_KMER_SIZE) {
                        valid = false;
                        break;
                    }
                    chars.insert(c[pos]);
                }
                if (!valid || chars.size() != 1) break;
                ++t;
            }

            int repeat_len = GROUP_KMER_SIZE + t;

            for (const auto& c : group_cycles) {
                if ((int)c.size() > repeat_len + TAIL_KMER_SIZE) {
                    std::string repeat = c.substr(0, repeat_len);
                    std::string spacer = c.substr(repeat_len, c.size() - repeat_len - TAIL_KMER_SIZE);
                    if ((int)spacer.size() >= MIN_SPACER_LEN) {
                        spacer_data.push_back({spacer, repeat, c});
                    }
                }
            }
        }

        std::cout << "  ▸ Extracted " << spacer_data.size() << " spacers" << std::endl;
        return spacer_data;
    }

    struct RepeatFamilies {
        std::vector<std::string> repeats;
        std::vector<std::vector<int>> members;
    };

    // Every member of a repeat family must lie within the edit radius of every
    // other member. Support and lexical ties make greedy assignment reproducible.
    RepeatFamilies cluster_repeat_families(
        const std::map<std::string, std::map<std::string, std::size_t>>& repeat_to_spacer_support) {
        std::vector<std::string> all_repeats;
        for (const auto& [rep, _] : repeat_to_spacer_support)
            all_repeats.push_back(rep);

        int nr = (int)all_repeats.size();

        std::vector<std::size_t> repeat_raw_support(nr, 0);
        std::vector<std::size_t> repeat_distinct_spacers(nr, 0);
        for (int i = 0; i < nr; ++i) {
            const auto& spacers = repeat_to_spacer_support.at(all_repeats[i]);
            repeat_distinct_spacers[i] = spacers.size();
            for (const auto& [spacer, support] : spacers)
                repeat_raw_support[i] += support;
        }

        std::vector<int> repeat_order(nr);
        std::iota(repeat_order.begin(), repeat_order.end(), 0);
        std::sort(repeat_order.begin(), repeat_order.end(), [&](int lhs, int rhs) {
            if (repeat_raw_support[lhs] != repeat_raw_support[rhs])
                return repeat_raw_support[lhs] > repeat_raw_support[rhs];
            if (repeat_distinct_spacers[lhs] != repeat_distinct_spacers[rhs])
                return repeat_distinct_spacers[lhs] > repeat_distinct_spacers[rhs];
            return all_repeats[lhs] < all_repeats[rhs];
        });

        struct RepeatCluster {
            std::vector<int> members;
            std::size_t aggregate_support = 0;
        };
        std::vector<RepeatCluster> bounded_clusters;

        for (int repeat_idx : repeat_order) {
            std::size_t best_cluster = bounded_clusters.size();
            int best_max_distance = MAX_EDIT_DIST + 1;
            std::size_t best_cluster_support = 0;
            std::string best_anchor;

            for (std::size_t cluster_idx = 0;
                 cluster_idx < bounded_clusters.size(); ++cluster_idx) {
                int max_distance = 0;
                bool compatible = true;
                for (int member : bounded_clusters[cluster_idx].members) {
                    int distance = edit_distance(all_repeats[repeat_idx],
                                                 all_repeats[member],
                                                 MAX_EDIT_DIST);
                    if (distance > MAX_EDIT_DIST) {
                        compatible = false;
                        break;
                    }
                    max_distance = std::max(max_distance, distance);
                }
                if (!compatible) continue;

                const std::string& anchor =
                    all_repeats[bounded_clusters[cluster_idx].members.front()];
                if (max_distance < best_max_distance ||
                    (max_distance == best_max_distance &&
                     (bounded_clusters[cluster_idx].aggregate_support >
                          best_cluster_support ||
                      (bounded_clusters[cluster_idx].aggregate_support ==
                           best_cluster_support &&
                       (best_cluster == bounded_clusters.size() ||
                        anchor < best_anchor))))) {
                    best_cluster = cluster_idx;
                    best_max_distance = max_distance;
                    best_cluster_support =
                        bounded_clusters[cluster_idx].aggregate_support;
                    best_anchor = anchor;
                }
            }

            if (best_cluster == bounded_clusters.size()) {
                bounded_clusters.push_back(
                    {{repeat_idx}, repeat_raw_support[repeat_idx]});
            } else {
                bounded_clusters[best_cluster].members.push_back(repeat_idx);
                bounded_clusters[best_cluster].aggregate_support +=
                    repeat_raw_support[repeat_idx];
            }
        }

        std::sort(bounded_clusters.begin(), bounded_clusters.end(),
                  [&](const RepeatCluster& lhs, const RepeatCluster& rhs) {
                      return all_repeats[lhs.members.front()] <
                             all_repeats[rhs.members.front()];
                  });

        std::vector<std::vector<int>> repeat_clusters;
        repeat_clusters.reserve(bounded_clusters.size());
        for (auto& cluster : bounded_clusters)
            repeat_clusters.push_back(std::move(cluster.members));

        return {std::move(all_repeats), std::move(repeat_clusters)};
    }

    // Optional views are computed only when requested. They reuse accepted
    // rows and never alter the primary groups or their identifiers.
    std::vector<std::size_t> append_alternative_views(
        std::size_t primary_array_count, std::size_t primary_row_count,
        const Zuck::supported_calls::RecruitmentEvidence& recruitment = {}) {
        if (!settings.experimental_alternative_views) return {};
        std::set<std::pair<std::size_t, std::size_t>> evidence_only_rows;
        for (const auto& [row, proof] : recruitment) if (proof.added_row) {
            if (row.first >= primary_array_count || row.second == 0 ||
                row.second > consensus_arrays[row.first].entries.size())
                throw std::invalid_argument("Supplemental view references an absent recruited row");
            evidence_only_rows.insert(row);
        }
        if (evidence_only_rows.size() > primary_row_count)
            throw std::logic_error("Recruited rows exceed primary row count");
        const auto eligible_primary_rows = primary_row_count - evidence_only_rows.size();
        using CanonicalSpacerMap = std::map<std::string, std::string>;
        std::vector<std::set<std::string>> parent_spacers(
            primary_array_count);
        std::map<
            std::string,
            std::map<std::size_t, CanonicalSpacerMap>>
            allele_observations;
        std::map<std::string, std::set<std::size_t>>
            canonical_repeat_parents;

        for (std::size_t parent = 0;
             parent < primary_array_count; ++parent) {
            const auto& array = consensus_arrays[parent];
            canonical_repeat_parents[
                alternative_repeat_view::canonical_sequence(
                    array.consensus)].insert(parent);
            for (std::size_t entry = 0; entry < array.entries.size(); ++entry) {
                if (evidence_only_rows.count({parent, entry + 1})) continue;
                const auto& [repeat, spacer] = array.entries[entry];
                canonical_repeat_parents[
                    alternative_repeat_view::canonical_sequence(
                        repeat)].insert(parent);
                const std::string canonical =
                    alternative_repeat_view::canonical_spacer(spacer);
                parent_spacers[parent].insert(canonical);
                if (repeat == array.consensus) continue;

                // Retain an exact observed orientation. If both orientations
                // were accepted upstream, use the lexical observation rather
                // than synthesizing a new sequence.
                auto& observations =
                    allele_observations[repeat][parent];
                if (spacer.size() <
                    static_cast<std::size_t>(MIN_SPACER_AFTER_TRIM))
                    continue;
                auto [found, inserted] =
                    observations.emplace(canonical, spacer);
                if (!inserted && spacer < found->second)
                    found->second = spacer;
            }
        }

        struct SupplementalView {
            std::size_t parent = 0;
            std::string consensus;
            std::vector<std::pair<std::string, std::string>> entries;
            alternative_repeat_view::CyclicMatch match;
        };
        std::vector<SupplementalView> eligible_views;

        for (const auto& [allele, observations_by_parent] :
             allele_observations) {
            // A repeat allele shared by multiple accepted parents is not a
            // label-free discriminator and cannot create a supplemental view.
            if (observations_by_parent.size() != 1) continue;
            const auto& [parent, observations] =
                *observations_by_parent.begin();
            const auto& parent_array = consensus_arrays[parent];

            const auto canonical_parents = canonical_repeat_parents.find(
                alternative_repeat_view::canonical_sequence(allele));
            if (canonical_parents == canonical_repeat_parents.end() ||
                canonical_parents->second.size() != 1 ||
                !canonical_parents->second.count(parent))
                continue;

            if (!alternative_repeat_view::eligible_allele(
                    allele, parent_array.consensus, observations.size(),
                    parent_spacers[parent].size()))
                continue;

            SupplementalView view;
            view.parent = parent;
            view.consensus = allele;
            view.match = alternative_repeat_view::oriented_cyclic_match(
                allele, parent_array.consensus);
            std::set<std::string> emitted_spacers;
            for (std::size_t entry = 0; entry < parent_array.entries.size(); ++entry) {
                if (evidence_only_rows.count({parent, entry + 1})) continue;
                const auto& [observed_repeat, observed_spacer] = parent_array.entries[entry];
                if (observed_repeat != allele) continue;
                const std::string canonical =
                    alternative_repeat_view::canonical_spacer(
                        observed_spacer);
                if (!observations.count(canonical) ||
                    !emitted_spacers.insert(canonical).second)
                    continue;
                view.entries.emplace_back(allele, observed_spacer);
            }
            eligible_views.push_back(std::move(view));
        }

        std::sort(eligible_views.begin(), eligible_views.end(),
                  [](const SupplementalView& lhs,
                     const SupplementalView& rhs) {
                      if (lhs.parent != rhs.parent)
                          return lhs.parent < rhs.parent;
                      return lhs.consensus < rhs.consensus;
                  });

        // Integer floors make the caps strict for every primary output size.
        const std::size_t supplemental_array_budget =
            primary_array_count / 100;        // <= 1.00%
        const std::size_t supplemental_row_budget =
            eligible_primary_rows / 400;     // <= 0.25% of original eligible rows
        std::size_t supplemental_rows = 0;
        std::size_t supplemental_arrays = 0;
        std::vector<std::size_t> supplemental_parents;
        std::set<std::string> supplemental_spacers;

        for (const auto& view : eligible_views) {
            if (!settings.experimental_alternative_views) break;
            if (supplemental_arrays >= supplemental_array_budget) break;
            if (supplemental_rows + view.entries.size() >
                supplemental_row_budget)
                continue; // Never truncate a biologically defined subset.

            bool duplicates_another_view = false;
            for (const auto& [repeat, spacer] : view.entries) {
                (void)repeat;
                if (supplemental_spacers.count(
                        alternative_repeat_view::canonical_spacer(spacer))) {
                    duplicates_another_view = true;
                    break;
                }
            }
            if (duplicates_another_view) continue;

            consensus_arrays.push_back(
                {view.consensus, view.entries});
            supplemental_parents.push_back(view.parent);
            ++supplemental_arrays;
            supplemental_rows += view.entries.size();
            for (const auto& [repeat, spacer] : view.entries) {
                (void)repeat;
                supplemental_spacers.insert(
                    alternative_repeat_view::canonical_spacer(spacer));
            }

            std::ostringstream provenance;
            provenance << "  ▸ Supplemental Array_"
                       << consensus_arrays.size() << " <- Array_"
                       << (view.parent + 1)
                       << "  allele=" << view.consensus
                       << "  spacers=" << view.entries.size()
                       << "  subset=" << view.entries.size() << "/"
                       << parent_spacers[view.parent].size()
                       << "  identity=" << std::fixed
                       << std::setprecision(6) << view.match.identity
                       << "  coverage=" << view.match.reciprocal_coverage;
            std::cout << provenance.str() << std::endl;
        }

        if (settings.experimental_alternative_views) {
            std::cout << "  ▸ Experimental alternative-repeat sidecar: "
                      << supplemental_arrays << " arrays / "
                      << supplemental_rows << " rows"
                      << " (caps " << supplemental_array_budget << " / "
                      << supplemental_row_budget << ")" << std::endl;
        }

        return supplemental_parents;
    }

    // The human-readable files and structured tables share the same array
    // order. Supplemental views start in a separate numbered file.
    int write_text_results(const fs::path& result_directory,
                           const std::string& output_dir,
                           std::size_t primary_array_count,
                           std::size_t primary_row_count,
                           std::size_t supplemental_arrays,
                           std::size_t supplemental_rows) const {
        // ==========================================================
        // Step 8: Output
        // ==========================================================
        // Generate timestamp string for file header
        auto now = std::chrono::system_clock::now();
        auto now_t = std::chrono::system_clock::to_time_t(now);
        std::ostringstream ts;
        ts << std::put_time(std::localtime(&now_t), "%Y-%m-%d %H:%M:%S");
        std::string timestamp_str = ts.str();

        int file_index = 1;
        int current_lines = 0;
        std::ofstream out;
        std::string output_filename;

        auto check_output = [&]() {
            if (!out)
                throw std::runtime_error("Error writing arrays file: " + output_filename);
        };
        auto close_output = [&]() {
            if (!out.is_open()) return;
            check_output();
            out.close();
            check_output();
        };

        auto open_new_file = [&](std::size_t header_arrays,
                                 std::size_t header_spacers,
                                 bool experimental_sidecar = false) {
            close_output();
            output_filename = (result_directory / ("CRISPR_Arrays_" +
                              std::to_string(file_index) + ".txt")).string();
            out.open(output_filename);
            if (!out.is_open() || !out)
                throw std::runtime_error("Cannot open arrays file for writing: " + output_filename);
            ++file_index;
            current_lines = 0;
            // Write file header
            out << "# Zuck — CRISPR Array Output\n";
            out << "# Generated : " << timestamp_str << "\n";
            out << "# Arrays    : " << header_arrays << "\n";
            out << "# Spacers   : " << header_spacers << "\n";
            if (settings.branch_aware_boundaries)
                out << "# Boundary  : branch-aware (experimental)\n";
            else if (settings.two_sided_boundaries)
                out << "# Boundary  : two-sided\n";
            if (settings.repeat_mixture_refinement)
                out << "# Mixtures  : nested repeat refinement (experimental; at most two partition levels)\n";
            if (settings.two_sided_boundaries)
                out << "# Rows      : observed orientations; unordered\n";
            if (experimental_sidecar)
                out << "# Kind      : experimental alternative-repeat views\n";
            out << "\n";
            check_output();
            std::cout << "  ▸ Writing " << (fs::path(output_dir) /
                fs::path(output_filename).filename()).string() << std::endl;
        };

        int array_num = 0;
        auto write_array = [&](const ConsensusArray& arr,
                               std::size_t header_arrays,
                               std::size_t header_spacers,
                               bool experimental_sidecar) {
            ++array_num;
            int lines_needed = 2 + (int)arr.entries.size() + 1;
            if (current_lines + lines_needed > MAX_LINES_PER_FILE && current_lines > 0)
                open_new_file(header_arrays, header_spacers, experimental_sidecar);

            out << ">Array_" << array_num << "  spacers=" << arr.entries.size() << "\n";
            out << arr.consensus << "\n";
            current_lines += 2;

            for (const auto& [repeat_var, spacer] : arr.entries) {
                const std::string& rv = (repeat_var == arr.consensus)
                    ? std::string(arr.consensus.size(), '-')
                    : repeat_var;
                out << "        " << rv << "\t" << spacer << "\n";
                ++current_lines;
            }

            out << "\n";
            ++current_lines;
            check_output();
        };

        // Preserve the release byte structure of all primary files: the
        // primary global counts and line splitting are exactly the same as if
        // supplemental views did not exist.
        open_new_file(primary_array_count, primary_row_count);
        for (std::size_t index = 0;
             index < primary_array_count; ++index) {
            write_array(consensus_arrays[index],
                        primary_array_count, primary_row_count, false);
        }

        // Supplemental views always start in the next numbered file, even if
        // the final primary file has room. The same rule naturally follows
        // any number of primary files and supports supplemental file splits.
        if (supplemental_arrays > 0) {
            close_output();
            open_new_file(supplemental_arrays, supplemental_rows, true);
            for (std::size_t index = primary_array_count;
                 index < consensus_arrays.size(); ++index) {
                write_array(consensus_arrays[index],
                            supplemental_arrays, supplemental_rows, true);
            }
        }

        close_output();

        return file_index - 1;
    }

public:
    explicit PostProcessor(Settings& s) : settings(s) {}

    std::map<std::string, std::vector<std::string>> getSystems() const {
        return crispr_arrays;
    }

    void run_analysis() {
        const auto started = std::chrono::steady_clock::now();
        if (settings.repeat_mixture_refinement)
            std::cout << "  ▸ Experimental repeat-mixture refinement: two partition levels, five/six-base signatures\n";
        consensus_arrays.clear();
        crispr_arrays.clear();
        Zuck::candidate_evidence::Ledger candidate_ledger;
        Zuck::supported_calls::RecruitmentEvidence recruited_evidence;
        std::unique_ptr<Zuck::boundary_candidates::Search> boundary_search;
        std::map<std::string, std::size_t> counts = {
            {"loaded_cycle_records", 0}, {"anchor_groups", 0},
            {"extracted_spacer_observations", 0}, {"exact_repeat_spacer_alleles", 0},
            {"repeat_variants", 0}, {"bounded_repeat_clusters", 0},
            {"spacer_alleles_in_families", 0}, {"spacer_components", 0},
            {"largest_spacer_component", 0}, {"filtered_fewer_than_two_components", 0},
            {"filtered_fewer_than_two_spacers_after_trim", 0},
            {"filtered_repeat_length", 0}, {"filtered_spacer_repeat_ratio", 0},
            {"filtered_short_spacer_singletons", 0},
            {"terminal_variant_rows_rescued", 0}, {"read_linkage_applied", 0},
            {"coalesced_strand_equivalent_arrays", 0},
            {"arrays_before_coalescing", 0}, {"primary_arrays", 0},
            {"primary_entries", 0}, {"supplemental_arrays", 0}, {"supplemental_entries", 0}
        };
        std::string input_path = settings.cycles_input.empty()
            ? settings.output_folder + "/cycles.txt"
            : settings.cycles_input;
        std::string output_dir = settings.output_folder;
        fs::create_directories(output_dir);

        auto cycles = load_rotated_cycles(input_path);
        std::cout << "  ▸ Loaded " << cycles.size() << " cycles" << std::endl;
        counts["loaded_cycle_records"] = cycles.size();
        const std::string input_checksum = result_bundle::file_checksum(input_path);
        auto additional_inputs = settings.effective_input_reads();
        const auto grouping_reads = settings.effective_grouping_reads();
        additional_inputs.insert(additional_inputs.end(), grouping_reads.begin(), grouping_reads.end());
        if (!settings.array_database.empty()) additional_inputs.push_back(settings.array_database);
        result_bundle::Transaction result_transaction(output_dir, input_path,
            additional_inputs);
        const auto result_directory = result_transaction.directory();
        Zuck::array_database::Result database_result;
        const std::string database_checksum = settings.array_database.empty() ? "" :
            result_bundle::file_checksum(settings.array_database);
        std::unique_ptr<Zuck::array_database::Database> database;
        if (!settings.array_database.empty()) {
            database = std::make_unique<Zuck::array_database::Database>(settings.array_database, settings.database_include_possible);
            database_result = database->infer(cycles, settings.effective_grouping_reads(),
                settings.database_read_validation, settings.database_novel_spacers);
        }
        // Apply explicit output constraints before assigning row IDs or writing
        // linkage/database membership tables. No constraints means no change.
        const auto filter_lengths = [&](auto& arrays) {
            const auto& limits = settings.dna_sequence_settings;
            if (!limits.active()) return;
            for (auto& array : arrays) {
                const auto before = array.entries.size();
                array.entries.erase(std::remove_if(array.entries.begin(), array.entries.end(),
                    [&](const auto& pair) { return !limits.pair_allowed(pair.first, pair.second); }),
                    array.entries.end());
                counts["length_filtered_rows"] += before - array.entries.size();
            }
            arrays.erase(std::remove_if(arrays.begin(), arrays.end(), [&](const auto& array) {
                std::set<std::string> spacers;
                for (const auto& pair : array.entries) spacers.insert(pair.second);
                const bool reject = !limits.repeat_allowed(array.consensus.size()) || spacers.size() < 2;
                if (reject) ++counts["length_filtered_groups"];
                return reject;
            }), arrays.end());
        };
        // Database evidence is indexed alongside entries. Filter both together.
        if (settings.dna_sequence_settings.active()) {
            for (auto& group : database_result.groups) {
                for (std::size_t i = group.entries.size(); i-- > 0;) {
                    const auto& pair = group.entries[i];
                    if (settings.dna_sequence_settings.pair_allowed(pair.first, pair.second)) continue;
                    group.entries.erase(group.entries.begin() + i);
                    group.evidence.erase(group.evidence.begin() + i);
                    ++counts["length_filtered_database_rows"];
                }
            }
            filter_lengths(database_result.groups);
        }
        const auto append_database_groups = [&]() {
            if (settings.array_database.empty()) return;
            if (result_bundle::file_checksum(settings.array_database) != database_checksum)
                throw std::runtime_error("Array database changed during inference");
            Zuck::array_database::append_groups(consensus_arrays, database_result,
                (result_directory / "database_evidence.tsv").string(), counts);
            if (settings.database_groups && settings.database_read_validation) {
                const auto grouped = Zuck::database_groups::infer(consensus_arrays,
                    counts.at("database_de_novo_groups"), database->records(), settings.database_include_possible);
                Zuck::database_groups::write_tables(result_directory.string(), grouped);
                counts.insert(grouped.counts.begin(), grouped.counts.end());
            }
            crispr_arrays.clear();
            for (const auto& array : consensus_arrays)
                for (const auto& entry : array.entries)
                    crispr_arrays[array.consensus].push_back(entry.second);
        };
        std::optional<Zuck::read_linked::Result> read_linked_evidence;
        std::vector<std::size_t> read_linked_row_counts;
        const auto write_read_linked_views = [&]() {
            const auto reads = settings.effective_grouping_reads();
            if (!settings.two_sided_boundaries || !settings.read_linked_groups || reads.empty()) return;
            read_linked_row_counts.clear();
            for (const auto& array : consensus_arrays)
                read_linked_row_counts.push_back(array.entries.size());
            read_linked_evidence.emplace(Zuck::read_linked::infer(consensus_arrays, reads));
            Zuck::read_linked::write_tables(result_directory.string(), *read_linked_evidence);
            counts.insert(read_linked_evidence->counts.begin(), read_linked_evidence->counts.end());
            counts["read_linkage_applied"] = 1;
        };
        std::set<std::size_t> exact_boundary_owners;
        const auto resolve_local_boundaries = [&]() {
            if (!settings.local_boundary_graph) return;
            const auto resolved = Zuck::local_boundary_graph::apply(consensus_arrays,
                settings.effective_grouping_reads(), settings.dna_sequence_settings,
                result_directory);
            counts.insert(resolved.counts.begin(), resolved.counts.end());
            for (const auto& decision : resolved.decisions)
                if (decision.selected_geometry)
                    exact_boundary_owners.insert(decision.key.owner);
            std::cout << "  ▸ Local boundary graph: "
                      << counts["local_boundary_graph_replaced_rows"]
                      << " existing rows recut in "
                      << counts["local_boundary_graph_selected_families"] << " families\n";
        };
        const auto resolve_motif_cliffs = [&]() {
            if (!settings.motif_cliff_boundary) return;
            const auto resolved = Zuck::motif_cliff_boundary::apply(consensus_arrays,
                exact_boundary_owners, settings.dna_sequence_settings, result_directory);
            counts.insert(resolved.counts.begin(), resolved.counts.end());
            std::cout << "  ▸ Motif-cliff boundary: "
                      << counts["motif_cliff_boundary_replaced_rows"]
                      << " rows recut in "
                      << counts["motif_cliff_boundary_selected_haplotypes"]
                      << " exact oriented haplotypes\n";
        };
        const auto write_call_evidence = [&]() {
            // One read scan supplies both optional confidence annotations and
            // the selected-call view. Original rows and memberships stay intact.
            if (!settings.confidence_scores && settings.call_profile == "off") return;
            std::vector<std::pair<std::string, std::string>> pairs;
            for (const auto& array : consensus_arrays)
                pairs.insert(pairs.end(), array.entries.begin(), array.entries.end());
            const std::vector<Zuck::array_database::Record> empty_records;
            std::set<std::pair<std::size_t, std::size_t>> evidence_only_rows;
            for (const auto& [row, proof] : recruited_evidence)
                if (proof.added_row) evidence_only_rows.insert(row);
            const auto next_repeat_prefixes = settings.family_repeat_context
                ? Zuck::element_confidence::family_repeat_prefixes(consensus_arrays, evidence_only_rows)
                : std::vector<std::vector<std::string>>{};
            const auto confidence = Zuck::element_confidence::annotate(pairs,
                settings.effective_grouping_reads(), database ? database->records() : empty_records,
                next_repeat_prefixes);
            if (settings.confidence_scores)
                Zuck::element_confidence::write_tables(result_directory.string(), consensus_arrays, confidence);
            counts.insert(confidence.counts.begin(), confidence.counts.end());
            if (settings.call_profile != "off") {
                const auto calls = Zuck::supported_calls::infer(consensus_arrays, confidence, settings.call_profile,
                    recruited_evidence, settings.read_recruitment,
                    settings.singleton_precision_prune, settings.circular_unit_dedup);
                Zuck::supported_calls::write_tables(result_directory.string(), consensus_arrays, calls);
                if (settings.read_first)
                    Zuck::read_first::write_inventory(result_directory, calls, settings, counts);
                counts.insert(calls.counts.begin(), calls.counts.end());
                if (settings.directed_primary_groups) {
                    std::set<Zuck::directed_haplotype_paths::RowKey> selected_original;
                    for (const auto& row : calls.rows)
                        if (row.selected && !row.recruitment.added_row)
                            selected_original.insert({row.array_index, row.entry_index});
                    const auto directed = Zuck::directed_haplotype_paths::infer(
                        consensus_arrays, settings.effective_grouping_reads(), selected_original);
                    Zuck::directed_haplotype_paths::write_tables(result_directory, directed);
                    counts.insert(directed.counts.begin(), directed.counts.end());
                    if (settings.cross_variant_scaffolds) {
                        if (!read_linked_evidence)
                            throw std::logic_error(
                                "Cross-variant scaffolds require native read-link evidence");
                        const auto scaffolded = Zuck::cross_variant_scaffolds::infer(
                            directed, *read_linked_evidence, consensus_arrays,
                            read_linked_row_counts, selected_original);
                        Zuck::cross_variant_scaffolds::write_table(
                            result_directory, scaffolded);
                        counts.insert(scaffolded.counts.begin(), scaffolded.counts.end());
                        std::cout << "  ▸ Cross-variant scaffolds: "
                                  << scaffolded.counts.at("cross_variant_scaffold_groups")
                                  << " complete row-disjoint groups after "
                                  << scaffolded.counts.at(
                                         "cross_variant_scaffold_component_joins")
                                  << " native read-link joins" << std::endl;
                    }
                    std::cout << "  ▸ Directed primary partition: "
                              << directed.counts.at("directed_path_path_groups")
                              << " read paths plus "
                              << directed.counts.at("directed_path_singleton_groups")
                              << " residual singletons" << std::endl;
                }
                std::cout << "  ▸ Supported calls (" << settings.call_profile << "): "
                          << calls.spacers.size() << " canonical spacers in "
                          << calls.groups.size() << " original sequence groups" << std::endl;
            }
            if (database && result_bundle::file_checksum(settings.array_database) != database_checksum)
                throw std::runtime_error("Array database changed during confidence annotation");
        };
        const auto recruit_reads = [&]() {
            if (!settings.read_recruitment) return;
            const auto recruited = Zuck::read_recruitment::infer(consensus_arrays,
                settings.effective_grouping_reads());
            counts.insert(recruited.counts.begin(), recruited.counts.end());
            recruited_evidence = Zuck::read_recruitment_output::apply(consensus_arrays, recruited,
                settings.dna_sequence_settings, result_directory, counts);
            std::cout << "  ▸ Graph-seeded read recruitment: " << counts["read_recruitment_added_rows"]
                      << " new rows; " << recruited_evidence.size() << " rows with explicit read evidence\n";
        };
        const auto write_candidates = [&]() {
            if (settings.candidate_evidence)
                candidate_ledger.write(result_directory, consensus_arrays,
                    settings.dna_sequence_settings, counts);
        };
        const auto append_boundaries = [&]() {
            if (!boundary_search) return;
            boundary_search->append(consensus_arrays, settings.dna_sequence_settings,
                database && settings.database_read_validation ? settings.boundary_policy : "off",
                result_directory / "boundary_evidence.tsv", counts);
            if (settings.candidate_evidence) boundary_search->write(result_directory / "boundary_candidates.tsv");
            crispr_arrays.clear();
            for (const auto& array : consensus_arrays)
                for (const auto& pair : array.entries) crispr_arrays[array.consensus].push_back(pair.second);
        };

        auto spacer_data = extract_spacer_observations(
            std::move(cycles), counts["anchor_groups"]);
        counts["extracted_spacer_observations"] = spacer_data.size();
        if (settings.two_sided_boundaries && (settings.candidate_evidence ||
            (database && settings.database_read_validation && settings.boundary_policy != "off"))) {
            std::set<Zuck::boundary_candidates::Pair> seeds;
            for (const auto& row : spacer_data) seeds.emplace(row.repeat,row.spacer);
            boundary_search = std::make_unique<Zuck::boundary_candidates::Search>(settings.boundary_radius);
            const std::vector<Zuck::array_database::Record> empty_records;
            boundary_search->prepare(seeds, settings.effective_grouping_reads(),
                database ? database->records() : empty_records);
        }
        if (settings.candidate_evidence) {
            for (const auto& row : spacer_data)
                candidate_ledger.observe(row.repeat, row.spacer, row.cycle);
            candidate_ledger.score(settings.effective_grouping_reads());
        }
        if (spacer_data.empty() && database_result.groups.empty()) {
            write_read_linked_views();
            recruit_reads();
            append_database_groups();
            append_boundaries();
            write_call_evidence();
            write_candidates();
            write_structured_results(result_directory, input_path, input_checksum, counts, {},
                std::chrono::duration<double>(std::chrono::steady_clock::now() - started).count());
            result_transaction.commit();
            std::cout << "  No spacers found." << std::endl;
            return;
        }

        // ==========================================================
        // Step 5: Preserve exact spacer evidence by repeat variant.
        //
        // The upstream global shared-23-mer union could transitively join
        // unrelated spacers and then assign the entire component to only one
        // repeat. It also could not deduplicate spacers shorter than 23 nt.
        // Aggregate exact observations here; compatible repeat variants are
        // combined below, without allowing one family to steal another's
        // spacer evidence.
        // ==========================================================
        std::map<std::string, std::map<std::string, std::size_t>>
            repeat_to_spacer_support;
        for (const auto& sd : spacer_data)
            ++repeat_to_spacer_support[sd.repeat][sd.spacer];
        std::vector<SpacerData>().swap(spacer_data);

        std::size_t exact_alleles = 0;
        for (const auto& [repeat, spacers] : repeat_to_spacer_support)
            exact_alleles += spacers.size();
        counts["exact_repeat_spacer_alleles"] = exact_alleles;
        counts["repeat_variants"] = repeat_to_spacer_support.size();

        std::cout << "  ▸ Exact repeat/spacer alleles: " << exact_alleles
                  << " across " << repeat_to_spacer_support.size()
                  << " repeat variants" << std::endl;

        // ==========================================================
        // Step 6: Cluster similar repeats by edit distance,
        //         SPOA consensus on front repeats
        // ==========================================================
        const auto repeat_families = cluster_repeat_families(repeat_to_spacer_support);
        const auto& all_repeats = repeat_families.repeats;
        const auto& repeat_clusters = repeat_families.members;

        std::cout << "  ▸ Repeat groups (post-SPOA merge): " << repeat_clusters.size() << std::endl;
        counts["bounded_repeat_clusters"] = repeat_clusters.size();

        // ==========================================================
        // Step 7: Per consensus group — front SPOA + validated tail SPOA
        // ==========================================================
        consensus_arrays.clear();
        crispr_arrays.clear();

        // Rescue evidence stays separate until every original acceptance and
        // grouping decision is complete, so enabling it cannot replace a row.
        std::map<std::string, std::vector<std::pair<std::string, std::string>>>
            terminal_variants_by_consensus;
        std::size_t total_spacer_alleles = 0;
        std::size_t total_spacer_components = 0;
        std::size_t largest_spacer_component = 0;

        for (const auto& rep_indices : repeat_clusters) {
            std::vector<std::pair<std::string, std::string>> family_terminal_variants;
            std::vector<std::pair<std::string, int>> weighted_repeats;
            std::map<std::string, std::map<std::string, std::size_t>>
                family_spacer_support;

            for (int ri : rep_indices) {
                const std::string& rep = all_repeats[ri];
                const auto& spacers = repeat_to_spacer_support[rep];
                weighted_repeats.emplace_back(rep, (int)spacers.size());
                for (const auto& [spacer, support] : spacers)
                    family_spacer_support[spacer][rep] += support;
            }

            std::vector<std::string> alleles;
            alleles.reserve(family_spacer_support.size());
            for (const auto& [spacer, repeat_support] : family_spacer_support)
                alleles.push_back(spacer);
            total_spacer_alleles += alleles.size();

            // Keep the established compact collapse, but scope it to this
            // bounded repeat family so components cannot steal evidence from
            // unrelated repeats.
            std::vector<std::size_t> allele_support(alleles.size(), 0);
            for (int idx = 0; idx < static_cast<int>(alleles.size()); ++idx)
                for (const auto& [repeat, support] :
                     family_spacer_support[alleles[idx]])
                    allele_support[idx] += support;

            const auto allele_clusters =
                Zuck::post_processing::shared_kmer_components(alleles, GROUP_KMER_SIZE);

            total_spacer_components += allele_clusters.size();
            for (const auto& cluster : allele_clusters)
                largest_spacer_component =
                    std::max(largest_spacer_component, cluster.size());

            // Start from the best-supported observation, then remove at most
            // eight exact terminal overhang bases by selecting the shortest
            // contained member. This corrects conserved repeat sequence that
            // was shifted into a spacer without chasing transitive substrings.
            std::vector<std::pair<std::string, std::string>> merged_entries;
            std::map<std::string, std::map<std::string, std::size_t>>
                emitted_spacer_support;
            for (const auto& cluster : allele_clusters) {
                int representative = cluster.front();
                for (int candidate_idx : cluster) {
                    if (allele_support[candidate_idx] >
                            allele_support[representative] ||
                        (allele_support[candidate_idx] ==
                             allele_support[representative] &&
                         (alleles[candidate_idx].size() >
                              alleles[representative].size() ||
                          (alleles[candidate_idx].size() ==
                               alleles[representative].size() &&
                           alleles[candidate_idx] <
                               alleles[representative]))))
                        representative = candidate_idx;
                }

                if (settings.two_sided_boundaries && settings.terminal_variant_rescue) {
                    for (int candidate : Zuck::co_dominant_terminal_variants(
                             alleles, allele_support, cluster, representative)) {
                        const auto& observed_repeats = family_spacer_support.at(alleles[candidate]);
                        std::string observed_repeat;
                        std::size_t observed_support = 0;
                        for (const auto& [repeat, support] : observed_repeats) {
                            if (support > observed_support ||
                                (support == observed_support &&
                                 (observed_repeat.empty() || repeat < observed_repeat))) {
                                observed_repeat = repeat;
                                observed_support = support;
                            }
                        }
                        family_terminal_variants.emplace_back(observed_repeat, alleles[candidate]);
                    }
                }

                const int support_anchor = representative;
                for (int candidate_idx : cluster) {
                    // Corrected extraction already assigns both boundaries.
                    // Reapplying terminal shortening can choose an error allele.
                    if (settings.two_sided_boundaries) break;
                    if (alleles[candidate_idx].size() >=
                            alleles[representative].size() ||
                        alleles[support_anchor].size() -
                                alleles[candidate_idx].size() >
                            8)
                        continue;
                    const auto& shorter = alleles[candidate_idx];
                    const auto& longer = alleles[support_anchor];
                    std::string shorter_rc = revcomp(shorter);
                    auto is_terminal_allele = [&](const std::string& allele) {
                        if (allele.size() > longer.size()) return false;
                        return longer.compare(0, allele.size(), allele) == 0 ||
                               longer.compare(longer.size() - allele.size(),
                                              allele.size(), allele) == 0;
                    };
                    if (is_terminal_allele(shorter) ||
                        is_terminal_allele(shorter_rc))
                        representative = candidate_idx;
                }

                const std::string& spacer = alleles[representative];
                auto& repeat_support = emitted_spacer_support[spacer];
                for (int member : cluster)
                    for (const auto& [repeat, support] :
                         family_spacer_support[alleles[member]])
                        repeat_support[repeat] += support;

                std::string best_repeat;
                std::size_t best_support = 0;
                const auto& chosen_repeat_support = settings.two_sided_boundaries
                    ? family_spacer_support.at(spacer) : repeat_support;
                for (const auto& [repeat, support] : chosen_repeat_support) {
                    if (support > best_support ||
                        (support == best_support &&
                         (best_repeat.empty() || repeat < best_repeat))) {
                        best_repeat = repeat;
                        best_support = support;
                    }
                }
                merged_entries.emplace_back(best_repeat, spacer);
            }

            std::sort(merged_entries.begin(), merged_entries.end(),
                      [](const auto& lhs, const auto& rhs) {
                          if (lhs.second != rhs.second) return lhs.second < rhs.second;
                          return lhs.first < rhs.first;
                      });

            if (merged_entries.size() < 2) {
                ++counts["filtered_fewer_than_two_components"];
                continue;
            }

            // --- Front consensus via SPOA ---
            std::string front_consensus = spoa_consensus(weighted_repeats);

            // --- Tail detection + validation ---
            std::vector<std::string> all_spacers;
            for (const auto& [rv, sp] : merged_entries)
                all_spacers.push_back(sp);

            int tail_len = settings.two_sided_boundaries ? 0 : detect_tail_length(all_spacers);
            std::string tail_consensus;
            bool tail_valid = false;

            if (tail_len > 0) {
                // Extract tail portions, SPOA consensus
                std::vector<std::pair<std::string, int>> tail_seqs;
                for (const auto& sp : all_spacers) {
                    std::string tail = sp.substr(sp.size() - tail_len);
                    tail_seqs.emplace_back(tail, 1);
                }
                tail_consensus = spoa_consensus(tail_seqs);

                // VALIDATE: tail must match a prefix of front_consensus
                tail_valid = validate_tail_matches_front(tail_consensus, front_consensus);

                if (!tail_valid) {
                    tail_consensus.clear();
                }
            }

            // Trimming can make previously distinct boundary variants exact
            // duplicates. Trim only rows whose own suffix supports the
            // family-level tail and carry their original observation support
            // through the transformation.
            std::map<std::string, std::map<std::string, std::size_t>>
                trimmed_spacer_support;
            for (const auto& [repeat, spacer] : merged_entries) {
                std::string output_spacer = spacer;
                if (tail_valid &&
                    static_cast<int>(spacer.size()) - tail_len >=
                        MIN_SPACER_AFTER_TRIM) {
                    const std::string row_tail =
                        spacer.substr(spacer.size() - tail_len);
                    if (edit_distance(row_tail, tail_consensus,
                                      MAX_TAIL_PREFIX_EDIT_DIST) <=
                        MAX_TAIL_PREFIX_EDIT_DIST)
                        output_spacer =
                            spacer.substr(0, spacer.size() - tail_len);
                }

                std::size_t support = 1;
                auto spacer_it = emitted_spacer_support.find(spacer);
                if (spacer_it != emitted_spacer_support.end()) {
                    auto repeat_it = spacer_it->second.find(repeat);
                    if (repeat_it != spacer_it->second.end())
                        support = repeat_it->second;
                }
                trimmed_spacer_support[output_spacer][repeat] += support;
            }

            merged_entries.clear();
            for (const auto& [spacer, repeat_support] : trimmed_spacer_support) {
                std::string best_repeat;
                std::size_t best_support = 0;
                for (const auto& [repeat, support] : repeat_support) {
                    if (support > best_support ||
                        (support == best_support &&
                         (best_repeat.empty() || repeat < best_repeat))) {
                        best_repeat = repeat;
                        best_support = support;
                    }
                }
                merged_entries.emplace_back(best_repeat, spacer);
            }

            if (merged_entries.size() < 2) {
                ++counts["filtered_fewer_than_two_spacers_after_trim"];
                continue;
            }

            // Full consensus = tail_consensus + front_consensus
            std::string full_consensus = tail_consensus + front_consensus;

            // --- Sanity filter: repeat length ---
            if ((int)full_consensus.size() > MAX_REPEAT_LEN ||
                (int)full_consensus.size() < MIN_REPEAT_LEN) {
                ++counts["filtered_repeat_length"];
                std::cout << "  [filtered] repeat len=" << full_consensus.size()
                          << "bp (allowed " << MIN_REPEAT_LEN << "–" << MAX_REPEAT_LEN << "): "
                          << full_consensus.substr(0, 30) << "\n";
                continue;
            }
            // --- Sanity filter: median spacer / repeat ratio ---
            std::vector<int> spacer_lens;
            for (const auto& [rv, sp] : merged_entries)
                spacer_lens.push_back((int)sp.size());

            double med_spacer = median(spacer_lens);
            double ratio = med_spacer / (double)full_consensus.size();
            if (ratio < MIN_MEDIAN_SPACER_REPEAT_RATIO) {
                ++counts["filtered_spacer_repeat_ratio"];
                std::cout << "  [filtered] spacer/repeat ratio="
                          << std::defaultfloat << std::setprecision(4) << ratio
                          << " (min " << MIN_MEDIAN_SPACER_REPEAT_RATIO << "): "
                          << full_consensus.substr(0, 30) << "\n";
                continue;
            }

            consensus_arrays.push_back({full_consensus, merged_entries});
            auto& accepted_terminal_variants = terminal_variants_by_consensus[full_consensus];
            accepted_terminal_variants.insert(accepted_terminal_variants.end(),
                family_terminal_variants.begin(), family_terminal_variants.end());

            std::vector<std::string> spacers;
            for (const auto& [rv, sp] : merged_entries)
                spacers.push_back(sp);
            crispr_arrays[full_consensus] = spacers;
        }

        std::cout << "  ▸ Spacer components (transitive): "
                  << total_spacer_components << " from "
                  << total_spacer_alleles << " alleles; largest="
                  << largest_spacer_component << std::endl;
        counts["spacer_components"] = total_spacer_components;
        counts["spacer_alleles_in_families"] = total_spacer_alleles;
        counts["largest_spacer_component"] = largest_spacer_component;
        counts["arrays_before_coalescing"] = consensus_arrays.size();
        decltype(repeat_to_spacer_support)().swap(repeat_to_spacer_support);

        // Different bounded raw-repeat families can converge to the same
        // final repeat consensus. Coalesce those blocks deterministically so
        // one biological family is not reported twice.
        std::map<std::string, std::map<std::string, std::string>>
            coalesced_arrays;
        for (const auto& array : consensus_arrays) {
            auto& spacer_to_repeat = coalesced_arrays[array.consensus];
            for (const auto& [repeat, spacer] : array.entries) {
                auto [it, inserted] = spacer_to_repeat.emplace(spacer, repeat);
                if (!inserted && repeat < it->second)
                    it->second = repeat;
            }
        }

        consensus_arrays.clear();
        crispr_arrays.clear();
        for (const auto& [consensus, spacer_to_repeat] : coalesced_arrays) {
            if (spacer_to_repeat.size() < 2) continue;
            ConsensusArray array;
            array.consensus = consensus;
            std::vector<std::string> spacers;
            for (const auto& [spacer, repeat] : spacer_to_repeat) {
                array.entries.emplace_back(repeat, spacer);
                spacers.push_back(spacer);
            }
            if (settings.two_sided_boundaries && settings.short_spacer_filter &&
                Zuck::has_only_short_terminal_variants(array.entries)) {
                ++counts["filtered_short_spacer_singletons"];
                continue;
            }
            consensus_arrays.push_back(std::move(array));
            crispr_arrays[consensus] = std::move(spacers);
        }

        // Record groups that survived all original filters. An absent exact
        // consensus may use its reverse-complement owner only if strand
        // coalescing removed it, not if an earlier filter rejected it.
        std::set<std::string> consensuses_before_strand_coalescing;
        for (const auto& array : consensus_arrays)
            consensuses_before_strand_coalescing.insert(array.consensus);

        // A whole-group strand reversal identifies an equivalent unordered
        // family, but both observed orientations carry literal paired sequence
        // evidence. Coalesce only exact full-row-equivalent groups and preserve
        // every original oriented row and its multiplicity. Legacy is unchanged.
        if (settings.two_sided_boundaries) {
            counts["coalesced_strand_equivalent_arrays"] =
                Zuck::coalesce_strand_equivalent_arrays(consensus_arrays, revcomp);
            crispr_arrays.clear();
            for (const auto& array : consensus_arrays) {
                auto& spacers = crispr_arrays[array.consensus];
                for (const auto& entry : array.entries)
                    spacers.push_back(entry.second);
            }
        }

        // Add tied terminal alleles only to the original accepted family's
        // group. Prefer the exact consensus: opposite-strand groups with the
        // same canonical consensus can have different content. A reverse-
        // complement fallback is valid only after exact strand coalescing has
        // removed the original group's otherwise identical representation.
        if (settings.two_sided_boundaries && settings.terminal_variant_rescue) {
            std::map<std::string, std::size_t> group_for_consensus;
            for (std::size_t i = 0; i < consensus_arrays.size(); ++i)
                group_for_consensus[consensus_arrays[i].consensus] = i;
            for (const auto& [consensus, variants] : terminal_variants_by_consensus) {
                auto group = group_for_consensus.find(consensus);
                if (group == group_for_consensus.end() &&
                    consensuses_before_strand_coalescing.count(consensus))
                    group = group_for_consensus.find(revcomp(consensus));
                if (group == group_for_consensus.end()) continue;
                auto& array = consensus_arrays[group->second];
                std::set<std::pair<std::string, std::string>> observed(
                    array.entries.begin(), array.entries.end());
                for (const auto& pair : variants) {
                    if (!observed.insert(pair).second) continue;
                    array.entries.push_back(pair);
                    crispr_arrays[array.consensus].push_back(pair.second);
                    ++counts["terminal_variant_rows_rescued"];
                }
            }
        }

        // ==========================================================
        // Experimental supplemental alternative-repeat-allele views
        // ==========================================================
        // These are deterministic views over rows already accepted by the
        // frozen primary pipeline. They are deliberately created only after
        // final coalescing, so they cannot influence primary consensus,
        // clustering, filtering, row order, or Array IDs.
        // Linking emits row-index sidecars only, after primary evidence is final.
        // The transaction also covers read parse failures and stale sidecar removal.
        filter_lengths(consensus_arrays);
        // Boundary repair precedes every grouping and selected-call view.
        resolve_local_boundaries();
        resolve_motif_cliffs();
        // Link only the original graph inventory. Evidence-only novel rows must
        // not join, split or suppress components in the selected linked view.
        write_read_linked_views();
        recruit_reads();
        crispr_arrays.clear();
        for (const auto& array : consensus_arrays)
            for (const auto& pair : array.entries) crispr_arrays[array.consensus].push_back(pair.second);
        // Database groups carry their own cycle evidence and do not change
        // any de novo parent rows, IDs, or read-linkage candidate set.
        append_database_groups();
        append_boundaries();
        write_candidates();
        const std::size_t primary_array_count = consensus_arrays.size();
        std::size_t primary_row_count = 0;
        for (const auto& array : consensus_arrays)
            primary_row_count += array.entries.size();

        const auto supplemental_parents =
            append_alternative_views(primary_array_count, primary_row_count, recruited_evidence);
        const std::size_t supplemental_arrays = supplemental_parents.size();
        std::size_t supplemental_rows = 0;
        for (std::size_t index = primary_array_count; index < consensus_arrays.size(); ++index)
            supplemental_rows += consensus_arrays[index].entries.size();

        const int output_files = write_text_results(result_directory, output_dir,
            primary_array_count, primary_row_count, supplemental_arrays, supplemental_rows);

        counts["primary_arrays"] = primary_array_count;
        counts["primary_entries"] = primary_row_count;
        counts["supplemental_arrays"] = supplemental_arrays;
        counts["supplemental_entries"] = supplemental_rows;
        write_call_evidence();
        write_structured_results(result_directory, input_path, input_checksum, counts,
            supplemental_parents,
            std::chrono::duration<double>(std::chrono::steady_clock::now() - started).count());
        result_transaction.commit();

        const std::size_t total_spacers =
            primary_row_count + supplemental_rows;

        std::cout << "  ▸ Arrays: " << consensus_arrays.size()
                  << "   Spacers: " << total_spacers
                  << "   Files: " << output_files << std::endl;
    }
};

#endif // POST_PROCESSING_H
