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
#include "spoa/spoa.hpp"
#include "core/settings.h"
#include "crispr_array/alternative_repeat_view.h"

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
    };

    std::string glue_kmers(const std::string& line) {
        std::istringstream iss(line);
        std::vector<std::string> kmers;
        std::string kmer;
        while (iss >> kmer) kmers.push_back(kmer);
        if (kmers.empty()) return "";

        std::string result = kmers[0];
        for (size_t i = 1; i < kmers.size(); ++i) {
            const std::string& next = kmers[i];
            bool found = false;
            for (int k = (int)result.size(); k > 0; --k) {
                if ((int)next.size() >= k &&
                    result.substr(result.size() - k) == next.substr(0, k)) {
                    result += next.substr(k);
                    found = true;
                    break;
                }
            }
            if (!found) result += next;
        }
        return result;
    }

    int edit_distance(const std::string& a, const std::string& b, int max_dist) {
        int m = (int)a.size(), n = (int)b.size();
        if (std::abs(m - n) > max_dist) return max_dist + 1;

        std::vector<int> prev(n + 1), curr(n + 1);
        for (int j = 0; j <= n; ++j) prev[j] = j;

        for (int i = 1; i <= m; ++i) {
            curr[0] = i;
            int row_min = curr[0];
            for (int j = 1; j <= n; ++j) {
                int cost = (a[i-1] == b[j-1]) ? 0 : 1;
                curr[j] = std::min({prev[j] + 1, curr[j-1] + 1, prev[j-1] + cost});
                row_min = std::min(row_min, curr[j]);
            }
            if (row_min > max_dist) return max_dist + 1;
            std::swap(prev, curr);
        }
        return prev[n];
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

    std::vector<SpacerData> extract_two_sided_group(
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
                core});
        }
        return extracted;
    }

public:
    explicit PostProcessor(Settings& s) : settings(s) {}

    std::map<std::string, std::vector<std::string>> getSystems() const {
        return crispr_arrays;
    }

    void run_analysis() {
        std::string input_path = settings.cycles_input.empty()
            ? settings.output_folder + "/cycles.txt"
            : settings.cycles_input;
        std::string output_dir = settings.output_folder;
        fs::create_directories(output_dir);

        // A replay can legitimately emit fewer files (or no arrays). Remove
        // only this tool's numbered result files after the input has been
        // parsed, so stale blocks can never be mixed into a later score.
        auto remove_numbered_outputs = [&]() {
            static const std::string prefix = "CRISPR_Arrays_";
            for (const auto& entry : fs::directory_iterator(output_dir)) {
                if (!entry.is_regular_file()) continue;
                const std::string name = entry.path().filename().string();
                if (name.size() <= prefix.size() + 4 ||
                    name.compare(0, prefix.size(), prefix) != 0 ||
                    name.compare(name.size() - 4, 4, ".txt") != 0)
                    continue;
                const std::string number =
                    name.substr(prefix.size(), name.size() - prefix.size() - 4);
                if (!number.empty() &&
                    std::all_of(number.begin(), number.end(),
                                [](unsigned char c) { return std::isdigit(c); }))
                    fs::remove(entry.path());
            }
        };

        // ==========================================================
        // ORIGINAL Step 1: Read all cycles and glue k-mers
        // ==========================================================
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
                            std::string glued = glue_kmers(line);
                            const bool valid_dna = std::all_of(
                                glued.begin(), glued.end(), [](unsigned char base) {
                                    base = static_cast<unsigned char>(std::toupper(base));
                                    return base == 'A' || base == 'C' ||
                                           base == 'G' || base == 'T' ||
                                           base == 'N';
                                });
                            if (!valid_dna)
                                throw std::runtime_error(
                                    "Invalid compact/legacy cycle record at " +
                                    input_path + ":" +
                                    std::to_string(line_number));
                            std::transform(glued.begin(), glued.end(), glued.begin(),
                                           [](unsigned char base) {
                                               return static_cast<char>(std::toupper(base));
                                           });
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
                            if (!glued.empty()) raw_cycles.push_back(glued);
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

            cycles.reserve(raw_cycles.size());
            for (const auto& seq : raw_cycles)
                cycles.push_back(settings.two_sided_boundaries
                    ? circular_freq_rotation(seq, kmer_freq)
                    : freq_rotation(seq, kmer_freq));
        }
        std::cout << "  ▸ Loaded " << cycles.size() << " cycles" << std::endl;

        // ==========================================================
        // ORIGINAL Step 2: Group by first 23-mer
        // ==========================================================
        std::unordered_map<std::string, std::vector<std::string>> groups;
        for (const auto& seq : cycles) {
            if ((int)seq.size() >= GROUP_KMER_SIZE)
                groups[seq.substr(0, GROUP_KMER_SIZE)].push_back(seq);
        }
        std::cout << "  ▸ Grouped into " << groups.size() << " repeat groups" << std::endl;

        // ==========================================================
        // Step 3 & 4: Unanimous extend, extract spacers
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
        remove_numbered_outputs();
        if (spacer_data.empty()) {
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

        std::size_t exact_alleles = 0;
        for (const auto& [repeat, spacers] : repeat_to_spacer_support)
            exact_alleles += spacers.size();

        std::cout << "  ▸ Exact repeat/spacer alleles: " << exact_alleles
                  << " across " << repeat_to_spacer_support.size()
                  << " repeat variants" << std::endl;

        // ==========================================================
        // Step 6: Cluster similar repeats by edit distance,
        //         SPOA consensus on front repeats
        // ==========================================================
        std::vector<std::string> all_repeats;
        for (const auto& [rep, _] : repeat_to_spacer_support)
            all_repeats.push_back(rep);

        int nr = (int)all_repeats.size();

        std::vector<std::size_t> repeat_raw_support(nr, 0);
        std::vector<std::size_t> repeat_distinct_spacers(nr, 0);
        for (int i = 0; i < nr; ++i) {
            const auto& spacers = repeat_to_spacer_support[all_repeats[i]];
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

        std::cout << "  ▸ Repeat groups (post-SPOA merge): " << repeat_clusters.size() << std::endl;

        // ==========================================================
        // Step 7: Per consensus group — front SPOA + validated tail SPOA
        // ==========================================================
        consensus_arrays.clear();
        crispr_arrays.clear();

        std::size_t total_spacer_alleles = 0;
        std::size_t total_spacer_components = 0;
        std::size_t largest_spacer_component = 0;

        for (const auto& rep_indices : repeat_clusters) {
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
            std::unordered_map<std::string, std::vector<int>> spacer_kmer_index;
            for (int idx = 0; idx < static_cast<int>(alleles.size()); ++idx) {
                std::unordered_set<std::string> seen_kmers;
                for (int pos = 0;
                     pos + GROUP_KMER_SIZE <= static_cast<int>(alleles[idx].size());
                     ++pos) {
                    std::string kmer = alleles[idx].substr(pos, GROUP_KMER_SIZE);
                    if (seen_kmers.insert(kmer).second)
                        spacer_kmer_index[kmer].push_back(idx);
                }
            }

            std::vector<std::size_t> allele_support(alleles.size(), 0);
            for (int idx = 0; idx < static_cast<int>(alleles.size()); ++idx)
                for (const auto& [repeat, support] :
                     family_spacer_support[alleles[idx]])
                    allele_support[idx] += support;

            // Frozen v1.0.0 collapse: any shared 23-mer joins the full
            // transitive component.
            std::vector<int> allele_parent(alleles.size());
            std::iota(allele_parent.begin(), allele_parent.end(), 0);
            std::function<int(int)> allele_find = [&](int x) -> int {
                if (allele_parent[x] != x)
                    allele_parent[x] = allele_find(allele_parent[x]);
                return allele_parent[x];
            };
            auto allele_unite = [&](int x, int y) {
                int px = allele_find(x), py = allele_find(y);
                if (px != py) allele_parent[px] = py;
            };
            for (const auto& [kmer, indices] : spacer_kmer_index)
                for (std::size_t i = 1; i < indices.size(); ++i)
                    allele_unite(indices.front(), indices[i]);

            std::map<int, std::vector<int>> transitive_clusters;
            for (int idx = 0; idx < static_cast<int>(alleles.size()); ++idx)
                transitive_clusters[allele_find(idx)].push_back(idx);

            std::vector<std::vector<int>> allele_clusters;
            allele_clusters.reserve(transitive_clusters.size());
            for (auto& [root, members] : transitive_clusters)
                allele_clusters.push_back(std::move(members));
            std::sort(allele_clusters.begin(), allele_clusters.end(),
                      [&](const auto& lhs, const auto& rhs) {
                          return alleles[lhs.front()] < alleles[rhs.front()];
                      });

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

            if (merged_entries.size() < 2) continue;

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

            if (merged_entries.size() < 2) continue;

            // Full consensus = tail_consensus + front_consensus
            std::string full_consensus = tail_consensus + front_consensus;

            // --- Sanity filter: repeat length ---
            if ((int)full_consensus.size() > MAX_REPEAT_LEN ||
                (int)full_consensus.size() < MIN_REPEAT_LEN) {
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
                std::cout << "  [filtered] spacer/repeat ratio="
                          << std::defaultfloat << std::setprecision(4) << ratio
                          << " (min " << MIN_MEDIAN_SPACER_REPEAT_RATIO << "): "
                          << full_consensus.substr(0, 30) << "\n";
                continue;
            }

            consensus_arrays.push_back({full_consensus, merged_entries});

            std::vector<std::string> spacers;
            for (const auto& [rv, sp] : merged_entries)
                spacers.push_back(sp);
            crispr_arrays[full_consensus] = spacers;
        }

        std::cout << "  ▸ Spacer components (transitive): "
                  << total_spacer_components << " from "
                  << total_spacer_alleles << " alleles; largest="
                  << largest_spacer_component << std::endl;

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
            consensus_arrays.push_back(std::move(array));
            crispr_arrays[consensus] = std::move(spacers);
        }

        // ==========================================================
        // Experimental supplemental alternative-repeat-allele views
        // ==========================================================
        // These are deterministic views over rows already accepted by the
        // frozen primary pipeline. They are deliberately created only after
        // final coalescing, so they cannot influence primary consensus,
        // clustering, filtering, row order, or Array IDs.
        const std::size_t primary_array_count = consensus_arrays.size();
        std::size_t primary_row_count = 0;
        for (const auto& array : consensus_arrays)
            primary_row_count += array.entries.size();

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
            for (const auto& [repeat, spacer] : array.entries) {
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
            for (const auto& [observed_repeat, observed_spacer] :
                 parent_array.entries) {
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
            primary_row_count / 400;          // <= 0.25%
        std::size_t supplemental_rows = 0;
        std::size_t supplemental_arrays = 0;
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
            output_filename = output_dir + "/CRISPR_Arrays_" +
                              std::to_string(file_index) + ".txt";
            out.open(output_filename);
            if (!out.is_open() || !out)
                throw std::runtime_error("Cannot open arrays file for writing: " + output_filename);
            ++file_index;
            current_lines = 0;
            // Write file header
            out << "# my_mcaat — CRISPR Array Output\n";
            out << "# Generated : " << timestamp_str << "\n";
            out << "# Arrays    : " << header_arrays << "\n";
            out << "# Spacers   : " << header_spacers << "\n";
            if (settings.two_sided_boundaries)
                out << "# Boundary  : two-sided (experimental)\n";
            if (experimental_sidecar)
                out << "# Kind      : experimental alternative-repeat views\n";
            out << "\n";
            check_output();
            std::cout << "  ▸ Writing " << output_filename << std::endl;
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

        const std::size_t total_spacers =
            primary_row_count + supplemental_rows;

        std::cout << "  ▸ Arrays: " << consensus_arrays.size()
                  << "   Spacers: " << total_spacers
                  << "   Files: " << (file_index - 1) << std::endl;
    }
};

#endif // POST_PROCESSING_H
