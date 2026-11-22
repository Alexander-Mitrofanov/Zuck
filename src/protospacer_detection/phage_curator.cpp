#include "protospacer_detection/phage_curator.h"
#include <fstream>
#include <deque>
#include <cstdio> // for std::remove
#include <functional>
#include <set>
#include <algorithm>
#include <limits>
#include <stdexcept>

// Decode sentinel and ambiguous graph bases without indexing outside the alphabet.
namespace {
inline char safe_base_from_code(uint8_t code) {
    static const char map[] = "NACGT"; // index 0..4
    if (code <= 4) return map[code];
    return 'N';
}
}

// constructor
PhageCurator::PhageCurator(SDBG& sdbg) : sdbg(sdbg) {
    const bool validation = RevalidateAllNodesButSingleton();
    if (validation) {
        std::cout << "Graph nodes have successfully been revalidated." << std::endl;
    }
}

PhageCurator::PhageCurator(
    SDBG& sdbg,
    const std::map<uint64_t, std::map<uint64_t, std::vector<std::vector<uint64_t>>>>& grouped_paths,
    const std::unordered_map<uint64_t, std::vector<std::vector<uint64_t>>>& cycles)
    : sdbg(sdbg), grouped_paths(grouped_paths) {

    const bool validation = RevalidateAllNodesButSingleton();
    if (validation) {
        std::cout << "Graph nodes have successfully been revalidated." << std::endl;
    }

    for (const auto& [id, group_cycles] : cycles) {
        for (const auto& path : group_cycles)
            cycle_nodes.insert(path.begin(), path.end());
    }
}

std::string PhageCurator::_ReconstructPath(const std::vector<uint64_t>& path) {
    if (path.empty()) return std::string();
    std::string result = _FetchFirstNode(path.front());
    for (size_t i = 1; i < path.size(); ++i) {
        result += _FetchNodeLastBase(path[i]);
    }
    return result;
}

std::string PhageCurator::_FetchFirstNode(size_t node) {
    if (node >= sdbg.size()) throw std::out_of_range("Path node is outside the graph");
    std::string label;
    const size_t k = sdbg.k();
    if (k == 0) return label;

    std::vector<uint8_t> seq(k, 0);
    sdbg.GetLabel(node, seq.data());

    for (const auto code : seq) label.push_back(safe_base_from_code(code));
    return label;
}

std::string PhageCurator::_FetchNodeLastBase(size_t node) {
    if (node >= sdbg.size()) throw std::out_of_range("Path node is outside the graph");
    const size_t k = sdbg.k();
    if (k == 0) return std::string();

    std::vector<uint8_t> seq(k, 0);
    sdbg.GetLabel(node, seq.data());
    char base = safe_base_from_code(seq.back());
    return std::string(1, base);
}

void PhageCurator::ReconstructPaths(const std::vector<std::vector<uint64_t>>& paths) {
    reconstructed_sequences.clear();
    reconstructed_sequences.reserve(paths.size());
    for (const auto& path : paths) {
        if (path.empty()) continue; // guard
        std::string result_path = _FetchFirstNode(path.front());
        for (size_t i = 1; i < path.size(); ++i) {
            result_path += _FetchNodeLastBase(path[i]);
        }
        reconstructed_sequences.push_back(std::move(result_path));
    }
}

void PhageCurator::WriteSequencesToFasta(const std::string& filename) {
    std::ofstream out(filename);
    if (!out) {
        throw std::runtime_error("Cannot open FASTA output: " + filename);
    }
    int count = 1;
    for (const auto& seq : reconstructed_sequences) {
        out << ">sequence_" << count << "\n" << seq << "\n";
        ++count;
    }
    out.close();
    if (!out) throw std::runtime_error("Cannot write FASTA output: " + filename);
}

bool PhageCurator::RevalidateAllNodesButSingleton() {
    // SetValidEdge updates packed bits. A serial pass avoids races and the
    // overhead of taking a critical-section lock separately for every node.
    for (uint64_t node = 0; node < sdbg.size(); ++node) {
        if (sdbg.EdgeMultiplicity(node) > 1 && !sdbg.IsValidEdge(node))
            sdbg.SetValidEdge(node);
    }
    return true;
}

std::vector<std::vector<uint64_t>> PhageCurator::BeamSearchPathsAvoiding(
    uint64_t start,
    int lower,
    int higher,
    const std::set<uint64_t>& forbidden,
    int beam_width,
    double min_mult,
    double max_mult,
    std::function<void(const std::vector<uint64_t>&)> path_callback) {

    std::vector<std::vector<uint64_t>> all_paths;
    if (!path_callback) {
        all_paths.reserve(static_cast<size_t>(std::max(0, beam_width)));
    }

    if (lower < 0 || higher < lower || beam_width <= 0 || min_mult > max_mult)
        throw std::invalid_argument("Invalid beam search depth, width or multiplicity bounds");
    if (start >= sdbg.size() || !sdbg.IsValidEdge(start)) return all_paths;

    // Expanding a branch appends paths while its siblings still read the
    // parent. deque keeps that parent reference valid across appends.
    std::deque<std::vector<uint64_t>> path_pool;

    auto comp = std::greater<std::pair<double, size_t>>();
    std::set<std::pair<double, size_t>, decltype(comp)> beam_set(comp);

    size_t unique_id = 0;

    // Initial multiplicity check
    const double initial_mult = sdbg.EdgeMultiplicity(start);
    if (initial_mult <= 1.0 || initial_mult < min_mult || initial_mult > max_mult) {
        return all_paths;
    }

    path_pool.emplace_back(std::vector<uint64_t>{start});
    beam_set.insert({initial_mult, unique_id++});

    while (!beam_set.empty()) {
        // Pop best
        const auto it = beam_set.begin();
        const double score = it->first;
        const size_t id = it->second;
        beam_set.erase(it);

        const std::vector<uint64_t>& path_ref = path_pool[id];
        const uint64_t v = path_ref.back();
        const int current_depth = static_cast<int>(path_ref.size()) - 1;

        if (current_depth >= lower && current_depth <= higher) {
            if (path_callback) {
                path_callback(path_ref);
            } else {
                all_paths.push_back(path_ref);
            }
            // Preserve original behavior: do not expand after collection
            continue;
        }

        if (current_depth >= higher || sdbg.EdgeOutdegreeZero(v)) {
            continue;
        }

        const int outdegree = sdbg.EdgeOutdegree(v);
        if (outdegree <= 0 || !sdbg.IsValidEdge(v)) {
            continue;
        }

        // SAFETY: cap iteration by buffer capacity
        const int MAX_EDGE_COUNT = 4;
        uint64_t neighbors[MAX_EDGE_COUNT] = {0, 0, 0, 0};
        const int flag = sdbg.OutgoingEdges(v, neighbors);
        if (flag == -1) {
            continue;
        }
        const int n_iter = std::min(outdegree, MAX_EDGE_COUNT);

        for (int i = 0; i < n_iter; ++i) {
            const uint64_t neighbor = neighbors[i];
            if (!sdbg.IsValidEdge(neighbor)) continue;

            // Cycle check (on current path)
            if (std::find(path_ref.begin(), path_ref.end(), neighbor) != path_ref.end()) {
                continue;
            }

            // Forbidden check
            if (forbidden.find(neighbor) != forbidden.end() && neighbor != start) {
                continue;
            }

            // Multiplicity guard
            const double mult = sdbg.EdgeMultiplicity(neighbor);
            if (mult <= 1.0 || mult < min_mult || mult > max_mult) {
                continue;
            }

            // Each branch owns its extension of the stable parent path.
            std::vector<uint64_t> new_path = path_ref;
            new_path.push_back(neighbor);

            double new_score = (score * path_ref.size() + mult) /
                               static_cast<double>(path_ref.size() + 1);
            uint64_t current = neighbor;

            // Extend along simple path (preserving original loop’s semantics)
            while (true) {
                const bool depth_limit = (static_cast<int>(new_path.size()) - 1) >= higher;
                if (depth_limit) break;

                const uint64_t next = sdbg.NextSimplePathEdge(current);
                if (next == SDBG::kNullID || !sdbg.IsValidEdge(next)) break;

                // Cycle check
                if (std::find(new_path.begin(), new_path.end(), next) != new_path.end()) {
                    break;
                }
                // Forbidden check
                if (forbidden.find(next) != forbidden.end() && next != start) {
                    break;
                }
                // Multiplicity check
                const double next_mult = sdbg.EdgeMultiplicity(next);
                if (next_mult <= 1.0 || next_mult < min_mult || next_mult > max_mult) {
                    break;
                }

                new_path.push_back(next);
                const int new_depth = static_cast<int>(new_path.size()) - 1;
                new_score = (new_score * new_depth + next_mult) / static_cast<double>(new_depth + 1);
                current = next;
            }

            // Store the candidate and retain only the best pending branches.
            const size_t new_id = unique_id++;
            path_pool.emplace_back(std::move(new_path));

            beam_set.insert({new_score, new_id});

            // Prune the lowest score if over beam width
            if (beam_width > 0 && beam_set.size() > static_cast<size_t>(beam_width)) {
                auto prune_it = beam_set.end();
                --prune_it;
                beam_set.erase(prune_it);
            }
        }
    }

    return all_paths;
}

// New function using beam search, mirroring FindQualityPathsDLSFromGroupedPaths
std::map<std::string, std::vector<std::string>>
PhageCurator::FindQualityPathsBeamSearchFromGroupedPaths(int min_length, int max_length, const std::string& filename, int beam_width) {
    std::ofstream out(filename, std::ios::app);  // Append mode to add to file
    if (!out) {
        throw std::runtime_error("Cannot open FASTA output: " + filename);
    }

    std::map<std::string, std::vector<std::string>> consensus_map;
    int path_count = 0;

    struct counts {
        int first;
        int second;
        int third;
    };
    counts counts_by_nodes = {0, 0, 0};

    for (const auto& [group_id, cycle_map] : grouped_paths) {
        std::vector<std::string> quality_paths;
        counts_by_nodes.first += 1;
        std::cout << "Processed groups: " << counts_by_nodes.first << "\r";

        for (const auto& [cycle_id, paths_vec] : cycle_map) {
            (void)cycle_id; // unused in current logic
            counts_by_nodes.second += 1;

            for (const auto& path : paths_vec) {
                counts_by_nodes.third += 1;

                if (path.empty()) continue; // guard

                // Start from last node as in original code
                const uint64_t start = path.back();

                // min/max multiplicity window around start
                double min_mult = 0.1 * sdbg.EdgeMultiplicity(start);
                if (min_mult < 1.0) min_mult = 1.0;
                double max_mult = 5.0 * sdbg.EdgeMultiplicity(start);
                if (max_mult < min_mult) max_mult = min_mult * 50.0;

                auto extended_paths = BeamSearchPathsAvoiding(
                    start, min_length, max_length, cycle_nodes, beam_width, min_mult, max_mult, nullptr);

                if (extended_paths.empty()) continue;

                std::vector<std::vector<uint64_t>> best_paths =
                    GetTopPathsFromBeamPaths(extended_paths, min_mult, max_mult, 2);

                if (best_paths.empty()) continue;

                for (const auto& ext_path : best_paths) {
                    if (ext_path.empty()) continue;

                    std::string result_path = _FetchFirstNode(ext_path.front());
                    for (size_t i = 1; i < ext_path.size(); ++i) {
                        result_path += _FetchNodeLastBase(ext_path[i]);
                    }

                    out << ">quality_path_" << path_count << "\n" << result_path << "\n";
                    ++path_count;
                    quality_paths.push_back(result_path);
                }
            }
        }

        std::cout << "Processed counts:" << counts_by_nodes.first << ","
                  << counts_by_nodes.second << "," << counts_by_nodes.third << ";";

        std::string group_id_str = _FetchFirstNode(group_id);
        consensus_map[group_id_str] = quality_paths;
    }

    std::cout << "\nSaved in " << filename << std::endl;
    out.close();
    if (!out) throw std::runtime_error("Cannot write FASTA output: " + filename);
    return consensus_map;
}

std::string PhageCurator::ComputeConsensusForCurrentGroup(std::vector<std::string> sequences) {
    if (sequences.empty()) return std::string();

    spoa::Graph graph{};
    auto alignment_engine = spoa::AlignmentEngine::Create(spoa::AlignmentType::kNW, 3, -5, -3);
    for (const auto& it : sequences) {
        auto alignment = alignment_engine->Align(it, graph);
        graph.AddAlignment(alignment, it);
    }
    auto consensus = graph.GenerateConsensus();
    return consensus;
}

std::vector<std::vector<uint64_t>> PhageCurator::GetTopPathsFromBeamPaths(
    const std::vector<std::vector<uint64_t>>& beam_paths,
    double min_mult,
    double max_mult,
    size_t top_n) {

    if (min_mult > max_mult) throw std::invalid_argument("Invalid multiplicity bounds");

    struct PathScore {
        std::vector<uint64_t> path;
        double delta;
    };

    std::vector<PathScore> scored_paths;
    scored_paths.reserve(beam_paths.size());

    for (const auto& path : beam_paths) {
        if (path.empty()) continue;

        double total = 0.0;
        for (uint64_t node : path) {
            if (node >= sdbg.size()) throw std::out_of_range("Path node is outside the graph");
            total += sdbg.EdgeMultiplicity(node);
        }
        const double avg = total / static_cast<double>(path.size());

        double delta = 0.0;
        if (avg < min_mult)
            delta = min_mult - avg;
        else if (avg > max_mult)
            delta = avg - max_mult;
        else
            delta = 0.0;

        scored_paths.push_back({path, delta});
    }

    std::sort(scored_paths.begin(), scored_paths.end(),
              [](const PathScore& a, const PathScore& b) {
                  return a.delta != b.delta ? a.delta < b.delta : a.path < b.path;
              });

    const size_t keep = std::min(top_n, scored_paths.size());
    std::vector<std::vector<uint64_t>> top_paths;
    top_paths.reserve(keep);
    for (size_t i = 0; i < keep; ++i) {
        top_paths.push_back(scored_paths[i].path);
    }
    return top_paths;
}
