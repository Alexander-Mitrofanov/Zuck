#include "crispr_array/cycle_finder.h"
#include "crispr_array/bounded_cycle_search.h"
#include "crispr_array/search_partitions.h"
#include "core/settings.h"
#include <algorithm>
#include <array>
#include <iostream>
#include <stdexcept>
#include <fstream>
#include <filesystem>
#include <chrono>
#include <stack>
#include <cmath>
#include <malloc.h>
#include <omp.h>

// Parallel hashmap for better performance in DLS
#include <parallel_hashmap/phmap.h>

using namespace std;

namespace {
void apply_adaptive_coverage(Settings& settings) {
    // Real SDBG edges describe k+1 bases. Exclude sentinel edges at both ends.
    const auto& graph = *settings.sdbg;
    std::array<uint64_t, 256> histogram{};
    uint64_t real_edges = 0;
    for (uint64_t edge = 0; edge < graph.size(); ++edge) {
        if (!graph.IsValidEdge(edge) || graph.IsTip(edge)) continue;
        const unsigned multiplicity = graph.EdgeMultiplicity(edge);
        if (multiplicity == 0) continue;
        ++histogram[std::min(multiplicity, 255u)];
        ++real_edges;
    }

    auto& support = settings.coverage_support;
    support.measured = true;
    support.real_edges = real_edges;
    support.repeated_edges = real_edges - histogram[1];
    const auto quartile = [&](unsigned numerator) -> unsigned {
        const uint64_t n = support.repeated_edges;
        if (n == 0) return 0;
        // Nearest rank ceil(n * numerator / 4), without multiplying a large n.
        const uint64_t rank = (n / 4) * numerator + ((n % 4) * numerator + 3) / 4;
        uint64_t cumulative = 0;
        for (unsigned m = 2; m < histogram.size(); ++m) {
            cumulative += histogram[m];
            if (cumulative >= rank) return m;
        }
        throw std::logic_error("Repeated edge histogram has inconsistent total");
    };
    support.q25 = quartile(1);
    support.q50 = quartile(2);
    support.q75 = quartile(3);
    support.low_support = support.repeated_edges >= 1000 && support.q75 <= 5;
    auto& search = settings.cycle_finder_settings;
    if (support.low_support) {
        if (!search.threshold_explicit)
            search.threshold_multiplicity = std::max(3u, 2 * support.q50);
        if (!search.low_abundance_explicit)
            search.low_abundance = true;
    }
    std::cout << "  ▸ Adaptive coverage: real edges=" << real_edges
              << ", singletons=" << histogram[1] << ", repeated=" << support.repeated_edges
              << ", repeated-edge Q25/Q50/Q75=" << support.q25 << '/' << support.q50 << '/' << support.q75
              << ", low support=" << (support.low_support ? "yes" : "no")
              << ", seed cutoff=" << search.threshold_multiplicity
              << ", retain singletons=" << (search.low_abundance ? "yes" : "no") << '\n';
}
} // namespace

bool CycleFinder::_HasSelfLoop(uint64_t node) {
    uint64_t incoming[MAX_EDGE_COUNT];
    const int count = settings.sdbg->IncomingEdges(node, incoming);
    return count > 0 && std::find(incoming, incoming + count, node) != incoming + count;
}

// Visited history is local to this search partition. Reject zero-multiplicity
// edges before division: imported graphs can still contain such edges.
bool CycleFinder::_BackgroundCheck(uint64_t original_node,
                                   size_t repeat_multiplicity, uint64_t neighbor) {
    if (original_node == neighbor || !settings.sdbg->IsValidEdge(neighbor))
        return false;
    const auto multiplicity = settings.sdbg->EdgeMultiplicity(neighbor);
    if (multiplicity == 0 || repeat_multiplicity / multiplicity > 500)
        return false;
    const size_t word = neighbor >> 6;
    const uint64_t mask = uint64_t{1} << (neighbor & 63);
    return (per_thread_visited[omp_get_thread_num()][word] & mask) == 0;
}

void CycleFinder::_GetOutgoings(uint64_t node, SmallNodeSet& neighbors,
                                size_t repeat_multiplicity) {
    if (!settings.sdbg->IsValidEdge(node)) return;
    uint64_t outgoing[MAX_EDGE_COUNT];
    const int count = settings.sdbg->OutgoingEdges(node, outgoing);
    for (int i = 0; i < count; ++i)
        if (_BackgroundCheck(node, repeat_multiplicity, outgoing[i]))
            neighbors.insert(outgoing[i]);
}

void CycleFinder::_GetIncomings(uint64_t node, SmallNodeSet& neighbors,
                                size_t repeat_multiplicity) {
    if (!settings.sdbg->IsValidEdge(node)) return;
    uint64_t incoming[MAX_EDGE_COUNT];
    const int count = settings.sdbg->IncomingEdges(node, incoming);
    for (int i = 0; i < count; ++i)
        if (_BackgroundCheck(node, repeat_multiplicity, incoming[i]))
            neighbors.insert(incoming[i]);
}

// Reachability checks deliberately ignore multiplicity and visited history.
void CycleFinder::_GetOutgoings(uint64_t node, SmallNodeSet& neighbors) {
    if (!settings.sdbg->IsValidEdge(node)) return;
    uint64_t outgoing[MAX_EDGE_COUNT];
    const int count = settings.sdbg->OutgoingEdges(node, outgoing);
    for (int i = 0; i < count; ++i)
        if (settings.sdbg->IsValidEdge(outgoing[i])) neighbors.insert(outgoing[i]);
}

void CycleFinder::_GetIncomings(uint64_t node, SmallNodeSet& neighbors) {
    if (!settings.sdbg->IsValidEdge(node)) return;
    uint64_t incoming[MAX_EDGE_COUNT];
    const int count = settings.sdbg->IncomingEdges(node, incoming);
    for (int i = 0; i < count; ++i)
        if (settings.sdbg->IsValidEdge(incoming[i])) neighbors.insert(incoming[i]);
}

/**
 * @brief Finds the cycles in the graph using settings configuration.
 * 
 * @param settings The application Settings containing SDBG pointer and cycle finder configuration.
 */
CycleFinder::CycleFinder(Settings& settings)
    : settings(settings) {
    // Ensure settings.sdbg is valid; this should be set by caller before constructing CycleFinder
    if (this->settings.sdbg == nullptr) {
        throw std::runtime_error("CycleFinder requires settings.sdbg to be set to a valid SDBG instance");
    }
    if (settings.adaptive_coverage) apply_adaptive_coverage(settings);
    this->FindApproximateCRISPRArrays();
}

vector<vector<uint64_t>> CycleFinder::FindCycle(uint64_t start_node, vector<uint64_t>& path, phmap::flat_hash_map<uint64_t,int>& lock, vector<SmallNodeSet>& stack, 
                                        vector<int>& backtrack_lengths) {
    int counter = 0;
    uint64_t current_node = start_node;
    vector<vector<uint64_t>> cycles;
    int steps_counter = 0;
    bool saturated = false;

    while (!stack.empty() && !saturated) {
        steps_counter += 1;
        if (steps_counter > 10000000) {
            step_limited_start_nodes.fetch_add(1, std::memory_order_relaxed);
            break;
        }
        
        SmallNodeSet neighbors = stack.back();
        bool flag = true;
        for (auto neighbor : neighbors) {
            current_node = neighbor;
            if (current_node == start_node ) {
                // Consume closure once, just like an explored outgoing edge.
                // Otherwise returning from another branch records this same
                // path again and spends the cycle cap on duplicate evidence.
                stack.back().erase(neighbor);
                backtrack_lengths.back() = 1;
                if (path.size() > static_cast<size_t>(this->settings.cycle_finder_settings.cycle_min_length)) {
                    cycles.push_back(path);
                    counter += 1;
                    if (counter >= this->cluster_bounds) {
                        // A rich start node is still evidence.  The original
                        // implementation erased all 500 retained cycles here,
                        // turning saturation into a complete false negative.
                        saturated = true;
                        saturated_start_nodes.fetch_add(1, std::memory_order_relaxed);
                        flag=false;
                        break;
                    }
                }
            } 
            else if (static_cast<int>(path.size()) < lock.try_emplace(neighbor, this->settings.cycle_finder_settings.cycle_max_length).first->second) {
                neighbors.erase(neighbor);
                path.push_back(neighbor);
                backtrack_lengths.push_back(this->settings.cycle_finder_settings.cycle_max_length);
                lock[neighbor] = path.size();
                stack.back().erase(neighbor);
                SmallNodeSet outgoings;
                this->_GetOutgoings(neighbor, outgoings, this->settings.sdbg->EdgeMultiplicity(start_node));
                stack.push_back(outgoings);
                flag = false;
                break;
            }
        }
        if (flag) {
            stack.pop_back();
            uint64_t v = path.back();
            path.pop_back();
            int backtrack_length = backtrack_lengths.back();
            backtrack_lengths.pop_back();

            if (!backtrack_lengths.empty()) {
                backtrack_lengths.back() = min(backtrack_lengths.back(), backtrack_length);
            }
            if (backtrack_length < this->settings.cycle_finder_settings.cycle_max_length) {
                vector<pair<int, uint64_t>> relax_stack;
                relax_stack.push_back(make_pair(backtrack_length, v));

                phmap::flat_hash_set<uint64_t> path_set(path.begin(), path.end());

                while (!relax_stack.empty()) {
                    int bl = relax_stack.back().first;
                    uint64_t u = relax_stack.back().second;
                    relax_stack.pop_back();
                    if (lock.try_emplace(u, this->settings.cycle_finder_settings.cycle_max_length).first->second < this->settings.cycle_finder_settings.cycle_max_length - bl + 1) {
                        lock[u] = this->settings.cycle_finder_settings.cycle_max_length - bl + 1;
                        SmallNodeSet incomings;
                        this->_GetIncomings(u, incomings, this->settings.sdbg->EdgeMultiplicity(start_node));
                        for (auto w : incomings)
                            if (path_set.find(w) == path_set.end())
                                relax_stack.push_back(make_pair(bl + 1, w));
                    }
                }
            }
        }
    }
    
    
    if (cycles.empty()) return {};

    // Mark visited nodes within this logical search partition only.
    for (const auto& cycle : cycles)
        for (const auto& node : cycle) {
            size_t idx = node >> 6;
            uint64_t mask = 1ULL << (node & 63);
            int tid = omp_get_thread_num();
            per_thread_visited[tid][idx] |= mask;
        }

    return cycles;
}


/**
 * @brief Utility function to initialize and start the cycle finding process from a given node.
 */
vector<vector<uint64_t>> CycleFinder::FindCycleUtil(uint64_t start_node) {
    vector<uint64_t> path;
    phmap::flat_hash_map<uint64_t, int> lock;
    vector<SmallNodeSet> stack;
    vector<int> backtrack_lengths;
    path.push_back(start_node);
    lock[start_node] = 0;
    SmallNodeSet outgoings;
    this->_GetOutgoings(start_node, outgoings, this->settings.sdbg->EdgeMultiplicity(start_node));
    stack.push_back(outgoings);
    backtrack_lengths.push_back(this->settings.cycle_finder_settings.cycle_max_length);
    return FindCycle(start_node, path, lock, stack, backtrack_lengths);
}
/**
 * @brief Performs a depth-limited search to determine if a path exists between two nodes within a certain depth.
 * Breadth-first visitation preserves shorter routes and includes closure
 * at the configured maximum cycle length.
 */
bool CycleFinder::DepthLevelSearch(uint64_t start, uint64_t target, int limit, int& reached_depth) {
    static thread_local Zuck::BoundedCycleWorkspace<
        phmap::flat_hash_set<uint64_t>> workspace;
    auto neighbors = [this](uint64_t node) {
        SmallNodeSet outgoings;
        this->_GetOutgoings(node, outgoings);
        return outgoings;
    };
    return Zuck::HasBoundedPath(start, target, limit, neighbors,
                                   workspace, reached_depth);
}


vector<uint64_t> CycleFinder::CollectTips() {
    size_t n = this->settings.sdbg->size();
    int threads = static_cast<int>(this->settings.threads);
    vector<vector<uint64_t>> local_tips(threads);

    #pragma omp parallel num_threads(threads)
    {
        int tid = omp_get_thread_num();
        #pragma omp for schedule(static)
        for (uint64_t node = 0; node < n; node++) {
            if (this->settings.sdbg->EdgeOutdegree(node) == 0 && this->settings.sdbg->IsValidEdge(node)) {
                local_tips[tid].push_back(node);
            }
        }
    }

    // Merge results while keeping memory usage minimal
    size_t total = 0;
    for (auto &v : local_tips) total += v.size();
    vector<uint64_t> tips;
    tips.reserve(total);
    for (auto &v : local_tips) {
        tips.insert(tips.end(), v.begin(), v.end());
    }
    return tips;
} 

void CycleFinder::RecursiveReduction(uint64_t tip) {
    std::stack<uint64_t> work;
    work.push(tip);
    while (!work.empty()) {
        uint64_t node = work.top();
        work.pop();
        if (this->settings.sdbg->EdgeOutdegree(node) > 0)
            continue;
        if (!this->settings.sdbg->IsValidEdge(node))
            continue;
        SmallNodeSet parents;
        this->_GetIncomings(node, parents);
        this->settings.sdbg->SetInvalidEdge(node);
        for (uint64_t parent : parents)
            if (this->settings.sdbg->IsValidEdge(parent))
                work.push(parent);
    }
}
void CycleFinder::InvalidateMultiplicityOneNodes() {
    uint64_t invalidated = 0;
    #pragma omp parallel for reduction(+:invalidated)
    for (uint64_t node = 0; node < this->settings.sdbg->size(); node++) {
        if (this->settings.sdbg->EdgeMultiplicity(node) <=1) {
            this->settings.sdbg->SetInvalidEdge(node);
            invalidated += 1;
        }
    }
    std::cout << "  ▸ Pre-filter: removed " << invalidated << " low-multiplicity edges" << std::endl;
}

/**
 * @brief Chunks the start nodes based on their multiplicity for parallel processing.
 */
size_t CycleFinder::ChunkStartNodes(map<int, vector<uint64_t>, greater<int>>& start_nodes_chunked) {
    const int chunk_size = 20000;
    size_t n = this->settings.sdbg->size();
    int threads = static_cast<int>(this->settings.threads);

    // Per-thread chunk collectors to avoid global critical regions
    vector<unordered_map<int, vector<uint64_t>>> local_chunks(threads);
    std::atomic<uint64_t> loaded(0);

    #pragma omp parallel num_threads(threads)
    {
        int tid = omp_get_thread_num();
        #pragma omp for schedule(dynamic, chunk_size)
        for (uint64_t node = 0; node < n; node++) {
            if(!this->settings.sdbg->IsValidEdge(node)) continue;
            size_t edge_indegree = this->settings.sdbg->EdgeIndegree(node);
            uint64_t my_loaded = loaded.fetch_add(1, std::memory_order_relaxed) + 1;
            if (my_loaded % 1000000 == 0) std::cout << "  ▸ Scanning: " << (my_loaded / 1000000) << "M nodes..." << std::endl;
            if (edge_indegree >= 2 && this->settings.sdbg->EdgeMultiplicity(node) > this->settings.cycle_finder_settings.threshold_multiplicity)
            {
                if(this->_HasSelfLoop(node)) continue;
                int reached_depth = 0;

                bool dls = this->DepthLevelSearch(node, node, this->settings.cycle_finder_settings.cycle_max_length, reached_depth);
                if(!dls) continue;
                int log2_mult = static_cast<int>(ceil(log2(double(this->settings.sdbg->EdgeMultiplicity(node)))));
                local_chunks[tid][log2_mult].push_back(node);
            }
        }
    }

    // Merge local chunks into the shared map (serial merge to avoid contention)
    for (int t = 0; t < threads; ++t) {
        for (auto &entry : local_chunks[t]) {
            auto &vec = start_nodes_chunked[entry.first];
            vec.insert(vec.end(), entry.second.begin(), entry.second.end());
        }
    }

    // Dynamic scanning assigns chunks to workers in an unpredictable order.
    // Logical partitions use this order independently of worker assignment.
    for (auto& bucket : start_nodes_chunked) {
        std::sort(bucket.second.begin(), bucket.second.end());
    }
    size_t sum_of_all_quantities_in_all_chunks = 0;
    for (const auto& [key, value] : start_nodes_chunked) {
        std::cout << "  ▸ Bucket mult\u00d7" << (1 << key) << " (log2=" << key << "): " << value.size() << " nodes" << std::endl;
        sum_of_all_quantities_in_all_chunks += value.size();
    }
    return sum_of_all_quantities_in_all_chunks;
} 


/**
 * @brief Finds all cycles in the graph by iterating over chunked start nodes and utilizing parallel processing.
 */
size_t CycleFinder::FindApproximateCRISPRArrays() {
    
    vector<uint64_t> tips = this->CollectTips();
    std::cout << "  ▸ Graph: " << this->settings.sdbg->size() << " nodes, tips: " << tips.size() << std::endl;
    
    if (!this->settings.cycle_finder_settings.low_abundance)
        this->InvalidateMultiplicityOneNodes();
    for ( uint64_t tip : tips) {
        this->RecursiveReduction(tip);
    }
    uint64_t valid_edges = 0;
    #pragma omp parallel for reduction(+:valid_edges)
    for (uint64_t node = 0; node < this->settings.sdbg->size(); node++) {
        if (this->settings.sdbg->IsValidEdge(node)) {
            valid_edges += 1;
        }
    }

    tips = this->CollectTips();
    std::cout << "  ▸ After tip pruning: " << valid_edges << " valid edges, " << tips.size() << " tips remain" << std::endl;
    std::size_t cumulative = 0;
    std::cout << "  ▸ Params: max-len=" << this->settings.cycle_finder_settings.cycle_max_length
              << "  min-len=" << this->settings.cycle_finder_settings.cycle_min_length
              << "  partitions=" << this->settings.cycle_finder_settings.traversal_partitions
              << "  threads=" << this->settings.threads << std::endl;

    map<int, vector<uint64_t>, greater<int>> start_nodes_chunked;
    size_t start_nodes_amount=this->ChunkStartNodes(start_nodes_chunked);
    std::cout << "  ▸ Start nodes: " << start_nodes_amount << std::endl;
    int max_threads = std::min(static_cast<int>(this->settings.threads),
                              this->settings.cycle_finder_settings.traversal_partitions);
    size_t words = (this->settings.sdbg->size() + 63) / 64;
    per_thread_visited.resize(max_threads);
    for (auto &v : per_thread_visited) {
        v.resize(words, 0);
    }
    for (auto nodes_iterator = start_nodes_chunked.begin(); nodes_iterator != start_nodes_chunked.end(); nodes_iterator++) {
        size_t cumulative_at_bucket_start = cumulative;
        auto thread_count = max_threads;
        thread_count = static_cast<int>(std::min(nodes_iterator->second.size(),
                                                 static_cast<size_t>(thread_count)));

        // Per-thread results to avoid concurrent modification of shared unordered_map
        vector<unordered_map<uint64_t, vector<vector<uint64_t>>>> local_results(thread_count);
        vector<size_t> local_cycle_counts(thread_count, 0);

        Zuck::visit_search_partitions(nodes_iterator->second.size(),
            this->settings.cycle_finder_settings.traversal_partitions, thread_count,
            [&](int tid) {
                std::fill(per_thread_visited[tid].begin(), per_thread_visited[tid].end(), 0);
            },
            [&](std::size_t start_node_index, int tid) {
                uint64_t start_node = nodes_iterator->second[start_node_index];
                vector<vector<uint64_t>> cycles = this->FindCycleUtil(start_node);
                local_cycle_counts[tid] += cycles.size();
                if (!cycles.empty()) {
                    local_results[tid][start_node] = std::move(cycles);
                }
            });

        // Merge local results into shared structures (serial merge)
        for (int t = 0; t < thread_count; ++t) {
            for (auto &entry : local_results[t]) {
                this->results[entry.first] = std::move(entry.second);
            }
            cumulative += local_cycle_counts[t];
        }

        malloc_trim(0);
        // Summarize work performed for this multiplicity bucket and avoid printing too often
        size_t cycles_in_bucket = cumulative - cumulative_at_bucket_start;
        std::cout << "  ▸ Bucket log2=" << nodes_iterator->first << ": "
                  << nodes_iterator->second.size() << " nodes, "
                  << cycles_in_bucket << " cycles (total " << cumulative << ")" << std::endl;
    }
        // Completed cycle enumeration
        std::cout << "  ▸ Total: " << cumulative << " cycles across "
              << this->results.size() << " start nodes" << std::endl;
        std::cout << "  ▸ Saturated starts (" << cluster_bounds << "-cycle cap): "
                  << saturated_start_nodes.load(std::memory_order_relaxed)
                  << std::endl;
        std::cout << "  ▸ Step-limited starts (10000000-step cap): "
                  << step_limited_start_nodes.load(std::memory_order_relaxed)
                  << std::endl;

    const std::string cycles_path = this->settings.output_folder + "/cycles.txt";
    this->writeMapToFile(this->results, cycles_path);

    return cumulative;
}

void CycleFinder::writeMapToFile(
    const std::unordered_map<uint64_t, std::vector<std::vector<uint64_t>>>& cycles,
    const std::string& filename)
{
    namespace fs = std::filesystem;
    const fs::path destination(filename);
    const auto check_destination = [&]() {
        const auto status = fs::symlink_status(destination);
        if (fs::is_symlink(status) || (fs::exists(status) && !fs::is_regular_file(status)))
            throw std::runtime_error("Cannot write cycles to: " + filename);
    };
    check_destination();

    // A unique sibling directory owns the temporary file. Closing and renaming
    // only after all serialization succeeds preserves previous cycle evidence
    // on write errors and never follows a destination symlink or hard link.
    struct StagingDirectory {
        fs::path path;
        ~StagingDirectory() {
            if (path.empty()) return;
            std::error_code ignored;
            fs::remove_all(path, ignored);
        }
    } staging;
    const auto parent = destination.has_parent_path() ? destination.parent_path() : fs::path(".");
    const auto nonce = std::chrono::steady_clock::now().time_since_epoch().count();
    for (unsigned int attempt = 0; attempt < 100; ++attempt) {
        const auto candidate = parent / (".Zuck-cycles-" + std::to_string(nonce) + "-" + std::to_string(attempt));
        if (fs::create_directory(candidate)) { staging.path = candidate; break; }
    }
    if (staging.path.empty())
        throw std::runtime_error("Cannot write cycles staging directory: " + filename);
    const auto temporary = staging.path / "cycles.txt";
    std::ofstream out(temporary);
    if (!out.is_open()) throw std::runtime_error("Cannot write cycles to: " + filename);
    // Lossless compact representation: the first k-mer followed by the one
    // new base contributed by each later graph edge. Keep it headerless so a
    // v1.0.0 parser, which already accepts one-token glued records, can replay
    // it too. If a future graph path violates the k-1 overlap invariant, fall
    // back to that record's legacy space-separated k-mers instead of silently
    // changing its sequence.

    const size_t k = this->settings.sdbg->k();
    std::vector<uint8_t> seq(k);
    std::vector<uint64_t> start_nodes;
    start_nodes.reserve(cycles.size());
    for (const auto& entry : cycles) start_nodes.push_back(entry.first);
    std::sort(start_nodes.begin(), start_nodes.end());

    std::string output_buffer;
    output_buffer.reserve(8 * 1024 * 1024);
    auto flush_buffer = [&]() {
        out.write(output_buffer.data(),
                  static_cast<std::streamsize>(output_buffer.size()));
        if (!out)
            throw std::runtime_error("Error writing cycles to: " + filename);
        output_buffer.clear();
    };
    std::size_t legacy_fallback_records = 0;

    auto node_label = [&](uint64_t node) {
        seq.assign(k, 0);
        this->settings.sdbg->GetLabel(node, seq.data());
        std::string label(k, 'N');
        for (size_t j = 0; j < k; ++j) {
            const uint8_t c = seq[j];
            label[j] = (c >= 1 && c <= 4) ? "ACGT"[c - 1] : 'N';
        }
        return label;
    };

    for (uint64_t start_node : start_nodes) {
        const auto& cycle_vec = cycles.at(start_node);
        std::vector<const std::vector<uint64_t>*> ordered_cycles;
        ordered_cycles.reserve(cycle_vec.size());
        for (const auto& cycle : cycle_vec) ordered_cycles.push_back(&cycle);
        std::sort(ordered_cycles.begin(), ordered_cycles.end(),
                  [](const auto* lhs, const auto* rhs) { return *lhs < *rhs; });

        for (const auto* cycle_ptr : ordered_cycles) {
            const auto& cycle = *cycle_ptr;
            if (cycle.empty()) continue;

            std::string previous = node_label(cycle.front());
            std::string record = previous;
            bool compact = true;
            for (size_t i = 1; i < cycle.size(); ++i) {
                std::string current = node_label(cycle[i]);
                if (!std::equal(previous.begin() + 1, previous.end(),
                                current.begin())) {
                    compact = false;
                    break;
                }
                record.push_back(current.back());
                previous = std::move(current);
            }

            if (!compact) {
                ++legacy_fallback_records;
                record.clear();
                for (std::size_t i = 0; i < cycle.size(); ++i) {
                    if (i != 0) record.push_back(' ');
                    record += node_label(cycle[i]);
                }
            }
            record.push_back('\n');
            output_buffer += record;
            if (output_buffer.size() >= 8 * 1024 * 1024) flush_buffer();
        }
    }
    if (!output_buffer.empty()) flush_buffer();
    out.close();
    if (!out)
        throw std::runtime_error("Error closing cycles file: " + filename);
    check_destination();
    fs::rename(temporary, destination);
    if (legacy_fallback_records != 0)
        std::cout << "  ▸ Compact-overlap fallbacks: "
                  << legacy_fallback_records << std::endl;
    std::cout << "  ▸ Cycles saved: " << filename << std::endl;
}
