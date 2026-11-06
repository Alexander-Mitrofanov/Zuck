#include "crispr_array/cycle_finder.h"
#include "crispr_array/bounded_cycle_search.h"
#include "crispr_array/filters.h"
#include "core/settings.h"
#include <stdexcept>
#include <fstream>
#include <stack>

// Parallel hashmap for better performance in DLS
#include <parallel_hashmap/phmap.h>

/**
 * @file cycle_finder.cpp
 * @brief Implementation of functions for cycle detection and analysis in a sequence graph.
 * - Checking if any incoming edge of a node is not equal to the node itself.
 * - Performing a background check on a neighbor node to determine if it meets certain criteria.
 * - Getting the outgoing edges of a node that pass the background check.
 * - Retrieving the incoming edges of a node that pass the background check.
 * - Finding cycles in the graph.
 * - Performing a depth-limited search to determine if a path exists between two nodes within a certain depth.
 * - Collecting tips (nodes with no outgoing edges) in the graph.
 * - Performing recursive reduction to invalidate edges.
 * - Invalidating nodes with edge multiplicity of one.
 * - Chunking start nodes based on their multiplicity for parallel processing.
 * - Finding all cycles in the graph by iterating over chunked start nodes and utilizing parallel processing.
 */


/**
 * @brief Checks if any incoming edge of a node is not equal to the node itself.
 */
bool CycleFinder::_IncomingNotEqualToCurrentNode(uint64_t node, size_t edge_indegree) {
    uint64_t incomings[edge_indegree];
    this->settings.sdbg->IncomingEdges(node, incomings);
    for (const auto& incoming : incomings)
        if (node==incoming)
            return true;
    return false;
}
/**
 * @brief Performs a background check on a neighbor node to determine if it meets certain criteria.
 */
bool CycleFinder::_BackgroundCheck(uint64_t original_node, size_t repeat_multiplicity, uint64_t neighbor_node) {
    auto neighbor_node_multiplicity = this->settings.sdbg->EdgeMultiplicity(neighbor_node);
    // Use lock-free visited test (atomic bitset) to avoid contention
    size_t idx = neighbor_node >> 6;
    uint64_t mask = 1ULL << (neighbor_node & 63);
    int tid = omp_get_thread_num();
    if ((per_thread_visited[tid][idx] & mask) != 0) {
        return false;
    }
    if (repeat_multiplicity / neighbor_node_multiplicity > 500) {
        return false;
    }
    if (original_node == neighbor_node) {
        return false;
    }
    return true; 
}


/**
 * @brief Gets the outgoing edges of a node that pass the background check.
 */
void CycleFinder::_GetOutgoings(uint64_t node, SmallNodeSet& outgoings_set, size_t repeat_multiplicity) {
   
    int edge_outdegree = this->settings.sdbg->EdgeOutdegree(node);
    if (edge_outdegree == 0 || !this->settings.sdbg->IsValidEdge(node)) {
        return;
    }
     uint64_t outgoings[edge_outdegree];
    int flag = this->settings.sdbg->OutgoingEdges(node, outgoings);
    if(flag!=-1)    
        for (const auto& outgoing : outgoings)
            if (this->_BackgroundCheck(node, repeat_multiplicity, outgoing) && this->settings.sdbg->IsValidEdge(outgoing))
                outgoings_set.insert(outgoing);
    
    
}
/**
 * @brief Retrieves the incoming edges of a node that pass the background check.
 */
void CycleFinder::_GetIncomings(uint64_t node, SmallNodeSet& incomings_set, size_t repeat_multiplicity) {
  
    int edge_indegree = this->settings.sdbg->EdgeIndegree(node);
    if (edge_indegree == 0 || !this->settings.sdbg->IsValidEdge(node)) {
        return;
    }
    uint64_t incomings[edge_indegree];
    int flag = this->settings.sdbg->IncomingEdges(node, incomings);
    if (flag!=-1)
        for (const auto& incoming : incomings)
            if (this->_BackgroundCheck(node, repeat_multiplicity, incoming) && this->settings.sdbg->IsValidEdge(incoming))
                incomings_set.insert(incoming);
}
// ## START: HELPER FUNCTIONS FOR DLS ##
/**
 * @brief Gets the outgoing edges of a node that pass the background check.
 */
void CycleFinder::_GetOutgoings(uint64_t node, SmallNodeSet& outgoings_set) {
   
    int edge_outdegree = this->settings.sdbg->EdgeOutdegree(node);
    if (edge_outdegree == 0 || !this->settings.sdbg->IsValidEdge(node)) {
        return;
    }
    uint64_t outgoings[edge_outdegree];
    int flag = this->settings.sdbg->OutgoingEdges(node, outgoings);
    if(flag!=-1)
        for (const auto& outgoing : outgoings)
            if (this->settings.sdbg->IsValidEdge(outgoing))
                outgoings_set.insert(outgoing);
    
    
}
/**
 * @brief Retrieves the incoming edges of a node that pass the background check.
 */
void CycleFinder::_GetIncomings(uint64_t node, SmallNodeSet& incomings_set) {
  
    int edge_indegree = this->settings.sdbg->EdgeIndegree(node);
    if (edge_indegree == 0 || !this->settings.sdbg->IsValidEdge(node)) {
        return;
    }
    uint64_t incomings[edge_indegree];
    int flag = this->settings.sdbg->IncomingEdges(node, incomings);
    if (flag!=-1)
        for (const auto& incoming : incomings)
            if (this->settings.sdbg->IsValidEdge(incoming))
                incomings_set.insert(incoming);
}
// ## END: HELPER FUNCTIONS FOR DLS ##

/**
 * @brief Finds the cycles in the graph using settings configuration.
 * 
 * @param settings The application Settings containing SDBG pointer and cycle finder configuration.
 */
CycleFinder::CycleFinder(Settings& settings)
    : settings(settings), cluster_bounds(500) {
    // Ensure settings.sdbg is valid; this should be set by caller before constructing CycleFinder
    if (this->settings.sdbg == nullptr) {
        throw std::runtime_error("CycleFinder requires settings.sdbg to be set to a valid SDBG instance");
    }
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

    // Mark visited nodes atomically (no global critical section)
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
    static thread_local my_mcaat::BoundedCycleWorkspace<
        phmap::flat_hash_set<uint64_t>> workspace;
    auto neighbors = [this](uint64_t node) {
        SmallNodeSet outgoings;
        this->_GetOutgoings(node, outgoings);
        return outgoings;
    };
    return my_mcaat::HasBoundedPath(start, target, limit, neighbors,
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
    if(!this->settings.cycle_finder_settings.low_abundance){
        this->InvalidateMultiplicityOneNodes();
    }
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
                if(this->_IncomingNotEqualToCurrentNode(node, edge_indegree)) continue;
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
    // Sort before the static traversal partition so a fixed thread count has
    // a reproducible start order and visited-state history.
    for (auto& bucket : start_nodes_chunked) {
        std::sort(bucket.second.begin(), bucket.second.end());
    }
   //writeStartNodesToFile(start_nodes_chunked, "start_nodes.txt");
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
int CycleFinder::FindApproximateCRISPRArrays()
 {
    
    vector<uint64_t> tips = this->CollectTips();
    std::cout << "  ▸ Graph: " << this->settings.sdbg->size() << " nodes, tips: " << tips.size() << std::endl;
    
    this->InvalidateMultiplicityOneNodes();
    for ( uint64_t tip : tips) {
        this->RecursiveReduction(tip);
    }
    int valid_edges = 0;
    #pragma omp parallel for reduction(+:valid_edges)
    for (uint64_t node = 0; node < this->settings.sdbg->size(); node++) {
        if (this->settings.sdbg->IsValidEdge(node)) {
            valid_edges += 1;
        }
    }

    tips = this->CollectTips();
    std::cout << "  ▸ After tip pruning: " << valid_edges << " valid edges, " << tips.size() << " tips remain" << std::endl;
    // struct mallinfo mem_info = mallinfo(); // deprecated
    // size_t graph_mem_info = mem_info.uordblks; // unused
    int cumulative = 0;
    std::cout << "  ▸ Params: max-len=" << this->settings.cycle_finder_settings.cycle_max_length
              << "  min-len=" << this->settings.cycle_finder_settings.cycle_min_length
              << "  threads=" << this->settings.threads << std::endl;
        std::unordered_map<uint64_t, std::vector<std::vector<uint64_t>>> all_cycles;

    map<int, vector<uint64_t>, greater<int>> start_nodes_chunked;
    size_t start_nodes_amount=this->ChunkStartNodes(start_nodes_chunked);
    std::cout << "  ▸ Start nodes: " << start_nodes_amount << std::endl;
    size_t counter = 0;
    int max_threads = static_cast<int>(this->settings.threads);
    size_t words = (this->settings.sdbg->size() + 63) / 64;
    per_thread_visited.resize(max_threads);
    for (auto &v : per_thread_visited) {
        v.resize(words, 0);
    }
    for (auto nodes_iterator = start_nodes_chunked.begin(); nodes_iterator != start_nodes_chunked.end(); nodes_iterator++) {
        // Reset visited bitmaps so nodes from a higher-multiplicity bucket don't
        // block traversal in lower-multiplicity buckets.
        for (auto &v : per_thread_visited)
            std::fill(v.begin(), v.end(), 0);
        size_t cumulative_at_bucket_start = cumulative;
        auto thread_count = static_cast<int>(this->settings.threads);
        if (static_cast<int>(nodes_iterator->second.size()) < thread_count)
            thread_count = nodes_iterator->second.size();

        // Per-thread results to avoid concurrent modification of shared unordered_map
        vector<unordered_map<uint64_t, vector<vector<uint64_t>>>> local_results(thread_count);
        vector<size_t> local_cycle_counts(thread_count, 0);
        vector<size_t> local_processed(thread_count, 0);

        #pragma omp parallel num_threads(thread_count)
        {
            int tid = omp_get_thread_num();
            #pragma omp for schedule(static)
            for (uint64_t start_node_index = 0; start_node_index < nodes_iterator->second.size(); start_node_index++) {
                uint64_t start_node = nodes_iterator->second[start_node_index];
                vector<vector<uint64_t>> cycles = this->FindCycleUtil(start_node);
                local_cycle_counts[tid] += cycles.size();
                local_processed[tid] += 1;
                if (!cycles.empty()) {
                    local_results[tid][start_node] = std::move(cycles);
                }
            }
        }

        // Merge local results into shared structures (serial merge)
        for (int t = 0; t < thread_count; ++t) {
            for (auto &entry : local_results[t]) {
                this->results[entry.first] = std::move(entry.second);
            }
            cumulative += local_cycle_counts[t];
            counter += local_processed[t];
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
    std::ofstream out(filename);
    if (!out.is_open()) {
        throw std::runtime_error("Cannot write cycles to: " + filename);
    }
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
    if (legacy_fallback_records != 0)
        std::cout << "  ▸ Compact-overlap fallbacks: "
                  << legacy_fallback_records << std::endl;
    std::cout << "  ▸ Cycles saved: " << filename << std::endl;
}

CycleFinder::~CycleFinder() {}
