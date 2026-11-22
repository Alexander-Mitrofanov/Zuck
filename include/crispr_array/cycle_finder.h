#ifndef CYCLE_FINDER_H
#define CYCLE_FINDER_H

#include <atomic>
#include <cstddef>
#include <cstdint>
#include <map>
#include <string>
#include <unordered_map>
#include <vector>
#include <parallel_hashmap/phmap.h>
#include "core/settings.h"

// DNA graph vertices have at most four neighbors. Erasure swaps in the last
// element, preserving the historical traversal order without heap allocation.
struct SmallNodeSet {
    uint64_t nodes[4];
    uint8_t  count = 0;
    void insert(uint64_t n)  { nodes[count++] = n; }
    bool empty()       const { return count == 0; }
    void erase(uint64_t n) {
        for (uint8_t i = 0; i < count; ++i)
            if (nodes[i] == n) { nodes[i] = nodes[--count]; return; }
    }
    uint64_t*       begin()       { return nodes; }
    uint64_t*       end()         { return nodes + count; }
    const uint64_t* begin() const { return nodes; }
    const uint64_t* end()   const { return nodes + count; }
};

// Enumerates bounded graph cycles. Construction runs the configured search;
// results retain graph-edge paths and cycles.txt stores their DNA sequences.
class CycleFinder {
private:
    Settings& settings;
    static constexpr std::size_t MAX_EDGE_COUNT = 4;
    std::uint16_t cluster_bounds = 500;
    // One scratch bitmap per worker, reset for each logical search partition.
    std::vector<std::vector<std::uint64_t>> per_thread_visited;
    std::atomic<std::uint64_t> saturated_start_nodes{0};
    std::atomic<std::uint64_t> step_limited_start_nodes{0};

    bool _HasSelfLoop(std::uint64_t node);
    bool _BackgroundCheck(std::uint64_t original_node,
                          std::size_t repeat_multiplicity, std::uint64_t neighbor);
    void _GetOutgoings(std::uint64_t node, SmallNodeSet& neighbors,
                       std::size_t repeat_multiplicity);
    void _GetIncomings(std::uint64_t node, SmallNodeSet& neighbors,
                       std::size_t repeat_multiplicity);
    void _GetOutgoings(std::uint64_t node, SmallNodeSet& neighbors);
    void _GetIncomings(std::uint64_t node, SmallNodeSet& neighbors);

public:
    explicit CycleFinder(Settings& settings);
    std::unordered_map<std::uint64_t, std::vector<std::vector<std::uint64_t>>> results;

    std::vector<std::uint64_t> CollectTips();
    void RecursiveReduction(std::uint64_t tip);
    void InvalidateMultiplicityOneNodes();
    std::size_t ChunkStartNodes(
        std::map<int, std::vector<std::uint64_t>, std::greater<int>>& buckets);
    bool DepthLevelSearch(std::uint64_t start, std::uint64_t target,
                          int limit, int& reached_depth);
    std::vector<std::vector<std::uint64_t>> FindCycle(
        std::uint64_t start, std::vector<std::uint64_t>& path,
        phmap::flat_hash_map<std::uint64_t, int>& locks,
        std::vector<SmallNodeSet>& stack, std::vector<int>& backtrack_lengths);
    std::vector<std::vector<std::uint64_t>> FindCycleUtil(std::uint64_t start);
    std::size_t FindApproximateCRISPRArrays();
    void writeMapToFile(
        const std::unordered_map<std::uint64_t, std::vector<std::vector<std::uint64_t>>>& cycles,
        const std::string& filename);
};
#endif
