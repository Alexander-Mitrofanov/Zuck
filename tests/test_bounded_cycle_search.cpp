#include "crispr_array/bounded_cycle_search.h"

#include <algorithm>
#include <functional>
#include <iostream>
#include <random>
#include <stdexcept>
#include <unordered_map>
#include <vector>

using Graph = std::unordered_map<uint64_t, std::vector<uint64_t>>;

static void Require(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}

static bool Search(const Graph& graph, uint64_t start, uint64_t target,
                   int limit, Zuck::BoundedCycleWorkspace<>& workspace) {
    const std::vector<uint64_t> empty;
    auto neighbors = [&](uint64_t node) -> const std::vector<uint64_t>& {
        auto found = graph.find(node);
        return found == graph.end() ? empty : found->second;
    };
    int reached_depth = -1;
    return Zuck::HasBoundedPath(start, target, limit, neighbors,
                                   workspace, reached_depth);
}

// Independent exhaustive simple-path oracle for small directed graphs.
static bool Oracle(const Graph& graph, uint64_t start, int limit) {
    std::unordered_set<uint64_t> path{start};
    std::function<bool(uint64_t, int)> dfs = [&](uint64_t node, int depth) {
        if (depth >= limit) return false;
        auto found = graph.find(node);
        if (found == graph.end()) return false;
        for (uint64_t next : found->second) {
            if (next == start && depth + 1 > 1) return true;
            if (!path.insert(next).second) continue;
            if (dfs(next, depth + 1)) return true;
            path.erase(next);
        }
        return false;
    };
    return dfs(start, 0);
}

int main() {
    Zuck::BoundedCycleWorkspace<> workspace;
    const Graph triangle{{0, {1}}, {1, {2}}, {2, {0}}};
    Require(Search(triangle, 0, 0, 3, workspace),
            "cycle exactly at maximum length must survive prefilter");
    Require(!Search(triangle, 0, 0, 2, workspace),
            "cycle exceeding maximum length must not pass");

    // LIFO visitation reaches node 5 by 0->2->3->4->5 before 0->1->5.
    // A node-only visited set then loses the valid four-edge return path.
    Graph merged{{0, {1, 2}}, {1, {5}}, {2, {3}}, {3, {4}},
                 {4, {5}}, {5, {6}}, {6, {0}}};
    Require(Search(merged, 0, 0, 5, workspace),
            "long first route must not suppress a shorter return path");
    std::reverse(merged[0].begin(), merged[0].end());
    Require(Search(merged, 0, 0, 5, workspace),
            "reachability must not depend on outgoing edge order");

    const Graph self_loop{{0, {0}}};
    Require(!Search(self_loop, 0, 0, 10, workspace),
            "a self loop alone is not an eligible cycle");
    const Graph sink{{0, {1}}};
    Require(Search(sink, 0, 1, 1, workspace),
            "a reachable target with no outgoing edges must be accepted");
    Require(!Search(sink, 0, 0, 10, workspace),
            "workspace reuse must not retain prior reachability state");

    constexpr uint64_t large = (uint64_t{1} << 40);
    const Graph wide_ids{{large, {large + 1}}, {large + 1, {large}}};
    Require(Search(wide_ids, large, large, 2, workspace),
            "node identifiers must retain all 64 bits");

    std::mt19937 rng(71843);
    std::size_t comparisons = 0;
    for (int trial = 0; trial < 500; ++trial) {
        Graph graph;
        for (uint64_t from = 0; from < 7; ++from)
            for (uint64_t to = 0; to < 7; ++to)
                if (rng() % 5 == 0) graph[from].push_back(to);
        for (uint64_t start = 0; start < 7; ++start) {
            for (int limit = 0; limit <= 8; ++limit) {
                Require(Search(graph, start, start, limit, workspace) ==
                            Oracle(graph, start, limit),
                        "bounded search disagrees with exhaustive oracle");
                ++comparisons;
            }
        }
    }
    std::cout << "Bounded cycle search passed " << comparisons
              << " oracle comparisons and targeted regressions\n";
}
