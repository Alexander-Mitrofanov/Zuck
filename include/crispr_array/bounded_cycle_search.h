#ifndef ZUCK_BOUNDED_CYCLE_SEARCH_H
#define ZUCK_BOUNDED_CYCLE_SEARCH_H

#include <algorithm>
#include <cstdint>
#include <unordered_set>
#include <utility>
#include <vector>

namespace Zuck {

template <typename VisitedSet = std::unordered_set<uint64_t>>
struct BoundedCycleWorkspace {
    std::vector<std::pair<uint64_t, int>> queue;
    VisitedSet visited;
};

// Decide whether target is reachable within limit edges. For start == target,
// require a directed cycle with at least two edges (ignore self loops).
// Breadth-first visitation is essential: visiting a node first along a longer
// DFS route must not suppress a shorter route with enough remaining depth.
// Neighbors must return only valid graph edges. Workspace capacity is reused.
template <typename Neighbors, typename VisitedSet>
bool HasBoundedPath(uint64_t start, uint64_t target, int limit,
                    Neighbors&& neighbors,
                    BoundedCycleWorkspace<VisitedSet>& workspace,
                    int& reached_depth) {
    workspace.queue.clear();
    workspace.visited.clear();
    reached_depth = 0;
    if (limit < 1 || (start == target && limit < 2)) return false;

    workspace.queue.emplace_back(start, 0);
    workspace.visited.insert(start);
    for (std::size_t head = 0; head < workspace.queue.size(); ++head) {
        const auto current = workspace.queue[head];
        const int next_depth = current.second + 1;
        if (next_depth > limit) continue;
        reached_depth = std::max(reached_depth, current.second);
        for (uint64_t next : neighbors(current.first)) {
            // Test closure before the depth bound and visited set. In
            // particular, a cycle of exactly limit edges is admissible.
            if (next == target && (start != target || next_depth > 1)) {
                reached_depth = next_depth;
                return true;
            }
            if (next_depth < limit && workspace.visited.insert(next).second)
                workspace.queue.emplace_back(next, next_depth);
        }
    }
    return false;
}

}  // namespace Zuck

#endif
