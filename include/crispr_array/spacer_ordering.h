#ifndef INCLUDE_SPACER_ORDERING_H_
#define INCLUDE_SPACER_ORDERING_H_

#include <cstddef>
#include <cstdint>
#include <functional>
#include <limits>
#include <stack>
#include <tuple>
#include <unordered_map>
#include <unordered_set>
#include <vector>

class SDBG;

// Optional graph-based spacer ordering helpers. The production array tables
// remain unordered sequence evidence; these heuristic functions do not infer
// genomic order for the main detection pipeline.

// Minimal adjacency graph used to select relevant reads and cycles.
struct Graph {
    std::unordered_map<std::uint64_t, std::vector<std::uint64_t>> adjacency_list;
    std::unordered_set<std::uint64_t> nodes;

    void add_edge(const std::uint64_t from, const std::uint64_t to) {
        adjacency_list[from].push_back(to);
        nodes.insert(from);
        nodes.insert(to);
    }

    std::size_t edge_count() const {
        std::size_t count = 0;
        for (const auto& kv : adjacency_list) {
            count += kv.second.size();
        }
        return count;
    }
};

inline constexpr std::uint32_t NOT_IN_ANY_CYCLE_INDEX = std::numeric_limits<std::uint32_t>::max();

struct TupleHash {
    std::size_t operator()(const std::tuple<std::uint32_t, std::uint32_t>& t) const {
        return std::hash<std::uint32_t>()(std::get<0>(t)) ^ (std::hash<std::uint32_t>()(std::get<1>(t)) << 1);
    }
};

void find_strongly_connected_components_dfs(
    const std::uint64_t node,
    const SDBG& sdbg,
    std::unordered_map<std::uint64_t, int>& index_map,
    std::unordered_map<std::uint64_t, int>& lowlink_map,
    std::unordered_set<std::uint64_t>& on_stack,
    std::stack<std::uint64_t>& stack,
    std::vector<std::vector<std::uint64_t>>& components,
    int& index_counter
);

std::vector<std::vector<std::uint64_t>> find_strongly_connected_components(const SDBG& sdbg);

void keep_crispr_regions_extended_by_k(
    SDBG& sdbg,
    const std::size_t& k,
    const std::vector<std::vector<std::uint64_t>>& cycles
);

std::vector<Graph> divide_graph_into_subgraphs(const SDBG& sdbg);

std::vector<Graph> get_crispr_regions_extended_by_k(
    SDBG& sdbg,
    const std::size_t& k,
    const std::vector<std::vector<std::uint64_t>>& cycles
);

// Return unique cycle IDs in ascending order.
std::vector<std::uint32_t> get_all_cycle_indices(
    const std::unordered_map<std::uint64_t, std::uint32_t>& node_to_cycle_map
);

// Select nonempty reads with at least one endpoint in the graph.
std::vector<std::vector<std::uint64_t>> get_relevant_reads(
    const Graph& graph,
    const std::vector<std::vector<std::uint64_t>>& all_reads
);

std::vector<std::vector<std::uint64_t>> get_relevant_cycles(
    const Graph& graph,
    const std::vector<std::vector<std::uint64_t>>& all_cycles
);

// Keep a covering subset in input order.
void get_minimum_cycles_for_full_coverage(std::vector<std::vector<std::uint64_t>>& cycles);

// Universe IDs must be dense in [0, universe.size()).
std::vector<std::size_t> solve_min_cover_problem(
    const std::unordered_set<std::uint32_t>& universe,
    const std::vector<std::vector<std::uint32_t>>& sets
);

// Map nodes present in exactly one cycle; repetition within that cycle is fine.
std::unordered_map<std::uint64_t, std::uint32_t> get_node_to_unique_cycle_map(
    const std::vector<std::vector<std::uint64_t>>& cycles
);

std::vector<std::tuple<std::uint32_t, std::uint32_t>> every_possible_combination(const std::vector<std::uint32_t>& v);

std::vector<std::tuple<std::uint32_t, std::uint32_t>> generate_constraints_from_read(
    const Graph& graph,
    const std::vector<std::uint64_t>& read,
    const std::unordered_map<std::uint64_t, std::uint32_t>& node_to_cycle_map
);

// Only adjacent transitions with exactly one unmapped endpoint are emitted.
std::vector<std::tuple<std::uint32_t, std::uint32_t>> generate_out_of_cycles_constraints_from_read(
    const Graph& graph,
    const std::vector<std::uint64_t>& read,
    const std::unordered_map<std::uint64_t, std::uint32_t>& node_to_cycle_map
);

std::vector<std::tuple<std::uint32_t, std::uint32_t>> generate_constraints(
    const Graph& graph,
    const std::vector<std::vector<std::uint64_t>>& reads,
    const std::unordered_map<std::uint64_t, std::uint32_t>& node_to_cycle_map
);

std::vector<std::tuple<std::uint32_t, std::uint32_t>> get_maximal_spanning_tree(
    const std::vector<std::tuple<std::uint32_t, std::uint32_t>>& edges
);

void resolve_cycles_greedy(
    std::vector<std::tuple<std::uint32_t, std::uint32_t>>& constraints,
    std::unordered_map<std::uint32_t, int>& heuristic_node_values
);

void apply_topological_sort(
    std::vector<std::uint32_t>& possible_start_nodes,
    const std::unordered_map<std::uint32_t, int>& node_affection_to_start,
    std::unordered_map<std::uint32_t, int>& heuristic_node_values,
    std::unordered_map<std::tuple<std::uint32_t, std::uint32_t>, int, TupleHash>& edges,
    std::vector<std::uint32_t>& total_order,
    float& confidence
);

// Constraints must be acyclic. Confidence is a heuristic evidence index;
// empty or unsupported choices have zero confidence. Ties use ascending IDs.
std::vector<std::uint32_t> solve_constraints_with_topological_sort(
    const std::vector<std::tuple<std::uint32_t, std::uint32_t>>& constraints,
    std::unordered_map<std::uint32_t, int>& heuristic_node_values,
    const std::vector<std::uint32_t>& nodes,
    float& confidence
);

std::vector<std::uint32_t> order_cycles(
    const Graph& graph,
    const std::vector<std::vector<std::uint64_t>>& reads,
    const std::vector<std::vector<std::uint64_t>>& cycles,
    float& confidence_cycle_resolution,
    float& confidence_topological_sort
);

std::vector<std::vector<std::uint64_t>> get_ordered_cycles(
    const std::vector<std::uint32_t>& cycle_order,
    const std::vector<std::vector<std::uint64_t>>& cycles
);

#endif
