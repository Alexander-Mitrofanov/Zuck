#include "crispr_array/spacer_ordering.h"

#include <cmath>
#include <iostream>
#include <stdexcept>

namespace {
void require(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}
}

int main() {
    Graph graph;
    graph.add_edge(1, 2);
    const std::vector<std::vector<uint64_t>> reads{{}, {8, 1}, {2, 9}, {7, 8}};
    require(get_relevant_reads(graph, reads) ==
        std::vector<std::vector<uint64_t>>{{8, 1}, {2, 9}}, "Ignore empty and unrelated reads");
    require(get_relevant_cycles(graph, {{}, {1, 2}, {1, 8}}) ==
        std::vector<std::vector<uint64_t>>{{1, 2}}, "Empty paths are not cycles");

    const auto owners = get_node_to_unique_cycle_map({{1, 1, 2}, {2, 3}, {2, 4}});
    require(owners.size() == 3 && owners.at(1) == 0 && owners.at(3) == 1 && owners.at(4) == 2,
            "Unique ownership counts cycles, not repeated node occurrences");
    require(get_all_cycle_indices(owners) == std::vector<uint32_t>{0, 1, 2},
            "Cycle IDs must have deterministic order");
    const std::unordered_map<uint64_t, uint32_t> mapped{{1, 0}, {2, 0}, {3, 1}, {4, 1}};
    require(generate_constraints_from_read(graph, {1, 2, 3, 4}, mapped) ==
        std::vector<std::tuple<uint32_t, uint32_t>>{{0, 1}},
        "A contiguous run from one cycle contributes one ordering observation");
    require(generate_out_of_cycles_constraints_from_read(graph, {}, mapped).empty(),
            "Empty reads have no outside transitions");
    require(generate_out_of_cycles_constraints_from_read(graph, {9, 1, 2, 3, 9}, mapped) ==
        std::vector<std::tuple<uint32_t, uint32_t>>{{NOT_IN_ANY_CYCLE_INDEX, 0}, {1, NOT_IN_ANY_CYCLE_INDEX}},
        "Keep actual transitions entering and leaving the mapped cycles");
    require(generate_out_of_cycles_constraints_from_read(graph, {1, 3}, mapped).empty(),
            "Inside-to-inside transitions are not outside evidence");

    std::unordered_map<uint32_t, int> heuristic;
    float confidence = -1.0f;
    require(solve_constraints_with_topological_sort({}, heuristic, {}, confidence).empty() && confidence == 0.0f,
            "Empty topological order has finite zero confidence");
    require(solve_constraints_with_topological_sort({}, heuristic, {2, 1, 0}, confidence) ==
        std::vector<uint32_t>{0, 1, 2} && confidence == 0.0f,
        "Unconstrained order is deterministic and has no confidence evidence");
    require(solve_constraints_with_topological_sort({{0, 1}, {1, 2}}, heuristic, {0, 1, 2}, confidence) ==
        std::vector<uint32_t>{0, 1, 2} && std::isfinite(confidence),
        "Topological ordering respects observed edges");
    float resolution = -1.0f;
    const auto ordered = order_cycles(graph, {}, {{1, 2}, {3, 4}}, resolution, confidence);
    require(ordered == std::vector<uint32_t>{0, 1} && resolution == 0.0f && confidence == 0.0f,
            "No ordering constraints must not produce NaN confidences");

    std::vector<std::vector<uint64_t>> redundant(6, std::vector<uint64_t>{1, 2, 3});
    get_minimum_cycles_for_full_coverage(redundant);
    require(redundant.size() == 1 && redundant.front() == std::vector<uint64_t>{1, 2, 3},
            "Coverage compaction must remove every unselected cycle");
    std::cout << "Spacer ordering, empty inputs, confidence and coverage regressions passed\n";
}
