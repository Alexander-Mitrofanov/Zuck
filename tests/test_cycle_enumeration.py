"""Exercise the production enumerator against a small in-memory graph oracle.

The production method and its fixed-capacity neighbor container are compiled
verbatim. Only SDBG adjacency/multiplicity and OpenMP thread lookup are replaced,
so the test does not need a gigabase graph or duplicate the search algorithm.
"""
import os
from pathlib import Path
import re
import shlex
import subprocess
import tempfile


def main() -> None:
    root = Path(__file__).resolve().parents[1]
    source = (root / "src/crispr_array/cycle_finder.cpp").read_text()
    header = (root / "include/crispr_array/cycle_finder.h").read_text()
    method = re.search(
        r"vector<vector<uint64_t>> CycleFinder::FindCycle\(.*?\n}\n",
        source, re.S,
    )
    nodes = re.search(r"struct SmallNodeSet \{.*?\n};", header, re.S)
    if method is None or nodes is None:
        raise AssertionError("production enumerator/container could not be located")
    harness = r'''
#include <algorithm>
#include <atomic>
#include <cstdint>
#include <functional>
#include <iostream>
#include <map>
#include <random>
#include <set>
#include <stdexcept>
#include <unordered_map>
#include <unordered_set>
#include <vector>
using namespace std;
namespace phmap {
template <typename K, typename V> using flat_hash_map = unordered_map<K, V>;
template <typename K> using flat_hash_set = unordered_set<K>;
}
int omp_get_thread_num() { return 0; }
__NODE_CONTAINER__
struct Graph {
    map<uint64_t, vector<uint64_t>> edges;
    size_t EdgeMultiplicity(uint64_t) const { return 30; }
};
struct Settings {
    Graph* sdbg;
    struct { int cycle_min_length = 0; int cycle_max_length = 90; }
        cycle_finder_settings;
};
class CycleFinder {
public:
    Settings& settings;
    uint16_t cluster_bounds = 500;
    // Sparse words allow testing 64-bit IDs without allocating a huge graph.
    vector<unordered_map<uint64_t, uint64_t>> per_thread_visited{1};
    atomic<uint64_t> saturated_start_nodes{0};
    atomic<uint64_t> step_limited_start_nodes{0};
    explicit CycleFinder(Settings& s) : settings(s) {}
    void _GetOutgoings(uint64_t n, SmallNodeSet& result, size_t) {
        for (uint64_t next : settings.sdbg->edges[n])
            if (next != n) result.insert(next);
    }
    void _GetIncomings(uint64_t n, SmallNodeSet& result, size_t) {
        for (const auto& row : settings.sdbg->edges)
            for (uint64_t next : row.second)
                if (next == n && row.first != n) result.insert(row.first);
    }
    vector<vector<uint64_t>> FindCycle(uint64_t, vector<uint64_t>&,
        phmap::flat_hash_map<uint64_t, int>&, vector<SmallNodeSet>&,
        vector<int>&);
    vector<vector<uint64_t>> Run(uint64_t start) {
        vector<uint64_t> path{start};
        phmap::flat_hash_map<uint64_t, int> locks{{start, 0}};
        SmallNodeSet next;
        _GetOutgoings(start, next, 30);
        vector<SmallNodeSet> stack{next};
        vector<int> backs{settings.cycle_finder_settings.cycle_max_length};
        return FindCycle(start, path, locks, stack, backs);
    }
};
__PRODUCTION_METHOD__
static void Require(bool ok, const char* message) {
    if (!ok) throw runtime_error(message);
}
static set<vector<uint64_t>> Oracle(const Graph& graph, uint64_t start,
                                    int maximum) {
    set<vector<uint64_t>> cycles;
    vector<uint64_t> path{start};
    function<void(uint64_t)> dfs = [&](uint64_t node) {
        const auto row = graph.edges.find(node);
        if (row == graph.edges.end()) return;
        for (uint64_t next : row->second) {
            if (next == start && path.size() > 1) cycles.insert(path);
            else if (path.size() < static_cast<size_t>(maximum) &&
                     find(path.begin(), path.end(), next) == path.end()) {
                path.push_back(next); dfs(next); path.pop_back();
            }
        }
    };
    dfs(start);
    return cycles;
}
static void Check(Graph graph, uint64_t start, int maximum, uint16_t cap) {
    Settings settings{&graph, {}};
    settings.cycle_finder_settings.cycle_max_length = maximum;
    CycleFinder finder(settings);
    finder.cluster_bounds = cap;
    const auto expected = Oracle(graph, start, maximum);
    Require(expected.size() <= cap, "oracle fixture exceeds test cap");
    const auto rows = finder.Run(start);
    const set<vector<uint64_t>> actual(rows.begin(), rows.end());
    Require(rows.size() == actual.size(), "enumerator emitted duplicate paths");
    Require(actual == expected, "enumerator disagrees with exhaustive oracle");
}
int main() {
    // Revisiting 1 after exploring 2 must not emit (0,1) a second time.
    Check(Graph{{{0, {1}}, {1, {0, 2}}, {2, {0}}}}, 0, 4, 500);
    // Duplicate closures previously filled the four-record cap before 0->2->0.
    Check(Graph{{{0, {1, 2}}, {1, {0, 3, 4}}, {2, {0}},
                 {3, {0}}, {4, {0}}}}, 0, 5, 4);
    // Published CYCLE_SEARCH counterexample: the current implementation differs
    // from that algorithm and must retain all four cycles in this graph.
    Check(Graph{{{0, {3, 4}}, {1, {3, 4}}, {2, {0, 1}},
                 {3, {0, 1}}, {4, {2}}}}, 0, 5, 500);
    // Revisiting a merge along a longer path requires relaxation of its real
    // 64-bit ID. An int temporary leaves the high-ID vertex incorrectly locked.
    constexpr uint64_t wide = uint64_t{1} << 40;
    Check(Graph{{{wide, {wide + 1, wide + 2}}, {wide + 1, {wide + 3}},
                 {wide + 2, {wide + 4}}, {wide + 4, {wide + 3}},
                 {wide + 3, {wide}}}}, wide, 5, 500);
    mt19937 rng(918572);
    size_t checked = 0;
    for (int trial = 0; trial < 10000; ++trial) {
        Graph graph;
        int incoming[7] = {};
        bool valid = true;
        for (uint64_t from = 0; from < 7; ++from) {
            for (uint64_t to = 0; to < 7; ++to) {
                if (from != to && rng() % 4 == 0) {
                    graph.edges[from].push_back(to);
                    if (++incoming[to] > 4) valid = false;
                }
            }
            if (graph.edges[from].size() > 4) valid = false;
        }
        if (!valid) continue;
        Check(graph, rng() % 7, 2 + rng() % 7, 500);
        ++checked;
    }
    cout << "Production enumeration passed " << checked
         << " oracle graphs and duplicate/cap regressions\n";
}
'''
    harness = harness.replace("__NODE_CONTAINER__", nodes.group(0)).replace(
        "__PRODUCTION_METHOD__", method.group(0)
    )
    with tempfile.TemporaryDirectory(prefix="zuck-enumeration-test-") as scratch:
        binary = str(Path(scratch) / "test-enumeration")
        compiler = shlex.split(os.environ.get("CXX", "c++"))
        subprocess.run(
            compiler + ["-std=c++17", "-O2", "-Wall", "-Wextra", "-Werror",
                        "-x", "c++", "-", "-o", binary],
            input=harness, text=True, check=True,
        )
        subprocess.run([binary], check=True)


if __name__ == "__main__":
    main()
