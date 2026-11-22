"""Production visited predicate/search under real, possibly reduced OpenMP teams."""
import os
from pathlib import Path
import re
import subprocess
import tempfile


def main():
    root = Path(__file__).resolve().parents[1]
    source = (root / 'src/crispr_array/cycle_finder.cpp').read_text()
    header = (root / 'include/crispr_array/cycle_finder.h').read_text()
    method = re.search(r'vector<vector<uint64_t>> CycleFinder::FindCycle\(.*?\n}\n', source, re.S).group()
    predicate = re.search(r'bool CycleFinder::_BackgroundCheck\(.*?\n}\n', source, re.S).group()
    container = re.search(r'struct SmallNodeSet \{.*?\n};', header, re.S).group()
    harness = r'''
#include <algorithm>
#include <atomic>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <map>
#include <set>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include "crispr_array/search_partitions.h"
using namespace std;
namespace phmap {
template<class K,class V> using flat_hash_map = unordered_map<K,V>;
template<class K> using flat_hash_set = unordered_set<K>;
}
__CONTAINER__
struct Graph {
    map<uint64_t,vector<uint64_t>> edges;
    set<uint64_t> invalid;
    set<uint64_t> zero_multiplicity;
    bool IsValidEdge(uint64_t node) const { return !invalid.count(node); }
    size_t EdgeMultiplicity(uint64_t node) const {
        assert(IsValidEdge(node)); // Invalid edges must be rejected before reading evidence.
        return zero_multiplicity.count(node) ? 0 : 30;
    }
};
struct Settings {
    Graph* sdbg;
    struct {int cycle_min_length=0; int cycle_max_length=90;} cycle_finder_settings;
};
class CycleFinder {
public:
    Settings& settings;
    uint16_t cluster_bounds=500;
    vector<vector<uint64_t>> per_thread_visited = vector<vector<uint64_t>>(8, vector<uint64_t>(1));
    atomic<uint64_t> saturated_start_nodes{0}, step_limited_start_nodes{0};
    explicit CycleFinder(Settings& s):settings(s){}
    bool _BackgroundCheck(uint64_t,size_t,uint64_t);
    void _GetOutgoings(uint64_t n,SmallNodeSet& out,size_t mult) {
        for(auto next:settings.sdbg->edges.at(n)) if(_BackgroundCheck(n,mult,next)) out.insert(next);
    }
    void _GetIncomings(uint64_t n,SmallNodeSet& out,size_t mult) {
        for(const auto& [from,edges]:settings.sdbg->edges)
            for(auto next:edges) if(next==n && _BackgroundCheck(n,mult,from)) out.insert(from);
    }
    vector<vector<uint64_t>> FindCycle(uint64_t,vector<uint64_t>&,
        phmap::flat_hash_map<uint64_t,int>&,vector<SmallNodeSet>&,vector<int>&);
    vector<vector<uint64_t>> Run(uint64_t start) {
        vector<uint64_t> path{start};
        phmap::flat_hash_map<uint64_t,int> locks{{start,0}};
        SmallNodeSet next; _GetOutgoings(start,next,30);
        vector<SmallNodeSet> stack{next}; vector<int> backs{90};
        return FindCycle(start,path,locks,stack,backs);
    }
};
__PREDICATE__
__METHOD__
int main() {
    Graph graph{{{0,{1,4}},{1,{0,2}},{2,{3}},{3,{1}},{4,{0}}}};
    Settings settings{&graph,{}};
    CycleFinder predicate_finder(settings);
    assert(predicate_finder._BackgroundCheck(0, 30, 1));
    graph.invalid.insert(1);
    assert(!predicate_finder._BackgroundCheck(0, 30, 1));
    graph.invalid.clear();
    graph.zero_multiplicity.insert(1);
    assert(!predicate_finder._BackgroundCheck(0, 30, 1));
    graph.zero_multiplicity.clear();
    assert(!predicate_finder._BackgroundCheck(1, 30, 1));
    using Results = vector<vector<vector<uint64_t>>>;
    auto collect = [&](const vector<uint64_t>& starts,int partitions,int workers,bool serial) {
        CycleFinder finder(settings);
        Results result(starts.size());
        auto reset = [&](int worker) { fill(finder.per_thread_visited[worker].begin(),
                                           finder.per_thread_visited[worker].end(),0); };
        auto visit = [&](size_t index,int) {
            result[index] = finder.Run(starts[index]);
            sort(result[index].begin(),result[index].end());
        };
        // A serial specification of historical static balanced partitions.
        if(serial) {
            size_t count=min(starts.size(),size_t(partitions)), cursor=0;
            for(size_t p=0;p<count;++p) {
                reset(0);
                size_t length=starts.size()/count+(p<starts.size()%count);
                for(size_t j=0;j<length;++j) visit(cursor++,0);
            }
        } else {
            // Reuse scratch across buckets: each invocation must reset histories.
            Zuck::visit_search_partitions(starts.size(),partitions,workers,reset,visit);
            Zuck::visit_search_partitions(starts.size(),partitions,workers,reset,visit);
        }
        return result;
    };
    assert(collect({0,1},1,1,true)!=collect({0,1},2,1,true));
    size_t checks=0;
    for(int count:{0,1,2,7,8,9,17,31}) {
        vector<uint64_t> starts;
        for(int i=0;i<count;++i) starts.push_back(i%2);
        for(int partitions:{1,2,8}) {
            const auto expected=collect(starts,partitions,1,true);
            for(int workers:{1,2,4,8}) {
                assert(collect(starts,partitions,workers,false)==expected);
                ++checks;
            }
        }
    }
    cout << checks << " search-history comparisons passed\n";
}
'''
    harness = harness.replace('__CONTAINER__', container).replace('__PREDICATE__', predicate).replace('__METHOD__', method)
    with tempfile.TemporaryDirectory(prefix='zuck-partitions-') as temporary:
        cpp = Path(temporary) / 'test.cpp'
        binary = Path(temporary) / 'test'
        cpp.write_text(harness)
        subprocess.run(['c++', '-std=c++17', '-O2', '-fopenmp', '-I', str(root / 'include'),
                        str(cpp), '-o', str(binary)], check=True)
        for overrides in ({'OMP_DYNAMIC': 'FALSE'}, {'OMP_THREAD_LIMIT': '1', 'OMP_DYNAMIC': 'FALSE'},
                          {'OMP_DYNAMIC': 'TRUE', 'OMP_THREAD_LIMIT': '2'}):
            result = subprocess.run([str(binary)], env={**os.environ, **overrides},
                                    capture_output=True, text=True, check=True, timeout=30)
            assert '96 search-history comparisons passed' in result.stdout
            print(overrides, result.stdout.strip())


if __name__ == '__main__':
    main()
