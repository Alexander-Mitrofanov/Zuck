"""Exercise the production experimental beam search on a deterministic graph.

Only the SDBG storage is replaced; the search, scoring and sequence reconstruction
methods compile unchanged from the production source.
"""
import os
from pathlib import Path
import shlex
import subprocess
import tempfile


def method(source, signature):
    start = source.index(signature)
    opening = source.index("{", start)
    depth = 1
    end = opening + 1
    while depth:
        depth += (source[end] == "{") - (source[end] == "}")
        end += 1
    return source[start:end]


def main():
    root = Path(__file__).resolve().parents[1]
    source = (root / "src/protospacer_detection/phage_curator.cpp").read_text()
    parser_source = (root / "src/protospacer_detection/isolate_protospacers.cpp").read_text()
    signatures = [
        "std::string PhageCurator::_FetchFirstNode",
        "std::string PhageCurator::_FetchNodeLastBase",
        "std::vector<std::vector<uint64_t>> PhageCurator::BeamSearchPathsAvoiding",
        "std::vector<std::vector<uint64_t>> PhageCurator::GetTopPathsFromBeamPaths",
    ]
    harness = r'''
#include <algorithm>
#include <cstdint>
#include <charconv>
#include <fstream>
#include <filesystem>
#include <map>
#include <optional>
#include <sstream>
#include <deque>
#include <functional>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>
struct SDBG {
    static constexpr uint64_t kNullID = UINT64_MAX;
    std::vector<std::vector<uint64_t>> edges{5461};
    std::vector<double> multiplicities = std::vector<double>(5461, 10);
    size_t size() const { return edges.size(); }
    bool IsValidEdge(uint64_t n) const { return n < size(); }
    int EdgeOutdegree(uint64_t n) const { return edges.at(n).size(); }
    bool EdgeOutdegreeZero(uint64_t n) const { return edges.at(n).empty(); }
    int OutgoingEdges(uint64_t n, uint64_t* out) const {
        std::copy(edges.at(n).begin(), edges.at(n).end(), out);
        return edges.at(n).size();
    }
    uint64_t NextSimplePathEdge(uint64_t) const { return kNullID; }
    double EdgeMultiplicity(uint64_t n) const { return multiplicities.at(n); }
    size_t k() const { return 3; }
    void GetLabel(uint64_t, uint8_t* out) const { out[0]=1; out[1]=2; out[2]=3; }
};
char safe_base_from_code(uint8_t code) { return code <= 4 ? "NACGT"[code] : 'N'; }
class PhageCurator {
public:
    SDBG& sdbg;
    std::string _FetchFirstNode(size_t);
    std::string _FetchNodeLastBase(size_t);
    std::vector<std::vector<uint64_t>> BeamSearchPathsAvoiding(uint64_t,int,int,
        const std::set<uint64_t>&,int,double,double,
        std::function<void(const std::vector<uint64_t>&)> = nullptr);
    std::vector<std::vector<uint64_t>> GetTopPathsFromBeamPaths(
        const std::vector<std::vector<uint64_t>>&,double,double,size_t);
};
using PathGroups = std::map<uint64_t, std::map<uint64_t, std::vector<std::vector<uint64_t>>>>;
class IsolateProtospacers {
public:
    SDBG& graph;
    std::map<uint64_t,uint64_t> cycleToGroup;
    std::map<uint64_t,std::set<uint64_t>> ReadCycles(const std::string&);
    PathGroups ReadPathsFromFile(const std::string&);
};
__PARSERS__
void require(bool ok) { if (!ok) throw std::runtime_error("experimental path regression"); }
__METHODS__
int main(int argc, char** argv) {
    SDBG graph;
    for (uint64_t parent=0; parent<1365; ++parent)
        for (uint64_t child=4*parent+1; child<=4*parent+4; ++child)
            graph.edges[parent].push_back(child);
    PhageCurator search{graph};
    require(search._FetchFirstNode(0) == "ACG");
    require(search._FetchNodeLastBase(0) == "G");
    bool bad_label = false;
    try { search._FetchFirstNode(graph.size()); } catch (const std::out_of_range&) { bad_label = true; }
    require(bad_label);
    auto paths = search.BeamSearchPathsAvoiding(0,6,6,{},5000,2,20);
    require(paths.size() == 4096);
    std::set<std::vector<uint64_t>> unique(paths.begin(), paths.end());
    require(unique.size() == paths.size());
    for (const auto& path : paths) {
        require(path.size() == 7 && path.front() == 0);
        for (size_t i=1; i<path.size(); ++i)
            require((path[i]-1)/4 == path[i-1]);
    }
    paths = search.BeamSearchPathsAvoiding(0,6,6,{1},5000,2,20);
    require(paths.size() == 3072);
    size_t callbacks = 0;
    require(search.BeamSearchPathsAvoiding(0,2,2,{},100,2,20,
        [&](const auto&) { ++callbacks; }).empty());
    require(callbacks == 16);
    bool rejected = false;
    try { search.BeamSearchPathsAvoiding(0,3,2,{},5,2,20); }
    catch (const std::invalid_argument&) { rejected = true; }
    require(rejected);
    graph.multiplicities[1] = 5;
    graph.multiplicities[2] = 10;
    graph.multiplicities[3] = 20;
    auto ranked = search.GetTopPathsFromBeamPaths({{1},{3},{2}},9,11,1);
    require(ranked == std::vector<std::vector<uint64_t>>{{2}});
    IsolateProtospacers reader{graph};
    const auto fixture = (std::filesystem::path(argv[1]) / "paths.txt").string();
    { std::ofstream file(fixture); file << "Group 8:\n  Cycle 3:\n1    [4 5 6]\n"; }
    require(reader.ReadPathsFromFile(fixture).at(8).at(3).at(0) == std::vector<uint64_t>({4,5,6}));
    for (const auto* bad : {"[4 5]\n", "Group 8:\n[4]\n", "Group 8oops:\n",
             "Group 8:\nCycle 3:\n[4 bad]\n", "Group 8:\nCycle 3:\n[-1]\n",
             "Group 8:\nCycle 3:\n[6000]\n"}) {
        { std::ofstream file(fixture); file << bad; }
        bool rejected = false;
        try { reader.ReadPathsFromFile(fixture); } catch (const std::runtime_error&) { rejected = true; }
        require(rejected);
    }
    { std::ofstream file(fixture); file << "\n3 4 5\n3 5 6\n"; }
    require(reader.ReadCycles(fixture).at(3) == std::set<uint64_t>({4,5,6}));
    for (const auto* bad : {"-1 4", "3 6000", "3 4junk", "6000 3"}) {
        { std::ofstream file(fixture); file << bad; }
        bool rejected = false;
        try { reader.ReadCycles(fixture); } catch (const std::runtime_error&) { rejected = true; }
        require(rejected);
    }
}
'''.replace("__METHODS__", "\n\n".join(method(source, s) for s in signatures))
    harness = harness.replace("__PARSERS__", "\n\n".join(method(parser_source, signature) for signature in (
        "uint64_t parse_node_id", "std::map<uint64_t, std::set<uint64_t>> IsolateProtospacers::ReadCycles",
        "std::map<uint64_t, std::map<uint64_t, std::vector<std::vector<uint64_t>>>> IsolateProtospacers::ReadPathsFromFile")))
    with tempfile.TemporaryDirectory(prefix="zuck-experimental-paths-") as temporary:
        directory = Path(temporary)
        fixture = directory / "test.cpp"
        fixture.write_text(harness)
        binary = directory / "test"
        subprocess.run([*shlex.split(os.environ.get("CXX", "c++")), "-std=c++17",
                        "-O1", "-D_GLIBCXX_ASSERTIONS", str(fixture), "-o", str(binary)],
                       check=True, capture_output=True, text=True)
        subprocess.run([str(binary), str(directory)], check=True, timeout=30)
    print("Experimental beam branching, bounds, ranking and path bases passed")


if __name__ == "__main__":
    main()
