#include "crispr_array/local_boundary_graph.h"
#include <chrono>
#include <filesystem>
#include <fstream>
#include <stdexcept>

namespace lbg = Zuck::local_boundary_graph;

struct Array {
    std::string consensus;
    std::vector<std::pair<std::string, std::string>> entries;
};
struct Limits {
    bool pair_allowed(const std::string& repeat, const std::string& spacer) const {
        return repeat.size() >= 20 && repeat.size() <= 55 &&
               spacer.size() >= 26 && spacer.size() <= 50;
    }
};
struct UnboundedLimits {
    bool pair_allowed(const std::string&, const std::string&) const { return true; }
};
void require(bool value, const char* message) {
    if (!value) throw std::runtime_error(message);
}
struct Temporary {
    std::filesystem::path path = std::filesystem::temp_directory_path() /
        ("zuck-local-boundary-" + std::to_string(
            std::chrono::steady_clock::now().time_since_epoch().count()));
    Temporary() { std::filesystem::create_directory(path); }
    ~Temporary() { std::error_code ignored; std::filesystem::remove_all(path, ignored); }
};

std::string spacer(char base, char tail) {
    return std::string(25, base) + tail;
}

void add_read(std::ofstream& out, int& id, const std::string& left_repeat,
              const std::string& spacer, const std::string& right_repeat,
              const std::string& prefix) {
    const auto sequence = prefix + left_repeat + spacer + right_repeat + "TGC";
    out << "@r" << ++id << "\n" << sequence << "\n+\n"
        << std::string(sequence.size(), 'I') << "\n";
}

std::vector<Array> fixture(const std::string& seed, const std::string& a,
                           const std::string& c) {
    return {{seed, {{seed, a + spacer('G', 'A')}, {seed, a + spacer('T', 'C')},
                    {seed, c + spacer('A', 'G')}, {seed, c + spacer('C', 'T')},
                    {seed, a + spacer('C', 'A')}}}};
}

int main() {
    const std::string seed = "GTTCACTGCCGTACAGGCAGCTA";
    const std::string a = "AAAAA", c = "CCCCC";
    const lbg::Config registered;
    require(lbg::family_candidate_allowed(
                {std::string(20, 'A'), std::string(50, 'C')},
                registered, UnboundedLimits{}),
            "The registered 50-base family spacer ceiling must remain inclusive");
    require(!lbg::family_candidate_allowed(
                {std::string(20, 'A'), std::string(51, 'C')},
                registered, UnboundedLimits{}),
            "A 51-base family spacer must be rejected without relying on CLI limits");
    require(lbg::family_candidate_allowed(
                {std::string(20, 'A'), std::string(20, 'C')},
                registered, UnboundedLimits{}) &&
            !lbg::family_candidate_allowed(
                {std::string(20, 'A'), std::string(19, 'C')},
                registered, UnboundedLimits{}),
            "The intentional 20-base family recut minimum must remain inclusive");
    Temporary temporary;
    const auto reads = temporary.path / "reads.fastq";
    {
        std::ofstream out(reads);
        int id = 0;
        for (const auto& item : std::vector<std::pair<std::string, std::string>>{
                {a, spacer('G', 'A')}, {a, spacer('T', 'C')},
                {c, spacer('A', 'G')}, {c, spacer('C', 'T')}}) {
            add_read(out, id, seed + item.first, item.second, seed + item.first, "");
            add_read(out, id, seed + item.first, item.second, seed + item.first, "G");
        }
    }
    auto arrays = fixture(seed, a, c);
    const auto result = lbg::apply(arrays, {reads.string()}, Limits{}, temporary.path);
    require(result.counts.at("local_boundary_graph_selected_families") == 1,
            "A supported two-haplotype geometry must be selected");
    require(result.counts.at("local_boundary_graph_replaced_rows") == 5,
            "Every family row, including a row without individual read support, must be recut");
    for (const auto& pair : arrays[0].entries)
        require(pair.first.size() == seed.size() + 5 && pair.second.size() == 26,
                "Resolution must move the five conserved repeat bases out of the spacer");
    require(std::filesystem::file_size(temporary.path / "local_boundary_graph.tsv") > 100,
            "The boundary evidence ledger must be written");

    // Standard production length settings are unconstrained. The boundary
    // module must still enforce its own registered spacer ceiling while
    // propagating one supported family geometry to rows without read support.
    arrays = fixture(seed, a, c);
    const lbg::Pair long_original{seed, std::string(66, 'G')};
    arrays[0].entries.push_back(long_original);
    const auto internally_bounded = lbg::apply(
        arrays, {reads.string()}, UnboundedLimits{}, temporary.path);
    require(internally_bounded.counts.at("local_boundary_graph_replaced_rows") == 5,
            "An out-of-range family row must not inflate the replacement count");
    require(arrays[0].entries.back() == long_original,
            "A family propagation yielding a spacer above 50 bases must abstain");
    for (const auto& decision : internally_bounded.decisions)
        if (decision.replaced)
            require(decision.key.pair.first.size() >= 20 &&
                    decision.key.pair.first.size() <= 55 &&
                    decision.key.pair.second.size() >= 20 &&
                    decision.key.pair.second.size() <= 50,
                    "Every replaced family ledger row must satisfy registered lengths");

    // Conserved spacer terminals can create the same apparent repeat alleles.
    // A single observed transition between supported terminal types makes the
    // cut non-identifiable and must preserve every original row.
    const auto cross_reads = temporary.path / "cross.fastq";
    {
        std::ofstream out(cross_reads);
        int id = 0;
        for (const auto& item : std::vector<std::pair<std::string, std::string>>{
                {a, spacer('G', 'A')}, {a, spacer('T', 'C')},
                {c, spacer('A', 'G')}, {c, spacer('C', 'T')}}) {
            add_read(out, id, seed + item.first, item.second, seed + item.first, "");
            add_read(out, id, seed + item.first, item.second, seed + item.first, "G");
        }
        add_read(out, id, seed + a, spacer('G', 'T'), seed + c, "A");
        add_read(out, id, seed + c, spacer('T', 'A'), seed + a, "C");
    }
    arrays = fixture(seed, a, c);
    const auto ambiguous = lbg::apply(arrays, {cross_reads.string()}, Limits{}, temporary.path);
    require(ambiguous.counts.at("local_boundary_graph_cross_transition_geometries") > 0,
            "Supported cross transitions must be measured");
    require(ambiguous.counts.at("local_boundary_graph_replaced_rows") == 0,
            "A cross-haplotype transition must force boundary abstention");
    for (const auto& pair : arrays[0].entries)
        require(pair.first == seed && pair.second.size() == 31,
                "Ambiguous conserved-terminal rows must remain unchanged");
}
