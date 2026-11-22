#include "crispr_array/motif_cliff_boundary.h"
#include "crispr_array/element_confidence.h"
#include <chrono>
#include <filesystem>
#include <stdexcept>

namespace mcb = Zuck::motif_cliff_boundary;

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
        ("zuck-motif-cliff-" + std::to_string(
            std::chrono::steady_clock::now().time_since_epoch().count()));
    Temporary() { std::filesystem::create_directory(path); }
    ~Temporary() { std::error_code ignored; std::filesystem::remove_all(path, ignored); }
};

std::string core(const mcb::Pair& pair) {
    return Zuck::element_confidence::circular_identity(pair.first + pair.second);
}

int main() {
    const std::string repeat = "GTTCACTGCCGTACAGGCAGCTA";
    require(!mcb::recut({repeat, std::string(55, 'A')}, 0, 5).first.empty() &&
            mcb::recut({repeat, std::string(55, 'A')}, 0, 5).second.size() == 50,
            "The registered 50-base motif spacer ceiling must remain inclusive");
    require(mcb::recut({repeat, std::string(56, 'A')}, 0, 5).first.empty(),
            "A motif recut yielding a 51-base spacer must be length-ineligible");
    // Two five-base terminal variants cover every distinct spacer. At six
    // bases each preceding spacer base splits support, creating the cliff.
    std::vector<Array> arrays{{repeat, {
        {repeat, std::string(26, 'A') + "GAAAAA"},
        {repeat, std::string(26, 'C') + "TAAAAA"},
        {repeat, std::string(26, 'G') + "GCCCCC"},
        {repeat, std::string(26, 'T') + "TCCCCC"},
    }}};
    const auto originals = arrays[0].entries;
    Temporary temporary;
    const auto result = mcb::apply(arrays, {}, Limits{}, temporary.path);
    require(result.counts.at("motif_cliff_boundary_selected_haplotypes") == 1 &&
            result.counts.at("motif_cliff_boundary_replaced_rows") == 4,
            "A bounded two-variant terminal cliff must recut its exact oriented haplotype");
    for (std::size_t index = 0; index < arrays[0].entries.size(); ++index) {
        require(arrays[0].entries[index].first.size() == repeat.size() + 5 &&
                arrays[0].entries[index].second.size() == 27,
                "The five-base motif must move from spacer to repeat");
        require(core(arrays[0].entries[index]) == core(originals[index]),
                "A recut must preserve the circular repeat-spacer content");
    }
    require(std::filesystem::file_size(temporary.path / "motif_cliff_boundary.tsv") > 100,
            "Every motif recut must have an explicit ledger");

    // Production defaults do not set optional CLI length limits. Internal
    // motif policy must still reject a supported cliff whose recut spacer is
    // above the registered 50-base ceiling.
    arrays = {{repeat, {
        {repeat, std::string(50, 'A') + "GAAAAA"},
        {repeat, std::string(50, 'C') + "TAAAAA"},
        {repeat, std::string(50, 'G') + "GCCCCC"},
        {repeat, std::string(50, 'T') + "TCCCCC"},
    }}};
    const auto too_long_original = arrays;
    const auto internally_bounded = mcb::apply(
        arrays, {}, UnboundedLimits{}, temporary.path);
    require(internally_bounded.counts.at("motif_cliff_boundary_replaced_rows") == 0 &&
            arrays[0].entries == too_long_original[0].entries,
            "Out-of-range motif rows must remain unchanged under unbounded CLI limits");
    require(internally_bounded.decisions.size() == 4,
            "Every out-of-range motif row must retain an audit decision");
    for (const auto& decision : internally_bounded.decisions)
        require(!decision.replaced && decision.status == "length_ineligible" &&
                decision.candidate.first.empty(),
                "An overlong motif recut must be explicitly length-ineligible");

    auto unchanged = arrays;
    const auto skipped = mcb::apply(unchanged, {0}, Limits{}, temporary.path);
    require(skipped.counts.at("motif_cliff_boundary_replaced_rows") == 0 &&
            unchanged[0].entries == arrays[0].entries,
            "Owners resolved by the exact graph must be skipped and remain idempotent");

    // Matching cliffs on both sides are not identifiable; preserve all rows.
    arrays = {{repeat, {
        {repeat, "AAAAA" + std::string(26, 'A') + "GAAAAA"},
        {repeat, "AAAAA" + std::string(26, 'C') + "TAAAAA"},
        {repeat, "CCCCC" + std::string(26, 'G') + "GCCCCC"},
        {repeat, "CCCCC" + std::string(26, 'T') + "TCCCCC"},
    }}};
    const auto ambiguous_original = arrays;
    const auto ambiguous = mcb::apply(arrays, {}, Limits{}, temporary.path);
    require(ambiguous.counts.at("motif_cliff_boundary_ambiguous_both_sides") == 1 &&
            ambiguous.counts.at("motif_cliff_boundary_replaced_rows") == 0 &&
            arrays[0].entries == ambiguous_original[0].entries,
            "A two-sided conserved-terminal pattern must force abstention");

    mcb::Config config;
    require(!mcb::find_cliff({"AAAA", "AAAA", "CCCC"}, false, config).present(),
            "Motifs shorter than five bases cannot create a cut");
    config.minimum_rows = 5;
    require(!mcb::find_cliff({"GAAAAA", "TAAAAA", "GCCCCC", "TCCCCC"},
                             true, config).present(),
            "Fewer than three configured distinct spacers cannot select a haplotype");
}
