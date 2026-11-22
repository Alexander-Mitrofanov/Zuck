#include <filesystem>
#include <fstream>
#include <random>
#include <set>
#include <stdexcept>
#include <string>
#include "crispr_array/cross_variant_scaffolds.h"

namespace cv = Zuck::cross_variant_scaffolds;
namespace dp = Zuck::directed_haplotype_paths;
struct Array { std::string consensus; std::vector<std::pair<std::string, std::string>> entries; };

static void require(bool condition, const char* reason) {
    if (!condition) throw std::runtime_error(reason);
}

int main() {
    const std::vector<std::string> spacers{
        "AAAACCCCGGGGTTTTAAAACCCCGGGGTTTA",
        "AAAACCCCGGGGTTTTAAAACCCCGGGGTTTC",
        "AAAACCCCGGGGTTTTAAAACCCCGGGGTTTG",
        "AAAACCCCGGGGTTTTAAAACCCCGGGGTTCA",
        "AAAACCCCGGGGTTTTAAAACCCCGGGGTTCC",
        "AAAACCCCGGGGTTTTAAAACCCCGGGGTTCT",
        "AAAACCCCGGGGTTTTAAAACCCCGGGGTTCG",
    };
    std::vector<Array> arrays{{"AAAA", {}}, {"TTTT", {}}};
    for (std::size_t i = 0; i < 6; ++i)
        arrays[0].entries.push_back({"AAAA", spacers[i]});
    arrays[1].entries.push_back({"TTTT", spacers[6]});
    const std::vector<std::size_t> linked_row_counts{6, 1};

    dp::Result directed;
    directed.groups = {
        {0, "AAAA", "path", {1, 2, 6}, {1, 2, 3}},
        {0, "CCCC", "singleton", {3}, {1}},
        {0, "GGGG", "singleton", {4}, {1}},
        {1, "TTTT", "singleton", {1}, {1}},
    };
    Zuck::read_linked::Result linked;
    Zuck::read_linked::Group bridge;
    bridge.parent = 0;
    bridge.entry_indices = {2, 3, 5}; // Row 5 is a valid unselected transit node.
    bridge.edges.push_back({spacers[1], spacers[4], 3, 3});
    bridge.edges.push_back({spacers[4], spacers[2], 3, 3});
    linked.scaffold_groups.push_back(bridge);
    const std::set<cv::RowKey> selected{
        {0, 1}, {0, 2}, {0, 3}, {0, 4}, {0, 6}, {1, 1}};

    const auto result = cv::infer(
        directed, linked, arrays, linked_row_counts, selected, cv::Config{1});
    require(result.groups.size() == 3,
            "cross-variant bridge did not merge exactly two base groups");
    require(result.counts.at("cross_variant_scaffold_component_joins") == 1 &&
            result.counts.at("cross_variant_scaffold_cross_haplotype_joins") == 1 &&
            result.counts.at("cross_variant_scaffold_memberships") == selected.size(),
            "cross-variant scaffold counts differ");
    const auto merged = std::find_if(result.groups.begin(), result.groups.end(),
        [](const cv::Group& group) { return group.kind == "cross_variant_scaffold"; });
    require(merged != result.groups.end() && merged->parent == 0 &&
            merged->base_groups == std::vector<std::size_t>({0, 1}) &&
            merged->entries == std::vector<std::size_t>({1, 2, 3, 6}) &&
            merged->base_for_entry.at(1) == 0 && merged->base_for_entry.at(3) == 1,
            "merged group lost row or base provenance");

    const auto gated = cv::infer(
        directed, linked, arrays, linked_row_counts, selected);
    require(gated.groups.size() == 2 &&
            gated.counts.at("cross_variant_scaffold_groups_before_parent_gate") == 3 &&
            gated.counts.at("cross_variant_scaffold_parent_gate_fallback_parents") == 2 &&
            gated.counts.at("cross_variant_scaffold_parent_gate_fallback_joins") == 1 &&
            gated.groups[0].kind == "content_family_fallback" &&
            gated.groups[0].base_groups == std::vector<std::size_t>({0, 1, 2}),
            "small parent partitions did not fall back to complete content families");

    // Evidence confined to one base group cannot manufacture a merge.
    auto internal = linked;
    internal.scaffold_groups[0].entry_indices = {1, 2, 6};
    internal.scaffold_groups[0].edges = {
        {spacers[0], spacers[1], 3, 3}, {spacers[1], spacers[5], 3, 3}};
    const auto unchanged = cv::infer(
        directed, internal, arrays, linked_row_counts, selected, cv::Config{1});
    require(unchanged.groups.size() == directed.groups.size() &&
            unchanged.counts.at("cross_variant_scaffold_component_joins") == 0,
            "within-component evidence changed the partition");

    // The selected-row equality is a fail-closed contract.
    auto incomplete = selected;
    incomplete.erase({0, 4});
    bool rejected = false;
    try { (void)cv::infer(
        directed, linked, arrays, linked_row_counts, incomplete, cv::Config{1}); }
    catch (const std::invalid_argument&) { rejected = true; }
    require(rejected, "incomplete base partition was accepted");

    auto must_reject = [&](Zuck::read_linked::Result evidence, const char* reason) {
        bool failed = false;
        try { (void)cv::infer(
            directed, evidence, arrays, linked_row_counts, selected, cv::Config{1}); }
        catch (const std::invalid_argument&) { failed = true; }
        require(failed, reason);
    };
    auto weak = linked;
    weak.scaffold_groups[0].edges[0].sequence_units = 1;
    must_reject(weak, "below-2/2 scaffold edge was accepted");
    auto missing_endpoint = linked;
    missing_endpoint.scaffold_groups[0].edges[0].spacer_a = spacers[0];
    must_reject(missing_endpoint, "nonmember scaffold endpoint was accepted");
    auto duplicate_row = linked;
    duplicate_row.scaffold_groups[0].entry_indices.push_back(2);
    must_reject(duplicate_row, "duplicate scaffold row was accepted");
    auto cycle = linked;
    cycle.scaffold_groups[0].edges.push_back({spacers[1], spacers[2], 3, 3});
    must_reject(cycle, "cyclic scaffold was accepted");
    auto branch = linked;
    branch.scaffold_groups[0].entry_indices = {1, 2, 5, 6};
    branch.scaffold_groups[0].edges = {
        {spacers[0], spacers[1], 3, 3},
        {spacers[0], spacers[4], 3, 3},
        {spacers[0], spacers[5], 3, 3}};
    must_reject(branch, "branched scaffold was accepted");
    auto wrong_parent = linked;
    wrong_parent.scaffold_groups[0].parent = 2;
    must_reject(wrong_parent, "out-of-range scaffold parent was accepted");

    const auto root = std::filesystem::temp_directory_path() /
        ("zuck-cross-variant-test-" + std::to_string(std::random_device{}()));
    std::filesystem::create_directory(root);
    cv::write_table(root, result);
    std::ifstream table(root / "cross_variant_scaffolds.tsv");
    std::string text((std::istreambuf_iterator<char>(table)), {});
    require(text.find("CrossVariant_1") != std::string::npos &&
            text.find("DirectedPath_1") != std::string::npos &&
            text.find("DirectedPath_2") != std::string::npos,
            "cross-variant provenance table differs");
    std::filesystem::remove_all(root);
}
