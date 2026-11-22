#ifndef ZUCK_MOTIF_CLIFF_BOUNDARY_H
#define ZUCK_MOTIF_CLIFF_BOUNDARY_H

#include <algorithm>
#include <filesystem>
#include <fstream>
#include <limits>
#include <map>
#include <set>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>
#include <vector>
#include "io/result_bundle.h"

namespace Zuck { namespace motif_cliff_boundary {

using Pair = std::pair<std::string, std::string>;

struct Config {
    std::size_t minimum_rows = 3;
    std::size_t maximum_variants = 2;
    std::size_t minimum_extension = 5;
    std::size_t maximum_extension = 16;
    std::size_t minimum_spacer = 26;
    std::size_t maximum_spacer = 50;
    std::size_t maximum_repeat = 55;
    std::size_t minimum_fraction_numerator = 4;
    std::size_t minimum_fraction_denominator = 5;
    std::size_t maximum_child_numerator = 7;
    std::size_t maximum_child_denominator = 10;
};

struct Decision {
    std::size_t owner = 0, entry = 0, left = 0, right = 0;
    std::size_t rows = 0, support = 0, child_support = 0;
    Pair original, candidate;
    bool replaced = false;
    std::string status;
};

struct Result {
    std::vector<Decision> decisions;
    std::map<std::string, std::size_t> counts;
};

inline std::size_t top_support(const std::set<std::string>& spacers,
        std::size_t length, bool suffix, std::size_t maximum_variants) {
    std::map<std::string, std::size_t> counts;
    for (const auto& spacer : spacers) {
        if (spacer.size() < length) continue;
        const auto motif = suffix ? spacer.substr(spacer.size() - length)
                                  : spacer.substr(0, length);
        ++counts[motif];
    }
    std::vector<std::size_t> support;
    for (const auto& item : counts) support.push_back(item.second);
    std::sort(support.begin(), support.end(), std::greater<std::size_t>());
    std::size_t result = 0;
    for (std::size_t index = 0;
         index < std::min(maximum_variants, support.size()); ++index)
        result += support[index];
    return result;
}

struct Cliff {
    std::size_t length = 0, support = 0, child_support = 0;
    bool present() const { return length != 0; }
};

inline Cliff find_cliff(const std::set<std::string>& spacers, bool suffix,
                        const Config& config = {}) {
    Cliff best;
    if (spacers.size() < config.minimum_rows) return best;
    for (std::size_t length = config.minimum_extension;
         length <= config.maximum_extension; ++length) {
        const auto support = top_support(
            spacers, length, suffix, config.maximum_variants);
        const auto child = top_support(
            spacers, length + 1, suffix, config.maximum_variants);
        if (support * config.minimum_fraction_denominator <
                spacers.size() * config.minimum_fraction_numerator ||
            child * config.maximum_child_denominator >
                support * config.maximum_child_numerator)
            continue;
        const auto rank = std::make_tuple(support, length,
            std::numeric_limits<std::size_t>::max() - child);
        const auto old_rank = std::make_tuple(best.support, best.length,
            std::numeric_limits<std::size_t>::max() - best.child_support);
        if (!best.present() || rank > old_rank)
            best = {length, support, child};
    }
    return best;
}

inline Pair recut(const Pair& original, std::size_t left, std::size_t right,
                  const Config& config = {}) {
    if ((!left && !right) || original.second.size() < left + right + config.minimum_spacer ||
        original.first.size() + left + right > config.maximum_repeat ||
        original.second.size() - left - right > config.maximum_spacer)
        return {};
    const auto end = original.second.size() - left;
    return {(left ? original.second.substr(end) : std::string()) + original.first +
                (right ? original.second.substr(0, right) : std::string()),
            original.second.substr(right, end - right)};
}

// Resolve a terminal motif only within an exact oriented repeat haplotype.
// The policy is replacement-only, skips owners already resolved by the exact
// read graph, and abstains when both spacer ends exhibit a cliff.
template<class Arrays, class Limits>
Result apply(Arrays& arrays, const std::set<std::size_t>& skip_owners,
             const Limits& limits, const std::filesystem::path& directory,
             const Config& config = {}) {
    if (!config.minimum_rows || !config.maximum_variants ||
        !config.minimum_extension || config.maximum_extension < config.minimum_extension ||
        !config.minimum_spacer || config.maximum_spacer < config.minimum_spacer ||
        !config.maximum_repeat ||
        !config.minimum_fraction_numerator || !config.minimum_fraction_denominator ||
        config.minimum_fraction_numerator > config.minimum_fraction_denominator ||
        !config.maximum_child_numerator || !config.maximum_child_denominator ||
        config.maximum_child_numerator > config.maximum_child_denominator)
        throw std::invalid_argument("Invalid motif-cliff boundary configuration");

    Result result;
    result.counts = {{"motif_cliff_boundary_applied", 1},
        {"motif_cliff_boundary_input_haplotypes", 0},
        {"motif_cliff_boundary_skipped_exact_graph_owners", skip_owners.size()},
        {"motif_cliff_boundary_ambiguous_both_sides", 0},
        {"motif_cliff_boundary_selected_haplotypes", 0},
        {"motif_cliff_boundary_replaced_rows", 0}};

    for (std::size_t owner = 0; owner < arrays.size(); ++owner) {
        if (skip_owners.count(owner)) continue;
        std::map<std::string, std::set<std::string>> spacers_by_repeat;
        for (const auto& pair : arrays[owner].entries)
            spacers_by_repeat[pair.first].insert(pair.second);
        result.counts["motif_cliff_boundary_input_haplotypes"] +=
            spacers_by_repeat.size();
        for (const auto& group : spacers_by_repeat) {
            if (group.second.size() < config.minimum_rows) continue;
            const auto prefix = find_cliff(group.second, false, config);
            const auto suffix = find_cliff(group.second, true, config);
            if (prefix.present() && suffix.present()) {
                ++result.counts["motif_cliff_boundary_ambiguous_both_sides"];
                continue;
            }
            const auto left = suffix.length, right = prefix.length;
            if (!left && !right) continue;
            bool selected = false;
            for (std::size_t entry = 0; entry < arrays[owner].entries.size(); ++entry) {
                if (arrays[owner].entries[entry].first != group.first) continue;
                Decision decision;
                decision.owner = owner; decision.entry = entry;
                decision.left = left; decision.right = right;
                decision.rows = group.second.size();
                decision.support = suffix.present() ? suffix.support : prefix.support;
                decision.child_support = suffix.present() ? suffix.child_support : prefix.child_support;
                decision.original = arrays[owner].entries[entry];
                decision.candidate = recut(decision.original, left, right, config);
                if (decision.candidate.first.empty() ||
                    !limits.pair_allowed(decision.candidate.first, decision.candidate.second)) {
                    decision.status = "length_ineligible";
                } else {
                    arrays[owner].entries[entry] = decision.candidate;
                    decision.replaced = true;
                    decision.status = "replaced";
                    selected = true;
                    ++result.counts["motif_cliff_boundary_replaced_rows"];
                }
                result.decisions.push_back(std::move(decision));
            }
            result.counts["motif_cliff_boundary_selected_haplotypes"] += selected;
        }
    }

    const auto path = directory / "motif_cliff_boundary.tsv";
    std::ofstream out(path);
    if (!out) throw std::runtime_error("Cannot write motif-cliff boundary ledger");
    out << "owner_array_id\tentry_index\toriginal_repeat\toriginal_spacer"
           "\tcandidate_repeat\tcandidate_spacer\tleft_extension\tright_extension"
           "\tdistinct_original_spacers\ttop_variant_support\tchild_variant_support"
           "\treplaced\tstatus\n";
    for (const auto& row : result.decisions)
        out << "Array_" << row.owner + 1 << '\t' << row.entry + 1 << '\t'
            << row.original.first << '\t' << row.original.second << '\t'
            << row.candidate.first << '\t' << row.candidate.second << '\t'
            << row.left << '\t' << row.right << '\t' << row.rows << '\t'
            << row.support << '\t' << row.child_support << '\t'
            << row.replaced << '\t' << row.status << '\n';
    result_bundle::close_checked(out, path);
    return result;
}

}} // namespace Zuck::motif_cliff_boundary

#endif
