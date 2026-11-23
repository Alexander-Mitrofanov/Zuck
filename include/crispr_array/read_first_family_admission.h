#ifndef ZUCK_READ_FIRST_FAMILY_ADMISSION_H
#define ZUCK_READ_FIRST_FAMILY_ADMISSION_H

#include "core/settings.h"
#include "crispr_array/read_first.h"
#include "crispr_array/supported_calls.h"

namespace Zuck { namespace read_first {

// Compare minor alleles only within a repeat family whose spacer diversity
// and independent fragment support are already well represented by the graph.
class FamilyAdmission {
    using Seed = std::pair<std::size_t, std::string>;
    using NeighbourIndex = std::map<Seed, std::vector<std::string>>;
    std::map<std::string, NeighbourIndex> represented_;

    // Two substitutions leave at least one of three disjoint blocks intact.
    static Seed seed(const std::string& s, std::size_t block) {
        const auto begin = block * s.size() / 3;
        const auto end = (block + 1) * s.size() / 3;
        return {3 * s.size() + block, s.substr(begin, end - begin)};
    }

public:
    void index(const supported_calls::Result& graph, const std::vector<Call>& strict,
               const Settings& settings) {
        std::map<std::string, std::vector<const Call*>> observed;
        for (const auto& call : strict) {
            if (settings.dna_sequence_settings.pair_allowed(call.repeat, call.spacer) &&
                supported_calls::length_allowed(graph.profile, call.spacer.size()))
                observed[call.repeat].push_back(&call);
        }
        std::map<std::string, std::set<std::string>> graph_spacers;
        for (const auto& row : graph.rows) {
            if (!row.selected) continue;
            const auto reverse = read_linked::reverse_complement(row.repeat_variant);
            const auto& repeat = std::min(row.repeat_variant, reverse);
            if (!observed.count(repeat)) continue;
            // Read-first calls follow their canonical repeat's orientation.
            const auto spacer = repeat == reverse && repeat != row.repeat_variant ?
                read_linked::reverse_complement(row.spacer) : row.spacer;
            graph_spacers[repeat].insert(reverse == row.repeat_variant ? row.canonical_spacer : spacer);
        }
        for (const auto& [repeat, calls] : observed) {
            const auto graph_family = graph_spacers.find(repeat);
            if (graph_family == graph_spacers.end()) continue;
            std::vector<const Call*> corroborated;
            std::size_t shared = 0;
            for (const auto* call : calls) {
                if (graph_family->second.count(call->spacer)) {
                    ++shared;
                    if (call->templates >= 3) corroborated.push_back(call);
                }
            }
            // Three shared interiors, graph support for at least half of the
            // read spacers and a lower median of at least three templates.
            if (shared < 3 || 2 * shared < calls.size() || 2 * corroborated.size() <= shared) continue;
            auto& index = represented_[repeat];
            for (const auto* call : corroborated)
                for (std::size_t block = 0; block != 3; ++block)
                    index[seed(call->spacer, block)].push_back(call->spacer);
        }
    }

    std::size_t family_count() const { return represented_.size(); }

    const char* decision(const Call& call) const {
        if (call.templates != 1) return "selected";
        const auto family = represented_.find(call.repeat);
        if (family == represented_.end()) return "selected";
        const std::string* parent = nullptr;
        std::uint64_t differences = 0;
        for (std::size_t block = 0; block != 3; ++block) {
            const auto hits = family->second.find(seed(call.spacer, block));
            if (hits == family->second.end()) continue;
            for (const auto& other : hits->second) {
                if (parent && *parent == other) continue;
                unsigned mismatch_count = 0;
                std::uint64_t mask = 0;
                for (std::size_t i = 0; i < other.size(); ++i) {
                    if (other[i] == call.spacer[i]) continue;
                    if (++mismatch_count > 2) break;
                    mask |= std::uint64_t{1} << i;
                }
                if (!mismatch_count || mismatch_count > 2) continue;
                if (parent) return "selected_dense_ambiguous_parent";
                parent = &other;
                differences = mask;
            }
        }
        if (!parent) return "selected";
        if ((call.whole_q30_bases & differences) == differences)
            return "selected_dense_whole_q30";
        return "dense_family_singleton_allele";
    }
};

}} // namespace Zuck::read_first
#endif
