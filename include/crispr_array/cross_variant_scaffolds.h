#ifndef ZUCK_CROSS_VARIANT_SCAFFOLDS_H
#define ZUCK_CROSS_VARIANT_SCAFFOLDS_H

#include <algorithm>
#include <cstddef>
#include <filesystem>
#include <fstream>
#include <map>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>
#include "crispr_array/directed_haplotype_paths.h"
#include "crispr_array/read_linked_groups.h"
#include "io/result_bundle.h"

namespace Zuck { namespace cross_variant_scaffolds {

using RowKey = directed_haplotype_paths::RowKey;

struct Group {
    std::size_t parent = 0;
    std::string kind;
    std::vector<std::size_t> base_groups; // Zero-based DirectedPath indexes.
    std::vector<std::size_t> entries;
    std::map<std::size_t, std::size_t> base_for_entry;
};

struct Result {
    std::vector<Group> groups;
    std::map<std::string, std::size_t> counts;
};

struct Config {
    // Preserve a read-supported partition only when it has enough independent
    // groups to distinguish a genuinely multi-locus repeat family from a
    // two-way or three-way allelic split.  Smaller partitions fall back to the
    // complete content family without changing any selected sequence call.
    std::size_t minimum_groups_to_preserve = 4;
};

class Components {
    std::vector<std::size_t> parent_;
public:
    explicit Components(std::size_t size) : parent_(size) {
        for (std::size_t i = 0; i < size; ++i) parent_[i] = i;
    }
    std::size_t find(std::size_t value) {
        auto root = value;
        while (parent_[root] != root) root = parent_[root];
        while (parent_[value] != value) {
            const auto next = parent_[value];
            parent_[value] = root;
            value = next;
        }
        return root;
    }
    bool join(std::size_t left, std::size_t right) {
        left = find(left); right = find(right);
        if (left == right) return false;
        if (right < left) std::swap(left, right);
        parent_[right] = left;
        return true;
    }
};

// Preserve every conservative directed component, then use the already
// validated native read-linked tree components only as cross-repeat-variant
// scaffolds.  No read-link component creates or deletes a selected row.
template<class Arrays>
Result infer(const directed_haplotype_paths::Result& directed,
             const read_linked::Result& linked, const Arrays& arrays,
             const std::vector<std::size_t>& linked_row_counts,
             const std::set<RowKey>& selected_original_rows,
             const Config& config = {}) {
    if (config.minimum_groups_to_preserve == 0)
        throw std::invalid_argument("Cross-variant parent gate cannot be zero");
    Result result;
    result.counts = {
        {"cross_variant_scaffold_base_groups", directed.groups.size()},
        {"cross_variant_scaffold_validated_read_link_groups", 0},
        {"cross_variant_scaffold_evidence_groups", 0},
        {"cross_variant_scaffold_component_joins", 0},
        {"cross_variant_scaffold_cross_haplotype_joins", 0},
        {"cross_variant_scaffold_groups_before_parent_gate", 0},
        {"cross_variant_scaffold_parent_gate_threshold", config.minimum_groups_to_preserve},
        {"cross_variant_scaffold_parent_gate_preserved_parents", 0},
        {"cross_variant_scaffold_parent_gate_fallback_parents", 0},
        {"cross_variant_scaffold_parent_gate_fallback_joins", 0},
        {"cross_variant_scaffold_groups", 0},
        {"cross_variant_scaffold_memberships", 0},
    };

    std::map<RowKey, std::size_t> base_for_row;
    for (std::size_t group_index = 0; group_index < directed.groups.size(); ++group_index) {
        const auto& group = directed.groups[group_index];
        if (group.entries.empty())
            throw std::invalid_argument("Directed base group is empty");
        for (auto entry : group.entries)
            if (!base_for_row.emplace(RowKey{group.parent, entry}, group_index).second)
                throw std::invalid_argument("Directed base partition overlaps a row");
    }
    std::set<RowKey> directed_rows;
    for (const auto& item : base_for_row) directed_rows.insert(item.first);
    if (directed_rows != selected_original_rows)
        throw std::invalid_argument("Directed base partition differs from selected rows");

    if (linked_row_counts.size() > arrays.size())
        throw std::invalid_argument("Read-link row-count inventory exceeds arrays");
    std::set<std::pair<std::size_t, std::string>> component_nodes;
    for (const auto& evidence : linked.scaffold_groups) {
        if (evidence.parent >= linked_row_counts.size() || evidence.edges.empty())
            throw std::invalid_argument("Read-link scaffold has no valid parent or edges");
        const auto row_count = linked_row_counts[evidence.parent];
        if (row_count > arrays[evidence.parent].entries.size())
            throw std::invalid_argument("Read-link row-count inventory exceeds parent rows");
        std::set<std::size_t> membership;
        std::set<std::string> nodes;
        for (auto entry : evidence.entry_indices) {
            if (!entry || entry > row_count || !membership.insert(entry).second)
                throw std::invalid_argument("Read-link scaffold has an invalid or duplicate row");
            const auto spacer = read_linked::canonical(
                arrays[evidence.parent].entries[entry - 1].second);
            if (spacer.empty())
                throw std::invalid_argument("Read-link scaffold row has non-DNA spacer");
            nodes.insert(spacer);
        }
        if (nodes.size() < 2)
            throw std::invalid_argument("Read-link scaffold has fewer than two nodes");
        std::set<std::size_t> expected_membership;
        for (std::size_t index = 0; index < row_count; ++index) {
            const auto spacer = read_linked::canonical(
                arrays[evidence.parent].entries[index].second);
            if (nodes.count(spacer)) expected_membership.insert(index + 1);
        }
        if (membership != expected_membership)
            throw std::invalid_argument("Read-link scaffold does not contain complete row fibers");

        std::map<std::string, std::set<std::string>> graph;
        std::set<std::pair<std::string, std::string>> edge_keys;
        std::set<std::string> endpoints;
        for (const auto& edge : evidence.edges) {
            if (edge.sequence_units < 2 || edge.start_offsets < 2)
                throw std::invalid_argument("Read-link scaffold edge is below 2/2 support");
            auto left = read_linked::canonical(edge.spacer_a);
            auto right = read_linked::canonical(edge.spacer_b);
            if (left.empty() || right.empty() || left == right ||
                !nodes.count(left) || !nodes.count(right))
                throw std::invalid_argument("Read-link scaffold edge endpoint differs from members");
            if (right < left) std::swap(left, right);
            if (!edge_keys.emplace(left, right).second)
                throw std::invalid_argument("Read-link scaffold duplicates an edge");
            graph[left].insert(right); graph[right].insert(left);
            endpoints.insert(left); endpoints.insert(right);
        }
        if (endpoints != nodes || edge_keys.size() + 1 != nodes.size())
            throw std::invalid_argument("Read-link scaffold is not a spanning tree");
        for (const auto& item : graph)
            if (item.second.size() > 2)
                throw std::invalid_argument("Read-link scaffold is branched");
        std::set<std::string> visited;
        std::vector<std::string> pending{*nodes.begin()};
        while (!pending.empty()) {
            const auto node = pending.back(); pending.pop_back();
            if (!visited.insert(node).second) continue;
            const auto found = graph.find(node);
            if (found != graph.end())
                pending.insert(pending.end(), found->second.begin(), found->second.end());
        }
        if (visited != nodes)
            throw std::invalid_argument("Read-link scaffold graph is disconnected");
        for (const auto& node : nodes)
            if (!component_nodes.insert({evidence.parent, node}).second)
                throw std::invalid_argument(
                    "Read-link scaffold node occurs in multiple components");
        ++result.counts["cross_variant_scaffold_validated_read_link_groups"];
    }

    Components components(directed.groups.size());
    for (const auto& evidence : linked.scaffold_groups) {
        std::set<std::size_t> bases;
        for (auto entry : evidence.entry_indices) {
            const auto found = base_for_row.find({evidence.parent, entry});
            if (found != base_for_row.end()) bases.insert(found->second);
        }
        if (bases.size() < 2) continue;
        ++result.counts["cross_variant_scaffold_evidence_groups"];
        const auto first = *bases.begin();
        for (auto base : bases) {
            if (directed.groups[base].parent != evidence.parent)
                throw std::logic_error("Read-link scaffold crosses parent families");
            if (!components.join(first, base)) continue;
            ++result.counts["cross_variant_scaffold_component_joins"];
            if (directed.groups[first].haplotype != directed.groups[base].haplotype)
                ++result.counts["cross_variant_scaffold_cross_haplotype_joins"];
        }
    }

    std::map<std::size_t, std::vector<std::size_t>> bases_by_parent;
    std::map<std::size_t, std::set<std::size_t>> roots_by_parent;
    for (std::size_t base = 0; base < directed.groups.size(); ++base) {
        const auto parent = directed.groups[base].parent;
        bases_by_parent[parent].push_back(base);
        roots_by_parent[parent].insert(components.find(base));
    }
    std::set<std::size_t> fallback_parents;
    for (const auto& item : bases_by_parent) {
        const auto parent = item.first;
        const auto& bases = item.second;
        const auto group_count = roots_by_parent.at(parent).size();
        result.counts["cross_variant_scaffold_groups_before_parent_gate"] +=
            group_count;
        if (group_count >= config.minimum_groups_to_preserve) {
            ++result.counts[
                "cross_variant_scaffold_parent_gate_preserved_parents"];
            continue;
        }
        fallback_parents.insert(parent);
        ++result.counts["cross_variant_scaffold_parent_gate_fallback_parents"];
        const auto first = bases.front();
        for (auto base : bases) {
            if (components.join(first, base))
                ++result.counts[
                    "cross_variant_scaffold_parent_gate_fallback_joins"];
        }
    }

    std::map<std::size_t, std::vector<std::size_t>> merged;
    for (std::size_t i = 0; i < directed.groups.size(); ++i)
        merged[components.find(i)].push_back(i);
    for (auto& item : merged) {
        auto& bases = item.second;
        std::sort(bases.begin(), bases.end());
        Group group;
        group.parent = directed.groups[bases.front()].parent;
        group.kind = fallback_parents.count(group.parent) ? "content_family_fallback" :
            bases.size() > 1 ? "cross_variant_scaffold" :
            directed.groups[bases.front()].kind == "path"
                ? "directed_path" : "residual_singleton";
        group.base_groups = bases;
        for (auto base : bases) {
            const auto& source = directed.groups[base];
            if (source.parent != group.parent)
                throw std::logic_error("Cross-variant scaffold crosses parent families");
            for (auto entry : source.entries) {
                group.entries.push_back(entry);
                if (!group.base_for_entry.emplace(entry, base).second)
                    throw std::logic_error("Cross-variant scaffold duplicates an entry");
            }
        }
        std::sort(group.entries.begin(), group.entries.end());
        result.counts["cross_variant_scaffold_memberships"] += group.entries.size();
        result.groups.push_back(std::move(group));
    }
    std::sort(result.groups.begin(), result.groups.end(), [](const Group& left, const Group& right) {
        return std::make_pair(left.parent, left.base_groups.front()) <
               std::make_pair(right.parent, right.base_groups.front());
    });

    std::set<RowKey> final_rows;
    for (const auto& group : result.groups) for (auto entry : group.entries)
        if (!final_rows.insert({group.parent, entry}).second)
            throw std::logic_error("Cross-variant scaffold partition overlaps a row");
    if (final_rows != selected_original_rows)
        throw std::logic_error("Cross-variant scaffold partition is incomplete");
    result.counts["cross_variant_scaffold_groups"] = result.groups.size();
    return result;
}

inline void write_table(const std::filesystem::path& directory,
                        const Result& result) {
    const auto path = directory / "cross_variant_scaffolds.tsv";
    std::ofstream out(path);
    if (!out) throw std::runtime_error("Cannot write cross-variant scaffold table");
    out << "group_id\tparent_array_id\tgroup_kind\tbase_group_id\tentry_index\n";
    for (std::size_t i = 0; i < result.groups.size(); ++i) {
        const auto& group = result.groups[i];
        for (auto entry : group.entries) {
            const auto owner = group.base_for_entry.at(entry);
            out << "CrossVariant_" << i + 1 << "\tArray_" << group.parent + 1
                << '\t' << group.kind << "\tDirectedPath_" << owner + 1
                << '\t' << entry << '\n';
        }
    }
    result_bundle::close_checked(out, path);
}

}} // namespace Zuck::cross_variant_scaffolds

#endif
