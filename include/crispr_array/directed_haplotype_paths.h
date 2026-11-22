#ifndef ZUCK_DIRECTED_HAPLOTYPE_PATHS_H
#define ZUCK_DIRECTED_HAPLOTYPE_PATHS_H

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <filesystem>
#include <fstream>
#include <map>
#include <set>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>
#include <vector>
#include "crispr_array/read_linked_groups.h"
#include "io/result_bundle.h"

namespace Zuck { namespace directed_haplotype_paths {

using RowKey = std::pair<std::size_t, std::size_t>; // family, one-based entry
using EdgeKey = std::tuple<std::size_t, std::string, std::size_t, std::size_t>;

struct Config {
    std::size_t minimum_anchor = 16;
    int minimum_repeat = 20;
    int maximum_repeat = 55;
    int minimum_quality = 10;
    double maximum_expected_errors = 1.0;
    std::size_t minimum_sequence_units = 2;
    std::size_t minimum_start_offsets = 2;
};

struct EdgeDecision {
    EdgeKey key;
    std::string spacer_from, spacer_to;
    std::size_t sequence_units = 0, start_offsets = 0;
    std::string status;
};

struct Group {
    std::size_t parent = 0;
    std::string haplotype, kind;
    std::vector<std::size_t> entries;
    std::vector<std::size_t> positions;
};

struct Result {
    std::vector<EdgeDecision> edges;
    std::vector<Group> groups;
    std::map<std::string, std::size_t> counts;
};

struct Endpoint {
    std::size_t entry = 0;
    std::string sequence;
};

struct HaplotypeIndex {
    std::map<std::string, std::vector<Endpoint>> forward_left, forward_right;
    std::map<std::string, std::vector<Endpoint>> reverse_left, reverse_right;
    std::map<std::size_t, std::string> spacers;
    std::map<std::size_t, std::vector<std::size_t>> fibers;
    std::map<std::string, std::vector<std::pair<std::size_t, std::string>>> node_rows;
};

struct PatternUse {
    std::size_t family = 0;
    std::string haplotype;
    bool reverse = false;
};

struct Support {
    std::set<std::size_t> sequence_units, start_offsets;
};

inline void add_endpoint(std::map<std::string, std::vector<Endpoint>>& index,
                         const std::string& anchor, std::size_t entry,
                         const std::string& sequence) {
    index[anchor].push_back({entry, sequence});
}

template<class Arrays>
Result infer(const Arrays& arrays, const std::vector<std::string>& files,
             const std::set<RowKey>& selected_original_rows,
             const Config& config = {}) {
    if (files.empty() || files.size() > 2)
        throw std::runtime_error("Directed haplotype paths require one or two read files");
    if (!config.minimum_anchor || config.minimum_repeat <= 0 ||
        config.maximum_repeat < config.minimum_repeat || config.minimum_quality < 0 ||
        config.minimum_quality > 93 || !std::isfinite(config.maximum_expected_errors) ||
        config.maximum_expected_errors < 0 || !config.minimum_sequence_units ||
        !config.minimum_start_offsets)
        throw std::runtime_error("Invalid directed haplotype path configuration");

    Result result;
    result.counts = {{"directed_path_input_files", files.size()},
        {"directed_path_selected_original_rows", selected_original_rows.size()},
        {"directed_path_records", 0}, {"directed_path_records_without_quality", 0},
        {"directed_path_input_units", 0}, {"directed_path_distinct_sequence_units", 0},
        {"directed_path_candidate_edges", 0}, {"directed_path_supported_edges", 0},
        {"directed_path_palindromic_haplotype_rows", 0},
        {"directed_path_reverse_conflict_edges", 0},
        {"directed_path_competing_edges", 0}, {"directed_path_cycle_edges", 0},
        {"directed_path_accepted_edges", 0}, {"directed_path_path_groups", 0},
        {"directed_path_singleton_groups", 0}, {"directed_path_groups", 0},
        {"directed_path_memberships", 0}};

    std::map<std::pair<std::size_t, std::string>, HaplotypeIndex> indexes;
    for (const auto& selected : selected_original_rows) {
        const auto family = selected.first, entry = selected.second;
        if (family >= arrays.size() || !entry || entry > arrays[family].entries.size())
            throw std::invalid_argument("Directed path selection references an absent row");
        const auto& row = arrays[family].entries[entry - 1];
        if (row.first.size() < static_cast<std::size_t>(config.minimum_repeat) ||
            row.first.size() > static_cast<std::size_t>(config.maximum_repeat) ||
            row.second.size() < config.minimum_anchor)
            throw std::invalid_argument("Directed path selected row is outside the registered geometry");
        const auto reverse_repeat = read_linked::reverse_complement(row.first);
        const auto reverse_spacer = read_linked::reverse_complement(row.second);
        if (reverse_repeat.empty() || reverse_spacer.empty())
            throw std::invalid_argument("Directed path selected row contains non-DNA sequence");
        const auto haplotype = std::min(row.first, reverse_repeat);
        const auto oriented_spacer = row.first == haplotype ? row.second : reverse_spacer;
        const auto canonical_spacer = read_linked::canonical(oriented_spacer);
        indexes[{family, haplotype}].node_rows[canonical_spacer].push_back(
            {entry, oriented_spacer});
    }
    for (auto& item : indexes) {
        auto& index = item.second;
        for (auto& node : index.node_rows) {
            std::sort(node.second.begin(), node.second.end());
            const auto representative = node.second.front().first;
            std::set<std::string> orientations;
            for (const auto& row : node.second) {
                index.fibers[representative].push_back(row.first);
                orientations.insert(row.second);
            }
            index.spacers[representative] = node.second.front().second;
            // Same canonical spacer under incompatible row orientations is
            // intrinsically unresolved: retain the full fiber as one residual.
            if (orientations.size() != 1) continue;
            const auto& oriented_spacer = *orientations.begin();
            add_endpoint(index.forward_left,
                oriented_spacer.substr(oriented_spacer.size() - config.minimum_anchor),
                representative, oriented_spacer);
            add_endpoint(index.forward_right,
                oriented_spacer.substr(0, config.minimum_anchor),
                representative, oriented_spacer);
            const auto reverse_spacer = read_linked::reverse_complement(oriented_spacer);
            add_endpoint(index.reverse_left,
                reverse_spacer.substr(reverse_spacer.size() - config.minimum_anchor),
                representative, reverse_spacer);
            add_endpoint(index.reverse_right,
                reverse_spacer.substr(0, config.minimum_anchor),
                representative, reverse_spacer);
        }
    }
    if (selected_original_rows.empty()) return result;

    std::map<std::string, std::vector<PatternUse>> use_by_pattern;
    for (const auto& item : indexes) {
        const auto family = item.first.first;
        const auto& haplotype = item.first.second;
        const auto reverse = read_linked::reverse_complement(haplotype);
        if (reverse.empty())
            throw std::invalid_argument("Directed path haplotype contains non-DNA sequence");
        if (reverse == haplotype) {
            result.counts["directed_path_palindromic_haplotype_rows"] +=
                item.second.spacers.size();
            continue;
        }
        use_by_pattern[haplotype].push_back({family, haplotype, false});
        use_by_pattern[reverse].push_back({family, haplotype, true});
    }

    read_evidence::ExactMatcher matcher;
    std::vector<std::string> patterns;
    std::vector<std::vector<PatternUse>> pattern_uses;
    for (const auto& item : use_by_pattern) {
        patterns.push_back(item.first);
        pattern_uses.push_back(item.second);
        matcher.add(item.first, patterns.size() - 1);
    }
    matcher.build();

    std::map<EdgeKey, Support> support;
    std::map<std::pair<std::string, std::string>, std::size_t> sequence_units;
    auto scan = [&](const klibpp::KSeq& read,
                    std::map<EdgeKey, std::set<std::size_t>>& links) {
        ++result.counts["directed_path_records"];
        if (read.qual.empty()) {
            ++result.counts["directed_path_records_without_quality"];
            return;
        }
        read_evidence::QualitySpans quality(read.qual, config.minimum_quality);
        matcher.scan(read.seq, [&](std::size_t pattern_id, std::size_t end) {
            const auto start = end - patterns[pattern_id].size();
            if (start < config.minimum_anchor || config.minimum_anchor > read.seq.size() - end)
                return;
            const auto left_anchor = read.seq.substr(
                start - config.minimum_anchor, config.minimum_anchor);
            const auto right_anchor = read.seq.substr(end, config.minimum_anchor);
            for (const auto& use : pattern_uses[pattern_id]) {
                const auto& index = indexes.at({use.family, use.haplotype});
                const auto& left_index = use.reverse ? index.reverse_left : index.forward_left;
                const auto& right_index = use.reverse ? index.reverse_right : index.forward_right;
                const auto left = left_index.find(left_anchor);
                const auto right = right_index.find(right_anchor);
                if (left == left_index.end() || right == right_index.end()) continue;
                std::set<std::size_t> left_entries, right_entries;
                for (const auto& endpoint : left->second) left_entries.insert(endpoint.entry);
                for (const auto& endpoint : right->second) right_entries.insert(endpoint.entry);
                if (left_entries.size() != 1 || right_entries.size() != 1 ||
                    *left_entries.begin() == *right_entries.begin()) continue;
                for (const auto& left_endpoint : left->second) {
                    for (const auto& right_endpoint : right->second) {
                        const auto left_length = std::min(start, left_endpoint.sequence.size());
                        const auto right_length = std::min(
                            read.seq.size() - end, right_endpoint.sequence.size());
                        if (read.seq.compare(start - left_length, left_length,
                                left_endpoint.sequence,
                                left_endpoint.sequence.size() - left_length,
                                left_length) != 0 ||
                            read.seq.compare(end, right_length,
                                right_endpoint.sequence, 0, right_length) != 0 ||
                            !quality.passes(start - left_length, end + right_length,
                                            config.maximum_expected_errors))
                            continue;
                        const auto from = use.reverse ? right_endpoint.entry : left_endpoint.entry;
                        const auto to = use.reverse ? left_endpoint.entry : right_endpoint.entry;
                        const auto offset = use.reverse ? read.seq.size() - end : start;
                        links[{use.family, use.haplotype, from, to}].insert(offset);
                    }
                }
            }
        });
    };

    read_linked::Reader first(files[0]);
    std::unique_ptr<read_linked::Reader> second;
    if (files.size() == 2) second = std::make_unique<read_linked::Reader>(files[1]);
    klibpp::KSeq a, b;
    while (first.next(a)) {
        if (second && (!second->next(b) || read_linked::read_name(a.name) !=
                       read_linked::read_name(b.name)))
            throw std::runtime_error("Directed path read pairs have different names or lengths");
        std::map<EdgeKey, std::set<std::size_t>> links;
        scan(a, links);
        if (second) scan(b, links);
        ++result.counts["directed_path_input_units"];
        if (links.empty()) continue;
        const auto unit = sequence_units.emplace(
            std::make_pair(a.seq, second ? b.seq : ""), sequence_units.size());
        for (const auto& item : links) {
            auto& evidence = support[item.first];
            evidence.sequence_units.insert(unit.first->second);
            evidence.start_offsets.insert(item.second.begin(), item.second.end());
        }
    }
    if (second && second->next(b))
        throw std::runtime_error("Directed path read files contain different record counts");
    result.counts["directed_path_distinct_sequence_units"] = sequence_units.size();
    result.counts["directed_path_candidate_edges"] = support.size();

    std::map<EdgeKey, std::size_t> edge_index;
    for (const auto& item : support) {
        const auto& [family, haplotype, from, to] = item.first;
        const auto& index = indexes.at({family, haplotype});
        EdgeDecision decision;
        decision.key = item.first;
        decision.spacer_from = index.spacers.at(from);
        decision.spacer_to = index.spacers.at(to);
        decision.sequence_units = item.second.sequence_units.size();
        decision.start_offsets = item.second.start_offsets.size();
        decision.status = decision.sequence_units < config.minimum_sequence_units
            ? "insufficient_sequence_units"
            : decision.start_offsets < config.minimum_start_offsets
            ? "insufficient_start_offsets" : "supported";
        edge_index[item.first] = result.edges.size();
        result.edges.push_back(std::move(decision));
        if (result.edges.back().status == "supported")
            ++result.counts["directed_path_supported_edges"];
    }

    for (auto& edge : result.edges) if (edge.status == "supported") {
        const auto& [family, haplotype, from, to] = edge.key;
        const auto reverse = support.find({family, haplotype, to, from});
        if (reverse != support.end() && !reverse->second.sequence_units.empty()) {
            edge.status = "reverse_direction_conflict";
            ++result.counts["directed_path_reverse_conflict_edges"];
        }
    }

    std::map<std::pair<std::size_t, std::string>,
             std::map<std::size_t, std::vector<std::size_t>>> outgoing, incoming;
    for (std::size_t i = 0; i < result.edges.size(); ++i) {
        const auto& edge = result.edges[i];
        if (edge.status != "supported") continue;
        const auto& [family, haplotype, from, to] = edge.key;
        outgoing[{family, haplotype}][from].push_back(i);
        incoming[{family, haplotype}][to].push_back(i);
    }
    auto rank = [&](std::size_t index) {
        const auto& edge = result.edges[index];
        return std::make_pair(std::min(edge.sequence_units, edge.start_offsets),
                              edge.sequence_units + edge.start_offsets);
    };
    auto unique_best = [&](const std::vector<std::size_t>& choices,
                           std::size_t candidate) {
        auto best = std::make_pair(std::size_t(0), std::size_t(0));
        std::size_t count = 0;
        for (auto index : choices) {
            const auto value = rank(index);
            if (value > best) { best = value; count = 1; }
            else if (value == best) ++count;
        }
        return count == 1 && rank(candidate) == best;
    };
    for (std::size_t i = 0; i < result.edges.size(); ++i) {
        auto& edge = result.edges[i];
        if (edge.status != "supported") continue;
        const auto& [family, haplotype, from, to] = edge.key;
        if (!unique_best(outgoing[{family, haplotype}][from], i) ||
            !unique_best(incoming[{family, haplotype}][to], i)) {
            edge.status = "competing_endpoint";
            ++result.counts["directed_path_competing_edges"];
        }
    }

    for (const auto& index_item : indexes) {
        const auto& identity = index_item.first;
        std::map<std::size_t, std::size_t> next;
        std::map<std::size_t, EdgeKey> key_by_from;
        for (const auto& edge : result.edges) {
            const auto& [family, haplotype, from, to] = edge.key;
            if (edge.status != "supported" ||
                std::make_pair(family, haplotype) != identity) continue;
            next[from] = to;
            key_by_from[from] = edge.key;
        }
        std::set<std::size_t> complete;
        for (const auto& vertex : index_item.second.spacers) {
            std::vector<std::size_t> chain;
            std::map<std::size_t, std::size_t> position;
            auto current = vertex.first;
            while (next.count(current) && !complete.count(current) &&
                   !position.count(current)) {
                position[current] = chain.size();
                chain.push_back(current);
                current = next.at(current);
            }
            if (position.count(current)) {
                for (std::size_t i = position[current]; i < chain.size(); ++i) {
                    auto& edge = result.edges.at(edge_index.at(key_by_from.at(chain[i])));
                    if (edge.status == "supported") {
                        edge.status = "cycle";
                        ++result.counts["directed_path_cycle_edges"];
                    }
                }
            }
            complete.insert(chain.begin(), chain.end());
        }
    }

    for (const auto& index_item : indexes) {
        const auto& identity = index_item.first;
        std::map<std::size_t, std::size_t> next;
        std::set<std::size_t> incoming_vertices;
        for (auto& edge : result.edges) {
            const auto& [family, haplotype, from, to] = edge.key;
            if (edge.status != "supported" ||
                std::make_pair(family, haplotype) != identity) continue;
            edge.status = "accepted";
            ++result.counts["directed_path_accepted_edges"];
            next[from] = to;
            incoming_vertices.insert(to);
        }
        std::set<std::size_t> assigned;
        for (const auto& item : next) {
            if (incoming_vertices.count(item.first)) continue;
            Group group{identity.first, identity.second, "path", {}};
            auto current = item.first;
            std::size_t path_position = 0;
            while (true) {
                if (!assigned.insert(current).second)
                    throw std::logic_error("Directed path vertex assigned twice");
                ++path_position;
                for (auto entry : index_item.second.fibers.at(current)) {
                    group.entries.push_back(entry);
                    group.positions.push_back(path_position);
                }
                const auto successor = next.find(current);
                if (successor == next.end()) break;
                current = successor->second;
            }
            if (path_position < 2)
                throw std::logic_error("Directed path has fewer than two vertices");
            result.groups.push_back(std::move(group));
            ++result.counts["directed_path_path_groups"];
        }
        for (const auto& vertex : index_item.second.spacers) {
            if (assigned.count(vertex.first)) continue;
            Group group{identity.first, identity.second, "singleton", {}, {}};
            for (auto entry : index_item.second.fibers.at(vertex.first)) {
                group.entries.push_back(entry);
                group.positions.push_back(1);
            }
            result.groups.push_back(std::move(group));
            ++result.counts["directed_path_singleton_groups"];
        }
    }

    std::set<RowKey> memberships;
    for (const auto& group : result.groups) for (auto entry : group.entries) {
        if (!memberships.insert({group.parent, entry}).second)
            throw std::logic_error("Directed path partition overlaps a row");
    }
    if (memberships != selected_original_rows)
        throw std::logic_error("Directed path partition does not cover selected original rows");
    result.counts["directed_path_groups"] = result.groups.size();
    result.counts["directed_path_memberships"] = memberships.size();
    return result;
}

inline void write_tables(const std::filesystem::path& directory,
                         const Result& result) {
    const auto groups_path = directory / "directed_haplotype_paths.tsv";
    const auto edges_path = directory / "directed_haplotype_edges.tsv";
    std::ofstream groups(groups_path), edges(edges_path);
    if (!groups || !edges)
        throw std::runtime_error("Cannot write directed haplotype path tables");
    groups << "group_id\tparent_array_id\trepeat_haplotype\tgroup_kind\tpath_position\tentry_index\n";
    edges << "parent_array_id\trepeat_haplotype\tfrom_entry_index\tto_entry_index"
             "\tspacer_from\tspacer_to\tdistinct_sequence_units\tdistinct_start_offsets\tstatus\n";
    for (std::size_t i = 0; i < result.groups.size(); ++i) {
        const auto& group = result.groups[i];
        if (group.entries.size() != group.positions.size())
            throw std::logic_error("Directed path position inventory differs");
        for (std::size_t position = 0; position < group.entries.size(); ++position)
            groups << "DirectedPath_" << i + 1 << "\tArray_" << group.parent + 1
                   << '\t' << group.haplotype << '\t' << group.kind << '\t'
                   << group.positions[position] << '\t' << group.entries[position] << '\n';
    }
    for (const auto& edge : result.edges) {
        const auto& [family, haplotype, from, to] = edge.key;
        edges << "Array_" << family + 1 << '\t' << haplotype << '\t'
              << from << '\t' << to << '\t' << edge.spacer_from << '\t'
              << edge.spacer_to << '\t' << edge.sequence_units << '\t'
              << edge.start_offsets << '\t' << edge.status << '\n';
    }
    result_bundle::close_checked(groups, groups_path);
    result_bundle::close_checked(edges, edges_path);
}

}} // namespace Zuck::directed_haplotype_paths

#endif
