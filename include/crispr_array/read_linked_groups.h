#ifndef ZUCK_READ_LINKED_GROUPS_H
#define ZUCK_READ_LINKED_GROUPS_H

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <fstream>
#include <map>
#include <memory>
#include <queue>
#include <set>
#include <stdexcept>
#include <string>
#include <tuple>
#include <unordered_map>
#include <utility>
#include <vector>
#include "kseq++/seqio.hpp"
#include "crispr_array/read_evidence.h"

namespace Zuck { namespace read_linked {

struct Config {
    std::size_t minimum_anchor = 16;
    int minimum_repeat = 20;
    int maximum_repeat = 55;
    int minimum_quality = 10;
    double maximum_expected_errors = 1.0;
    std::size_t minimum_sequence_units = 2;
    std::size_t minimum_start_offsets = 2;
    std::size_t minimum_component_spacers = 3;
    std::size_t scaffold_minimum_sequence_units = 2;
    std::size_t scaffold_minimum_start_offsets = 2;
};

struct Edge {
    std::string spacer_a, spacer_b;
    std::size_t sequence_units = 0, start_offsets = 0;
    bool selected = true;
    std::string selection_status = "selected_native_path";
};
struct Group {
    std::size_t parent = 0;
    std::vector<std::size_t> entry_indices; // One-based indexes in the original parent.
    std::vector<Edge> edges;
};
struct ComponentDecision {
    std::size_t parent = 0;
    std::vector<std::string> nodes;
    std::vector<Edge> edges;
    std::size_t entry_count = 0;
    std::size_t selected_group_count = 0;
    bool full_parent = false;
    std::string status;
};
struct Result {
    std::vector<Group> groups;
    // Complete evidence inventory for disjoint downstream grouping.  Unlike
    // the legacy overlapping diagnostic groups, this retains two-node and
    // full-parent nonbranched acyclic components.
    std::vector<Group> scaffold_groups;
    std::vector<ComponentDecision> scaffold_components;
    std::map<std::string, std::size_t> counts;
};

inline int base(char c) {
    return read_evidence::base(c);
}
inline std::string reverse_complement(const std::string& s) {
    std::string r;
    r.reserve(s.size());
    for (auto i = s.rbegin(); i != s.rend(); ++i) {
        const int b = base(*i);
        if (b < 0) return {};
        r.push_back("TGCA"[b]);
    }
    return r;
}
inline std::string canonical(const std::string& s) {
    const auto r = reverse_complement(s);
    return r.empty() ? std::string{} : std::min(s, r);
}
inline std::string read_name(std::string name) {
    const auto space = name.find_first_of(" \t");
    if (space != std::string::npos) name.resize(space);
    if (name.size() >= 2 && name[name.size() - 2] == '/' &&
        (name.back() == '1' || name.back() == '2')) name.resize(name.size() - 2);
    return name;
}

// The existing kseqpp parser supplies FASTA/FASTQ and gzip decoding. Inspect
// the first header separately so a FASTQ record lacking '+' cannot silently
// become a quality-free FASTA record. Parser and gzip errors remain fatal.
class Reader {
    std::string path_;
    bool fastq_ = false;
    using Stream = klibpp::KStreamIn<gzFile, int(*)(gzFile_s*, void*, unsigned int)>;
    std::unique_ptr<Stream> stream_;
    static int checked_read(gzFile file, void* buffer, unsigned int length) {
        const int size = gzread(file, buffer, length);
        int error = Z_OK;
        gzerror(file, &error);
        // zlib reports a missing gzip footer as EOF plus Z_BUF_ERROR. The
        // parser must receive an error, not a successful truncated stream.
        return error != Z_OK && error != Z_STREAM_END ? -1 : size;
    }
public:
    explicit Reader(const std::string& path) : path_(path) {
        gzFile file = gzopen(path.c_str(), "rb");
        if (!file) throw std::runtime_error("Cannot open grouping reads: " + path);
        char c = 0;
        int n = 0;
        do { n = gzread(file, &c, 1); }
        while (n == 1 && (c == '\n' || c == '\r' || c == ' ' || c == '\t'));
        const int close_status = gzclose(file);
        if (n < 0 || close_status != Z_OK)
            throw std::runtime_error("Cannot decode grouping reads: " + path);
        if (n == 1 && c != '@' && c != '>')
            throw std::runtime_error("Expected FASTA/FASTQ grouping reads: " + path);
        fastq_ = n == 1 && c == '@';
        gzFile input = gzopen(path.c_str(), "rb");
        if (!input) throw std::runtime_error("Cannot open grouping reads: " + path);
        stream_ = std::make_unique<Stream>(input, checked_read, gzclose);
    }
    bool next(klibpp::KSeq& record) {
        *stream_ >> record;
        if (stream_->err() || stream_->tqs())
            throw std::runtime_error("Malformed or unreadable grouping reads: " + path_);
        if (!*stream_) return false;
        if (record.name.empty() || record.seq.empty() ||
            (fastq_ && record.qual.size() != record.seq.size()) ||
            (!fastq_ && !record.qual.empty()))
            throw std::runtime_error("Malformed grouping read record: " + path_);
        for (char& c : record.seq)
            if (c >= 'a' && c <= 'z') c = static_cast<char>(c - 'a' + 'A');
        for (unsigned char c : record.qual)
            if (c < 33 || c > 126)
                throw std::runtime_error("Invalid grouping read quality: " + path_);
        return true;
    }
};

struct Endpoint { int spacer = 0; std::string sequence; };
using EdgeKey = std::tuple<std::size_t, int, int>;
struct Support {
    std::set<std::size_t> sequence_units;
    std::set<std::size_t> start_offsets;
};

// Infer supported content views only. The original arrays and every oriented
// row remain untouched; a view references original parent rows by index.
template <class Array>
Result infer(const std::vector<Array>& arrays, const std::vector<std::string>& files,
             const Config& config = {}) {
    if (files.empty() || files.size() > 2)
        throw std::runtime_error("Read-linked grouping requires one or two read files");
    if (config.minimum_anchor == 0 || config.minimum_component_spacers < 2 ||
        config.minimum_repeat <= 0 || config.maximum_repeat < config.minimum_repeat ||
        config.minimum_quality < 0 || config.minimum_quality > 93 ||
        !std::isfinite(config.maximum_expected_errors) || config.maximum_expected_errors < 0 ||
        config.minimum_sequence_units == 0 || config.minimum_start_offsets == 0 ||
        config.scaffold_minimum_sequence_units == 0 ||
        config.scaffold_minimum_start_offsets == 0)
        throw std::runtime_error("Invalid read-linked grouping configuration");
    Result result;
    result.counts = {{"read_linked_input_files", files.size()},
        {"read_linked_records", 0}, {"read_linked_records_without_quality", 0},
        {"read_linked_input_units", 0}, {"read_linked_distinct_support_sequence_units", 0},
        {"read_linked_duplicate_support_sequence_units", 0},
        {"read_linked_supported_edges", 0}, {"read_linked_branched_components", 0},
        {"read_linked_groups", 0}, {"read_linked_memberships", 0},
        {"read_linked_scaffold_groups", 0},
        {"read_linked_scaffold_memberships", 0},
        {"read_linked_scaffold_supported_edges", 0},
        {"read_linked_scaffold_components", 0},
        {"read_linked_scaffold_branched_components", 0},
        {"read_linked_scaffold_cycle_components", 0},
        {"read_linked_scaffold_forest_components", 0},
        {"read_linked_scaffold_forest_groups", 0},
        {"read_linked_scaffold_forest_selected_edges", 0},
        {"read_linked_scaffold_forest_rejected_degree_edges", 0},
        {"read_linked_scaffold_forest_rejected_cycle_edges", 0}};

    std::set<std::string> unique_spacers;
    for (const auto& array : arrays) for (const auto& row : array.entries) {
        const auto s = canonical(row.second);
        if (!s.empty()) unique_spacers.insert(s);
    }
    std::vector<std::string> spacers(unique_spacers.begin(), unique_spacers.end());
    std::map<std::string, int> spacer_ids;
    for (std::size_t i = 0; i < spacers.size(); ++i) spacer_ids[spacers[i]] = static_cast<int>(i);
    std::vector<std::set<int>> family_spacers(arrays.size());
    std::map<std::string, std::set<std::size_t>> repeat_families;
    using EndpointIndex = std::unordered_map<std::string, std::vector<Endpoint>>;
    std::vector<EndpointIndex> left(arrays.size()), right(arrays.size());
    for (std::size_t family = 0; family < arrays.size(); ++family) {
        std::set<std::string> repeats{canonical(arrays[family].consensus)};
        for (const auto& row : arrays[family].entries) {
            const auto s = canonical(row.second);
            if (!s.empty()) family_spacers[family].insert(spacer_ids.at(s));
            repeats.insert(canonical(row.first));
        }
        for (const auto& r : repeats)
            if (r.size() >= static_cast<std::size_t>(config.minimum_repeat) &&
                r.size() <= static_cast<std::size_t>(config.maximum_repeat))
                repeat_families[r].insert(family);
        for (int sid : family_spacers[family]) {
            const std::set<std::string> orientations{spacers[sid], reverse_complement(spacers[sid])};
            for (const auto& s : orientations) if (s.size() >= config.minimum_anchor) {
                left[family][s.substr(s.size() - config.minimum_anchor)].push_back({sid, s});
                right[family][s.substr(0, config.minimum_anchor)].push_back({sid, s});
            }
        }
    }

    read_evidence::ExactMatcher matcher;
    std::vector<std::string> repeats;
    std::vector<std::set<std::size_t>> owners;
    for (const auto& item : repeat_families) {
        const auto id = repeats.size();
        repeats.push_back(item.first); owners.push_back(item.second);
        matcher.add(item.first, id);
        const auto reverse = reverse_complement(item.first);
        if (reverse != item.first) matcher.add(reverse, id);
    }
    matcher.build();
    std::map<EdgeKey, Support> support;
    std::map<std::pair<std::string, std::string>, std::size_t> sequence_units;
    auto scan = [&](const klibpp::KSeq& read, std::map<EdgeKey, std::set<std::size_t>>& links) {
        ++result.counts["read_linked_records"];
        if (read.qual.empty()) { ++result.counts["read_linked_records_without_quality"]; return; }
        read_evidence::QualitySpans quality(read.qual, config.minimum_quality);
        matcher.scan(read.seq, [&](std::size_t repeat_id, std::size_t end) {
            const std::size_t start = end - repeats[repeat_id].size();
            if (start < config.minimum_anchor || config.minimum_anchor > read.seq.size() - end)
                return;
            const auto left_anchor = read.seq.substr(start - config.minimum_anchor, config.minimum_anchor);
            const auto right_anchor = read.seq.substr(end, config.minimum_anchor);
            for (std::size_t family : owners[repeat_id]) {
                const auto left_candidates = left[family].find(left_anchor);
                const auto right_candidates = right[family].find(right_anchor);
                if (left_candidates == left[family].end() || right_candidates == right[family].end())
                    continue;

                // A shared anchor cannot choose between different spacer types.
                // Multiple orientations of the same spacer remain eligible.
                std::set<int> left_ids, right_ids;
                for (const auto& endpoint : left_candidates->second) left_ids.insert(endpoint.spacer);
                for (const auto& endpoint : right_candidates->second) right_ids.insert(endpoint.spacer);
                if (left_ids.size() != 1 || right_ids.size() != 1 || *left_ids.begin() == *right_ids.begin())
                    continue;

                for (const auto& left_endpoint : left_candidates->second) {
                    for (const auto& right_endpoint : right_candidates->second) {
                        const auto left_length = std::min(start, left_endpoint.sequence.size());
                        const auto right_length = std::min(read.seq.size() - end, right_endpoint.sequence.size());
                        // Check every available flank base, not just the anchor.
                        const bool left_matches = read.seq.compare(start - left_length, left_length,
                            left_endpoint.sequence, left_endpoint.sequence.size() - left_length, left_length) == 0;
                        const bool right_matches = read.seq.compare(end, right_length,
                            right_endpoint.sequence, 0, right_length) == 0;
                        if (!left_matches || !right_matches ||
                            !quality.passes(start - left_length, end + right_length, config.maximum_expected_errors))
                            continue;
                        const auto first_spacer = std::min(left_endpoint.spacer, right_endpoint.spacer);
                        const auto second_spacer = std::max(left_endpoint.spacer, right_endpoint.spacer);
                        links[{family, first_spacer, second_spacer}].insert(start);
                    }
                }
            }
        });
    };

    Reader first(files[0]);
    std::unique_ptr<Reader> second;
    if (files.size() == 2) second = std::make_unique<Reader>(files[1]);
    klibpp::KSeq a, b;
    while (first.next(a)) {
        if (second && (!second->next(b) || read_name(a.name) != read_name(b.name)))
            throw std::runtime_error("Grouping read pairs have different names or lengths");
        std::map<EdgeKey, std::set<std::size_t>> links;
        scan(a, links); if (second) scan(b, links);
        ++result.counts["read_linked_input_units"];
        if (links.empty()) continue;
        const auto unit = sequence_units.emplace(std::make_pair(a.seq, second ? b.seq : ""), sequence_units.size());
        if (!unit.second) ++result.counts["read_linked_duplicate_support_sequence_units"];
        for (const auto& link : links) {
            auto& evidence = support[link.first];
            evidence.sequence_units.insert(unit.first->second);
            evidence.start_offsets.insert(link.second.begin(), link.second.end());
        }
    }
    if (second && second->next(b)) throw std::runtime_error("Grouping read files contain different record counts");
    result.counts["read_linked_distinct_support_sequence_units"] = sequence_units.size();

    std::vector<std::map<int, std::set<int>>> graphs(arrays.size());
    std::vector<std::map<int, std::set<int>>> scaffold_graphs(arrays.size());
    for (const auto& item : support) {
        const auto [family, a_id, b_id] = item.first;
        if (item.second.sequence_units.size() >= config.minimum_sequence_units &&
            item.second.start_offsets.size() >= config.minimum_start_offsets) {
            graphs[family][a_id].insert(b_id); graphs[family][b_id].insert(a_id);
            ++result.counts["read_linked_supported_edges"];
        }
        if (item.second.sequence_units.size() >= config.scaffold_minimum_sequence_units &&
            item.second.start_offsets.size() >= config.scaffold_minimum_start_offsets) {
            scaffold_graphs[family][a_id].insert(b_id);
            scaffold_graphs[family][b_id].insert(a_id);
            ++result.counts["read_linked_scaffold_supported_edges"];
        }
    }

    auto components = [](const std::map<int, std::set<int>>& graph) {
        std::vector<std::set<int>> answer;
        std::set<int> remaining;
        for (const auto& item : graph) remaining.insert(item.first);
        while (!remaining.empty()) {
            std::set<int> component;
            std::vector<int> pending{*remaining.begin()};
            while (!pending.empty()) {
                const int node = pending.back(); pending.pop_back();
                if (!component.insert(node).second) continue;
                const auto& neighbors = graph.at(node);
                pending.insert(pending.end(), neighbors.begin(), neighbors.end());
            }
            for (int node : component) remaining.erase(node);
            answer.push_back(std::move(component));
        }
        return answer;
    };
    auto make_group = [&](std::size_t family, const std::set<int>& component,
                          const std::map<int, std::set<int>>& graph) {
        Group group; group.parent = family;
        for (std::size_t i = 0; i < arrays[family].entries.size(); ++i) {
            const auto s = canonical(arrays[family].entries[i].second);
            const auto id = spacer_ids.find(s);
            if (id != spacer_ids.end() && component.count(id->second))
                group.entry_indices.push_back(i + 1);
        }
        for (int node : component) for (int neighbor : graph.at(node)) if (node < neighbor) {
            const auto& evidence = support.at({family, node, neighbor});
            group.edges.push_back({spacers[node], spacers[neighbor],
                evidence.sequence_units.size(), evidence.start_offsets.size()});
        }
        return group;
    };

    // Preserve the legacy overlapping diagnostic view byte-for-byte in
    // meaning: partial components with at least three spacer nodes at 2/2.
    for (std::size_t family = 0; family < arrays.size(); ++family) {
        const auto& graph = graphs[family];
        for (const auto& component : components(graph)) {
            if (component.size() < config.minimum_component_spacers ||
                component == family_spacers[family]) continue;
            std::size_t degree_sum = 0; bool branched = false;
            for (int node : component) { degree_sum += graph.at(node).size(); branched |= graph.at(node).size() > 2; }
            if (branched) { ++result.counts["read_linked_branched_components"]; continue; }
            if (degree_sum / 2 != component.size() - 1) continue; // No cyclic or arbitrarily cut views.
            auto group = make_group(family, component, graph);
            result.counts["read_linked_memberships"] += group.entry_indices.size();
            result.groups.push_back(std::move(group));
        }
    }

    // The primary disjoint scaffolder receives a separate 2/2 graph and an
    // explicit decision for every connected component, including full-parent
    // paths and two-node bridges that the diagnostic overlay intentionally omits.
    // Branched or cyclic components are resolved deterministically into a
    // maximum-evidence degree-at-most-two acyclic forest.  Every discarded edge
    // remains in the component ledger with its rejection reason.
    for (std::size_t family = 0; family < arrays.size(); ++family) {
        const auto& graph = scaffold_graphs[family];
        for (const auto& component : components(graph)) {
            auto group = make_group(family, component, graph);
            std::size_t degree_sum = 0; bool branched = false;
            for (int node : component) {
                degree_sum += graph.at(node).size();
                branched |= graph.at(node).size() > 2;
            }
            ComponentDecision decision;
            decision.parent = family;
            for (int node : component) decision.nodes.push_back(spacers[node]);
            decision.edges = group.edges;
            decision.entry_count = group.entry_indices.size();
            decision.full_parent = component == family_spacers[family];
            const bool cyclic = degree_sum / 2 != component.size() - 1;
            if (!branched && !cyclic) {
                decision.status = component.size() == 2 ? "accepted_two_node" :
                    decision.full_parent ? "accepted_full_parent_path" :
                    "accepted_partial_path";
                decision.selected_group_count = 1;
                result.counts["read_linked_scaffold_memberships"] +=
                    group.entry_indices.size();
                result.scaffold_groups.push_back(std::move(group));
            } else {
                if (branched) {
                    decision.status = "resolved_branched_forest";
                    ++result.counts["read_linked_scaffold_branched_components"];
                } else {
                    decision.status = "resolved_cycle_forest";
                    ++result.counts["read_linked_scaffold_cycle_components"];
                }
                ++result.counts["read_linked_scaffold_forest_components"];

                struct Candidate {
                    int left = 0, right = 0;
                    std::size_t units = 0, offsets = 0;
                };
                std::vector<Candidate> candidates;
                for (int node : component) {
                    for (int neighbor : graph.at(node)) if (node < neighbor) {
                        const auto& evidence = support.at({family, node, neighbor});
                        candidates.push_back({node, neighbor, evidence.sequence_units.size(),
                                              evidence.start_offsets.size()});
                    }
                }
                std::sort(candidates.begin(), candidates.end(), [](const Candidate& left,
                                                                    const Candidate& right) {
                    const auto left_rank = std::make_tuple(
                        std::min(left.units, left.offsets), left.units + left.offsets,
                        std::max(left.units, left.offsets), left.units, left.offsets);
                    const auto right_rank = std::make_tuple(
                        std::min(right.units, right.offsets), right.units + right.offsets,
                        std::max(right.units, right.offsets), right.units, right.offsets);
                    if (left_rank != right_rank) return left_rank > right_rank;
                    return std::make_pair(left.left, left.right) <
                           std::make_pair(right.left, right.right);
                });

                std::map<int, int> roots;
                std::map<int, std::size_t> degrees;
                for (int node : component) roots[node] = node;
                auto find_root = [&](int value) {
                    int root = value;
                    while (roots.at(root) != root) root = roots.at(root);
                    while (roots.at(value) != value) {
                        const int following = roots.at(value);
                        roots[value] = root;
                        value = following;
                    }
                    return root;
                };
                std::map<int, std::set<int>> selected_graph;
                std::map<std::pair<std::string, std::string>, std::string> edge_status;
                for (const auto& candidate : candidates) {
                    auto key = std::make_pair(spacers[candidate.left],
                                              spacers[candidate.right]);
                    if (key.second < key.first) std::swap(key.first, key.second);
                    if (degrees[candidate.left] >= 2 || degrees[candidate.right] >= 2) {
                        edge_status[key] = "rejected_degree_limit";
                        ++result.counts[
                            "read_linked_scaffold_forest_rejected_degree_edges"];
                        continue;
                    }
                    int left_root = find_root(candidate.left);
                    int right_root = find_root(candidate.right);
                    if (left_root == right_root) {
                        edge_status[key] = "rejected_cycle";
                        ++result.counts[
                            "read_linked_scaffold_forest_rejected_cycle_edges"];
                        continue;
                    }
                    if (right_root < left_root) std::swap(left_root, right_root);
                    roots[right_root] = left_root;
                    ++degrees[candidate.left]; ++degrees[candidate.right];
                    selected_graph[candidate.left].insert(candidate.right);
                    selected_graph[candidate.right].insert(candidate.left);
                    edge_status[key] = "selected_max_evidence_forest";
                    ++result.counts["read_linked_scaffold_forest_selected_edges"];
                }
                for (auto& edge : decision.edges) {
                    auto key = std::make_pair(edge.spacer_a, edge.spacer_b);
                    if (key.second < key.first) std::swap(key.first, key.second);
                    edge.selection_status = edge_status.at(key);
                    edge.selected = edge.selection_status ==
                        "selected_max_evidence_forest";
                }
                for (const auto& forest_component : components(selected_graph)) {
                    auto forest_group = make_group(
                        family, forest_component, selected_graph);
                    for (auto& edge : forest_group.edges)
                        edge.selection_status = "selected_max_evidence_forest";
                    result.counts["read_linked_scaffold_memberships"] +=
                        forest_group.entry_indices.size();
                    result.scaffold_groups.push_back(std::move(forest_group));
                    ++decision.selected_group_count;
                    ++result.counts["read_linked_scaffold_forest_groups"];
                }
            }
            result.scaffold_components.push_back(std::move(decision));
        }
    }
    result.counts["read_linked_groups"] = result.groups.size();
    result.counts["read_linked_scaffold_groups"] = result.scaffold_groups.size();
    result.counts["read_linked_scaffold_components"] = result.scaffold_components.size();
    return result;
}

inline void write_tables(const std::string& directory, const Result& result) {
    const std::string members_path = directory + "/read_linked_groups.tsv";
    const std::string edges_path = directory + "/read_linked_edges.tsv";
    std::ofstream members(members_path), edges(edges_path);
    if (!members || !edges) throw std::runtime_error("Cannot stage read-linked grouping tables");
    members << "group_id\tparent_array_id\tentry_index\n";
    edges << "group_id\tparent_array_id\tspacer_a\tspacer_b\tdistinct_sequence_units\tdistinct_start_offsets\tquality_mode\n";
    for (std::size_t i = 0; i < result.groups.size(); ++i) {
        const auto& group = result.groups[i];
        const auto id = "ReadLinked_" + std::to_string(i + 1);
        const auto parent = "Array_" + std::to_string(group.parent + 1);
        for (std::size_t index : group.entry_indices) members << id << '\t' << parent << '\t' << index << '\n';
        for (const auto& edge : group.edges)
            edges << id << '\t' << parent << '\t' << edge.spacer_a << '\t' << edge.spacer_b << '\t'
                  << edge.sequence_units << '\t' << edge.start_offsets << "\tee1_q10_exact_flanks\n";
    }
    members.close(); edges.close();
    if (!members || !edges) throw std::runtime_error("Cannot finish read-linked grouping tables");

    const std::string scaffold_members_path = directory + "/read_linked_scaffolds.tsv";
    const std::string scaffold_edges_path = directory + "/read_linked_scaffold_edges.tsv";
    std::ofstream scaffold_members(scaffold_members_path), scaffold_edges(scaffold_edges_path);
    if (!scaffold_members || !scaffold_edges)
        throw std::runtime_error("Cannot stage read-linked scaffold tables");
    scaffold_members << "group_id\tparent_array_id\tentry_index\n";
    scaffold_edges << "group_id\tparent_array_id\tspacer_a\tspacer_b"
                      "\tdistinct_sequence_units\tdistinct_start_offsets\tquality_mode\n";
    for (std::size_t i = 0; i < result.scaffold_groups.size(); ++i) {
        const auto& group = result.scaffold_groups[i];
        const auto id = "ReadScaffold_" + std::to_string(i + 1);
        const auto parent = "Array_" + std::to_string(group.parent + 1);
        for (std::size_t index : group.entry_indices)
            scaffold_members << id << '\t' << parent << '\t' << index << '\n';
        for (const auto& edge : group.edges)
            scaffold_edges << id << '\t' << parent << '\t' << edge.spacer_a << '\t'
                           << edge.spacer_b << '\t' << edge.sequence_units << '\t'
                           << edge.start_offsets << "\tee1_q10_exact_flanks\n";
    }
    scaffold_members.close(); scaffold_edges.close();
    if (!scaffold_members || !scaffold_edges)
        throw std::runtime_error("Cannot finish read-linked scaffold tables");

    const std::string components_path = directory + "/read_linked_scaffold_components.tsv";
    const std::string component_edges_path =
        directory + "/read_linked_scaffold_component_edges.tsv";
    std::ofstream component_rows(components_path), component_edges(component_edges_path);
    if (!component_rows || !component_edges)
        throw std::runtime_error("Cannot stage read-linked scaffold decision tables");
    component_rows << "component_id\tparent_array_id\tstatus\tfull_parent"
                      "\tnode_count\tentry_count\tedge_count\tselected_edge_count"
                      "\tselected_group_count"
                      "\tminimum_sequence_units\tminimum_start_offsets\n";
    component_edges << "component_id\tparent_array_id\tstatus\tspacer_a\tspacer_b"
                       "\tdistinct_sequence_units\tdistinct_start_offsets"
                       "\tselected\tselection_status\n";
    for (std::size_t i = 0; i < result.scaffold_components.size(); ++i) {
        const auto& component = result.scaffold_components[i];
        const auto id = "ReadComponent_" + std::to_string(i + 1);
        const auto parent = "Array_" + std::to_string(component.parent + 1);
        std::size_t minimum_units = 0, minimum_offsets = 0, selected_edges = 0;
        for (const auto& edge : component.edges) {
            minimum_units = !minimum_units ? edge.sequence_units :
                std::min(minimum_units, edge.sequence_units);
            minimum_offsets = !minimum_offsets ? edge.start_offsets :
                std::min(minimum_offsets, edge.start_offsets);
            component_edges << id << '\t' << parent << '\t' << component.status
                            << '\t' << edge.spacer_a << '\t' << edge.spacer_b << '\t'
                            << edge.sequence_units << '\t' << edge.start_offsets << '\t'
                            << edge.selected << '\t' << edge.selection_status << '\n';
            selected_edges += edge.selected;
        }
        component_rows << id << '\t' << parent << '\t' << component.status << '\t'
                       << component.full_parent << '\t' << component.nodes.size() << '\t'
                       << component.entry_count << '\t' << component.edges.size() << '\t'
                       << selected_edges << '\t' << component.selected_group_count << '\t'
                       << minimum_units << '\t' << minimum_offsets << '\n';
    }
    component_rows.close(); component_edges.close();
    if (!component_rows || !component_edges)
        throw std::runtime_error("Cannot finish read-linked scaffold decision tables");
}

} } // namespace Zuck::read_linked
#endif
