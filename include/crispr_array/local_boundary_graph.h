#ifndef ZUCK_LOCAL_BOUNDARY_GRAPH_H
#define ZUCK_LOCAL_BOUNDARY_GRAPH_H

#include "crispr_array/read_recruitment.h"
#include <filesystem>
#include <iomanip>

namespace Zuck { namespace local_boundary_graph {

using Pair = std::pair<std::string, std::string>;
using Unit = read_recruitment::detail::Unit;

struct Config {
    std::size_t minimum_extension = 5;
    std::size_t maximum_extension = 16;
    std::size_t minimum_repeat = 20;
    std::size_t maximum_repeat = 55;
    std::size_t minimum_spacer = 26;
    std::size_t maximum_spacer = 50;
    std::size_t minimum_units = 2;
    std::size_t minimum_offsets = 2;
    std::size_t minimum_haplotype_spacers = 2;
    std::size_t minimum_within_haplotype_transitions = 2;
    std::size_t minimum_between_haplotype_edits = 3;
    std::size_t minimum_family_entries = 3;
    std::size_t minimum_family_spacer_after_recut = 20;
    std::size_t maximum_child_support_numerator = 1;
    std::size_t maximum_child_support_denominator = 2;
    std::size_t maximum_hits_per_read = 128;
    std::size_t maximum_pending_per_read = 512;
    std::size_t maximum_candidates = 65536;
};

struct Key {
    std::size_t owner = 0, entry = 0, left = 0, right = 0;
    Pair pair;
    bool operator<(const Key& other) const {
        return std::tie(owner, entry, left, right, pair) <
               std::tie(other.owner, other.entry, other.left, other.right, other.pair);
    }
};
struct TransitionKey {
    std::size_t owner = 0, left = 0, right = 0;
    std::string first, second;
    bool operator<(const TransitionKey& other) const {
        return std::tie(owner, left, right, first, second) <
               std::tie(other.owner, other.left, other.right, other.first, other.second);
    }
};
struct Evidence {
    std::set<Unit> units;
    std::set<std::size_t> offsets;
};
struct Decision {
    Key key;
    Pair original;
    std::size_t units = 0, offsets = 0, haplotype_spacers = 0;
    std::size_t within_transitions = 0, cross_transitions = 0;
    std::size_t geometry_entries = 0, child_entries = 0;
    bool owner_oriented = false;
    bool supported_haplotype = false, selected_geometry = false, replaced = false;
    std::string status;
};
struct Result {
    std::vector<Decision> decisions;
    std::map<std::string, std::size_t> counts;
};

inline std::size_t edit_distance(const std::string& a, const std::string& b) {
    return read_recruitment::detail::edit_distance(a, b);
}
inline std::string canonical(const std::string& sequence) {
    return read_recruitment::detail::canonical(sequence);
}
inline std::string core(const Pair& pair) {
    return read_recruitment::detail::core(pair);
}

template<class Limits>
inline bool family_candidate_allowed(const Pair& candidate, const Config& config,
                                     const Limits& limits) {
    return candidate.first.size() >= config.minimum_repeat &&
           candidate.first.size() <= config.maximum_repeat &&
           candidate.second.size() >= config.minimum_family_spacer_after_recut &&
           candidate.second.size() <= config.maximum_spacer &&
           limits.pair_allowed(candidate.first, candidate.second);
}

// A read-supported cut is generated only when the same added repeat bases are
// observed at two consecutive copies.  The old circular unit is bound before
// any replacement, so this module cannot invent a spacer or a primary row.
template<class Arrays, class Limits>
Result apply(Arrays& arrays, const std::vector<std::string>& files,
             const Limits& limits, const std::filesystem::path& directory,
             const Config& config = {}) {
    if (files.empty() || files.size() > 2)
        throw std::invalid_argument("Local boundary graph requires one or two read files");
    if (!config.minimum_extension || config.maximum_extension < config.minimum_extension ||
        config.maximum_repeat < config.minimum_repeat || config.maximum_spacer < config.minimum_spacer ||
        !config.minimum_units || !config.minimum_offsets || !config.minimum_haplotype_spacers ||
        !config.minimum_within_haplotype_transitions || !config.minimum_family_entries ||
        !config.minimum_family_spacer_after_recut ||
        config.maximum_spacer < config.minimum_family_spacer_after_recut ||
        !config.maximum_child_support_numerator ||
        !config.maximum_child_support_denominator ||
        config.maximum_child_support_numerator > config.maximum_child_support_denominator ||
        !config.maximum_hits_per_read ||
        !config.maximum_pending_per_read || !config.maximum_candidates)
        throw std::invalid_argument("Invalid local boundary graph configuration");

    Result result;
    result.counts = {{"local_boundary_graph_applied", 1}, {"local_boundary_graph_input_units", 0},
        {"local_boundary_graph_records", 0}, {"local_boundary_graph_records_without_quality", 0},
        {"local_boundary_graph_resource_limited_units", 0}, {"local_boundary_graph_candidates", 0},
        {"local_boundary_graph_supported_haplotypes", 0}, {"local_boundary_graph_selected_families", 0},
        {"local_boundary_graph_ambiguous_rows", 0}, {"local_boundary_graph_replaced_rows", 0},
        {"local_boundary_graph_cross_transition_geometries", 0}};

    struct Seed { std::string sequence; std::set<std::size_t> owners; };
    std::map<std::string, std::set<std::size_t>> seed_owners;
    std::vector<std::map<std::string, std::vector<std::size_t>>> original_cores(arrays.size());
    for (std::size_t owner = 0; owner < arrays.size(); ++owner) {
        for (std::size_t entry = 0; entry < arrays[owner].entries.size(); ++entry) {
            const auto& pair = arrays[owner].entries[entry];
            if (!read_recruitment::detail::dna(pair.first) || !read_recruitment::detail::dna(pair.second))
                throw std::invalid_argument("Local boundary graph rows must contain only ACGT");
            original_cores[owner][core(pair)].push_back(entry);
            if (pair.first.size() >= config.minimum_repeat && pair.first.size() <= config.maximum_repeat) {
                seed_owners[pair.first].insert(owner);
                seed_owners[read_recruitment::detail::rc(pair.first)].insert(owner);
            }
        }
    }
    std::vector<Seed> seeds;
    read_evidence::ExactMatcher matcher;
    for (const auto& item : seed_owners) {
        matcher.add(item.first, seeds.size());
        seeds.push_back({item.first, item.second});
    }
    matcher.build();

    std::map<Key, Evidence> evidence;
    std::map<TransitionKey, std::set<Unit>> transitions;
    bool global_overflow = false;
    auto scan = [&](const klibpp::KSeq& read, std::map<Key, std::set<std::size_t>>& pending,
                    std::map<TransitionKey, bool>& pending_transitions) {
        if (read.qual.empty()) return true;
        read_evidence::QualitySpans quality(read.qual, read_recruitment::MINIMUM_QUALITY);
        std::map<std::size_t, std::vector<std::size_t>> hits;
        std::size_t total_hits = 0;
        matcher.scan(read.seq, [&](std::size_t seed, std::size_t end) {
            if (++total_hits <= config.maximum_hits_per_read)
                hits[seed].push_back(end - seeds[seed].sequence.size());
        });
        if (total_hits > config.maximum_hits_per_read) return false;
        for (const auto& hit_list : hits) {
            const auto& seed = seeds[hit_list.first];
            const auto repeat_length = seed.sequence.size();
            for (std::size_t h = 1; h < hit_list.second.size(); ++h) {
                const auto first = hit_list.second[h - 1], second = hit_list.second[h];
                if (second < first + repeat_length + config.minimum_spacer ||
                    second > first + repeat_length + config.maximum_spacer) continue;
                const auto old_spacer_length = second - first - repeat_length;
                for (std::size_t extension = config.minimum_extension;
                     extension <= config.maximum_extension; ++extension) {
                    if (repeat_length + extension > config.maximum_repeat ||
                        old_spacer_length < config.minimum_spacer + extension) continue;
                    for (std::size_t left = 0; left <= extension; ++left) {
                        const auto right = extension - left;
                        if (left > first || second + repeat_length + right > read.seq.size()) continue;
                        const auto span_start = first - left;
                        const auto span_end = second + repeat_length + right;
                        if (!quality.passes(span_start, span_end,
                                read_recruitment::MAXIMUM_EXPECTED_ERRORS)) continue;
                        const auto first_signature = read.seq.substr(
                            first - left, repeat_length + extension);
                        const auto second_signature = read.seq.substr(
                            second - left, repeat_length + extension);
                        for (auto owner : seed.owners) {
                            pending_transitions[{owner, left, right,
                                first_signature, second_signature}] = true;
                            if (first_signature != second_signature) continue;
                            Pair candidate{first_signature,
                                read.seq.substr(first + repeat_length + right,
                                    old_spacer_length - extension)};
                            if (!limits.pair_allowed(candidate.first, candidate.second)) continue;
                            const auto found = original_cores[owner].find(core(candidate));
                            if (found == original_cores[owner].end()) continue;
                            std::set<Pair> old_pairs;
                            for (auto entry : found->second) old_pairs.insert(arrays[owner].entries[entry]);
                            if (old_pairs.size() != 1) continue;
                            for (auto entry : found->second)
                                pending[{owner, entry, left, right, candidate}].insert(first - left);
                            if (pending.size() > config.maximum_pending_per_read) return false;
                        }
                    }
                }
            }
        }
        return true;
    };

    read_linked::Reader first(files[0]);
    std::unique_ptr<read_linked::Reader> second;
    if (files.size() == 2) second = std::make_unique<read_linked::Reader>(files[1]);
    klibpp::KSeq a, b;
    while (first.next(a)) {
        if (second && (!second->next(b) || read_linked::read_name(a.name) != read_linked::read_name(b.name)))
            throw std::runtime_error("Local boundary graph read pairs have different names or lengths");
        ++result.counts["local_boundary_graph_input_units"];
        result.counts["local_boundary_graph_records"] += second ? 2 : 1;
        result.counts["local_boundary_graph_records_without_quality"] += a.qual.empty() + (second && b.qual.empty());
        if (global_overflow) continue;
        std::map<Key, std::set<std::size_t>> pending;
        std::map<TransitionKey, bool> pending_transitions;
        bool bounded = scan(a, pending, pending_transitions);
        if (second) bounded = scan(b, pending, pending_transitions) && bounded;
        if (!bounded) { ++result.counts["local_boundary_graph_resource_limited_units"]; continue; }
        const auto unit = read_recruitment::detail::unit_identity(a.seq, second ? b.seq : "", bool(second));
        for (const auto& item : pending) {
            if (evidence.size() >= config.maximum_candidates && !evidence.count(item.first)) {
                global_overflow = true; break;
            }
            auto& destination = evidence[item.first];
            destination.units.insert(unit);
            destination.offsets.insert(item.second.begin(), item.second.end());
        }
        if (global_overflow) { evidence.clear(); transitions.clear(); continue; }
        for (const auto& item : pending_transitions) transitions[item.first].insert(unit);
    }
    if (second && second->next(b))
        throw std::runtime_error("Local boundary graph read files contain different record counts");
    result.counts["local_boundary_graph_global_resource_limit"] = global_overflow;
    result.counts["local_boundary_graph_candidates"] = evidence.size();

    using Geometry = std::pair<std::size_t, std::size_t>;
    using OwnerGeometry = std::pair<std::size_t, Geometry>;
    std::map<OwnerGeometry, std::map<std::string, std::set<std::string>>> haplotype_spacers;
    for (const auto& item : evidence) {
        if (item.second.units.size() < config.minimum_units ||
            item.second.offsets.size() < config.minimum_offsets) continue;
        haplotype_spacers[{item.first.owner, {item.first.left, item.first.right}}]
            [canonical(item.first.pair.first)].insert(canonical(item.first.pair.second));
    }

    struct GeometryState {
        std::set<std::string> supported;
        std::size_t within = 0, cross = 0, total_spacers = 0;
        bool selected = false;
    };
    std::map<OwnerGeometry, GeometryState> states;
    for (const auto& item : haplotype_spacers) {
        auto& state = states[item.first];
        for (const auto& haplotype : item.second)
            if (haplotype.second.size() >= config.minimum_haplotype_spacers) {
                state.supported.insert(haplotype.first);
                state.total_spacers += haplotype.second.size();
            }
        bool separated = false;
        for (auto one = state.supported.begin(); one != state.supported.end(); ++one)
            for (auto two = std::next(one); two != state.supported.end(); ++two)
                separated |= edit_distance(*one, *two) >= config.minimum_between_haplotype_edits;
        if (!separated) continue;
        for (const auto& transition : transitions) {
            if (transition.first.owner != item.first.first ||
                transition.first.left != item.first.second.first ||
                transition.first.right != item.first.second.second) continue;
            const auto first_repeat = canonical(transition.first.first);
            const auto second_repeat = canonical(transition.first.second);
            if (!state.supported.count(first_repeat) || !state.supported.count(second_repeat)) continue;
            if (first_repeat == second_repeat) state.within += transition.second.size();
            else state.cross += transition.second.size();
        }
        if (state.cross) ++result.counts["local_boundary_graph_cross_transition_geometries"];
        state.selected = !state.cross && state.within >=
            config.minimum_within_haplotype_transitions * state.supported.size();
    }
    // Learn one owner-relative boundary from the whole family.  Seeds are
    // indexed on both strands, so only an extension that contains the row's
    // original repeat in the same orientation owns coordinates for that row.
    // The boundary is the last broadly supported geometry before both one-base
    // children lose at least half their distinct-row support.
    struct FamilyState {
        std::set<std::size_t> entries;
        std::map<std::string, std::set<std::string>> haplotype_spacers;
        std::map<std::string, std::size_t> haplotype_units;
        std::string dominant;
        std::size_t dominant_spacers = 0, child_entries = 0, cross = 0;
        bool eligible = false;
    };
    std::map<OwnerGeometry, FamilyState> family_states;
    for (const auto& item : evidence) {
        if (item.second.units.empty() || item.second.offsets.empty()) continue;
        const auto& original = arrays[item.first.owner].entries[item.first.entry];
        if (item.first.pair.first.find(original.first) == std::string::npos) continue;
        const OwnerGeometry geometry{item.first.owner, {item.first.left, item.first.right}};
        auto& state = family_states[geometry];
        state.entries.insert(item.first.entry);
        const auto haplotype = canonical(item.first.pair.first);
        state.haplotype_spacers[haplotype].insert(canonical(item.first.pair.second));
        state.haplotype_units[haplotype] += item.second.units.size();
    }
    for (auto& item : family_states) {
        auto& state = item.second;
        for (const auto& haplotype : state.haplotype_spacers) {
            const auto rank = std::make_tuple(haplotype.second.size(),
                state.haplotype_units[haplotype.first], haplotype.first);
            const auto old_rank = state.dominant.empty()
                ? std::make_tuple(std::size_t(0), std::size_t(0), std::string())
                : std::make_tuple(state.haplotype_spacers[state.dominant].size(),
                    state.haplotype_units[state.dominant], state.dominant);
            if (state.dominant.empty() || rank > old_rank) state.dominant = haplotype.first;
        }
        if (!state.dominant.empty())
            state.dominant_spacers = state.haplotype_spacers[state.dominant].size();
        const auto old = states.find(item.first);
        if (old != states.end()) state.cross = old->second.cross;
    }
    for (auto& item : family_states) {
        const auto owner = item.first.first;
        const auto left = item.first.second.first, right = item.first.second.second;
        auto child_entries = [&](std::size_t child_left, std::size_t child_right) {
            const auto child = family_states.find({owner, {child_left, child_right}});
            return child == family_states.end() ? std::size_t(0) : child->second.entries.size();
        };
        item.second.child_entries = std::max(
            child_entries(left + 1, right), child_entries(left, right + 1));
        item.second.eligible = left + right < config.maximum_extension &&
            item.second.entries.size() >= config.minimum_family_entries && !item.second.cross &&
            item.second.child_entries * config.maximum_child_support_denominator <=
                item.second.entries.size() * config.maximum_child_support_numerator;
    }

    std::map<std::string, std::set<std::size_t>> cluster_owners;
    std::map<std::string, std::size_t> cluster_support;
    for (const auto& item : family_states) if (item.second.eligible) {
        cluster_owners[item.second.dominant].insert(item.first.first);
        cluster_support[item.second.dominant] += item.second.dominant_spacers;
    }
    std::map<std::size_t, Geometry> selected_geometry;
    for (const auto& item : family_states) if (item.second.eligible) {
        const auto owner = item.first.first, geometry = item.first.second;
        const auto rank = std::make_tuple(item.second.entries.size(),
            item.second.dominant_spacers, cluster_owners[item.second.dominant].size(),
            cluster_support[item.second.dominant], geometry.first + geometry.second,
            std::numeric_limits<std::size_t>::max() - geometry.first);
        const auto old = selected_geometry.find(owner);
        if (old == selected_geometry.end()) selected_geometry[owner] = geometry;
        else {
            const auto& old_state = family_states.at({owner, old->second});
            const auto old_rank = std::make_tuple(old_state.entries.size(),
                old_state.dominant_spacers, cluster_owners[old_state.dominant].size(),
                cluster_support[old_state.dominant], old->second.first + old->second.second,
                std::numeric_limits<std::size_t>::max() - old->second.first);
            if (rank > old_rank) old->second = geometry;
        }
    }
    result.counts["local_boundary_graph_selected_families"] = selected_geometry.size();

    for (const auto& item : evidence) {
        Decision decision;
        decision.key = item.first;
        decision.original = arrays[item.first.owner].entries[item.first.entry];
        decision.units = item.second.units.size();
        decision.offsets = item.second.offsets.size();
        decision.owner_oriented = item.first.pair.first.find(decision.original.first) !=
            std::string::npos;
        const OwnerGeometry geometry{item.first.owner, {item.first.left, item.first.right}};
        const auto state = states.find(geometry);
        const auto haplotypes = haplotype_spacers.find(geometry);
        if (haplotypes != haplotype_spacers.end()) {
            const auto hap = haplotypes->second.find(canonical(item.first.pair.first));
            if (hap != haplotypes->second.end()) decision.haplotype_spacers = hap->second.size();
        }
        if (state != states.end()) {
            decision.within_transitions = state->second.within;
            decision.cross_transitions = state->second.cross;
            decision.supported_haplotype = state->second.supported.count(canonical(item.first.pair.first));
        }
        const auto family = family_states.find(geometry);
        if (family != family_states.end()) {
            decision.geometry_entries = family->second.entries.size();
            decision.child_entries = family->second.child_entries;
        }
        const auto selected = selected_geometry.find(item.first.owner);
        decision.selected_geometry = selected != selected_geometry.end() &&
            selected->second == geometry.second && decision.owner_oriented;
        decision.status = !decision.owner_oriented ? "opposite_strand_geometry" :
            decision.units < 1 ? "insufficient_sequence_units" :
            decision.offsets < 1 ? "insufficient_start_offsets" :
            decision.geometry_entries < config.minimum_family_entries ? "insufficient_family_entries" :
            decision.cross_transitions ? "cross_haplotype_transition" :
            family == family_states.end() || !family->second.eligible ? "no_boundary_cliff" :
            !decision.selected_geometry ? "nonselected_geometry" : "eligible";
        result.decisions.push_back(std::move(decision));
    }
    for (const auto& item : family_states)
        result.counts["local_boundary_graph_supported_haplotypes"] +=
            item.second.haplotype_spacers.size();

    struct RecutRecord { Pair original, candidate; Geometry geometry; };
    std::map<std::pair<std::size_t, std::size_t>, RecutRecord> family_recuts;
    std::set<std::pair<std::size_t, std::size_t>> replaced;
    for (const auto& selected : selected_geometry) {
        const auto owner = selected.first;
        const auto left = selected.second.first, right = selected.second.second;
        for (std::size_t entry = 0; entry < arrays[owner].entries.size(); ++entry) {
            const auto original = arrays[owner].entries[entry];
            if (original.second.size() < left + right + config.minimum_family_spacer_after_recut)
                continue;
            const auto end = original.second.size() - left;
            Pair candidate{
                (left ? original.second.substr(end) : std::string()) + original.first +
                    (right ? original.second.substr(0, right) : std::string()),
                original.second.substr(right, end - right)};
            if (!family_candidate_allowed(candidate, config, limits)) continue;
            family_recuts[{owner, entry}] = {original, candidate, selected.second};
            arrays[owner].entries[entry] = std::move(candidate);
            replaced.insert({owner, entry});
            ++result.counts["local_boundary_graph_replaced_rows"];
        }
    }
    std::set<std::pair<std::size_t, std::size_t>> ledger_replaced;
    for (auto& decision : result.decisions) {
        decision.replaced = replaced.count({decision.key.owner, decision.key.entry}) &&
            arrays[decision.key.owner].entries[decision.key.entry] == decision.key.pair;
        if (decision.replaced) {
            decision.status = "replaced";
            ledger_replaced.insert({decision.key.owner, decision.key.entry});
        }
        else if (decision.status == "eligible") decision.status = "family_row_not_recut";
    }
    for (const auto& item : family_recuts) if (!ledger_replaced.count(item.first)) {
        Decision decision;
        decision.key = {item.first.first, item.first.second,
            item.second.geometry.first, item.second.geometry.second, item.second.candidate};
        decision.original = item.second.original;
        const auto& family = family_states.at({item.first.first, item.second.geometry});
        decision.geometry_entries = family.entries.size();
        decision.child_entries = family.child_entries;
        decision.cross_transitions = family.cross;
        decision.owner_oriented = true;
        decision.selected_geometry = true;
        decision.replaced = true;
        decision.status = "family_recut";
        result.decisions.push_back(std::move(decision));
    }

    // Keep the family label aligned with the most common resolved row repeat.
    for (auto& array : arrays) {
        std::map<std::string, std::size_t> frequencies;
        for (const auto& pair : array.entries) ++frequencies[pair.first];
        if (!frequencies.empty())
            array.consensus = std::max_element(frequencies.begin(), frequencies.end(),
                [](const auto& a, const auto& b) {
                    return a.second != b.second ? a.second < b.second : a.first > b.first;
                })->first;
    }

    const auto path = directory / "local_boundary_graph.tsv";
    std::ofstream out(path);
    if (!out) throw std::runtime_error("Cannot write local boundary graph evidence");
    out << "owner_array_id\tentry_index\toriginal_repeat\toriginal_spacer\tcandidate_repeat\tcandidate_spacer"
           "\tleft_extension\tright_extension\tsequence_units\tstart_offsets\thaplotype_spacers"
           "\twithin_transitions\tcross_transitions\towner_oriented\tgeometry_entries\tchild_entries"
           "\tsupported_haplotype\tselected_geometry\treplaced\tstatus\n";
    for (const auto& d : result.decisions)
        out << "Array_" << d.key.owner + 1 << '\t' << d.key.entry + 1 << '\t'
            << d.original.first << '\t' << d.original.second << '\t' << d.key.pair.first << '\t'
            << d.key.pair.second << '\t' << d.key.left << '\t' << d.key.right << '\t'
            << d.units << '\t' << d.offsets << '\t' << d.haplotype_spacers << '\t'
            << d.within_transitions << '\t' << d.cross_transitions << '\t'
            << d.owner_oriented << '\t' << d.geometry_entries << '\t' << d.child_entries << '\t'
            << d.supported_haplotype << '\t' << d.selected_geometry << '\t'
            << d.replaced << '\t' << d.status << '\n';
    return result;
}

}} // namespace Zuck::local_boundary_graph
#endif
