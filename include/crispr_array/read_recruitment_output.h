#ifndef ZUCK_READ_RECRUITMENT_OUTPUT_H
#define ZUCK_READ_RECRUITMENT_OUTPUT_H

#include "crispr_array/read_recruitment.h"
#include "crispr_array/supported_calls.h"
#include <iomanip>
#include <tuple>

namespace Zuck { namespace read_recruitment_output {

inline int kind_rank(const std::string& kind) {
    if (kind == "full_pair") return 3;
    if (kind == "full_partial") return 2;
    if (kind == "graph_junction") return 1;
    throw std::invalid_argument("Unknown read-recruitment evidence kind");
}

inline std::string selection_status(const read_recruitment::Candidate& candidate) {
    if (!candidate.accepted) return "absent";
    if (candidate.selection_eligible) return "supported";
    return candidate.tandem_ambiguous ? "tandem_ambiguous" : "error_ambiguous";
}

// Schema four constrains automatic promotion, not the frozen schema-three
// discovery/error/tandem evidence. Context is fixed before applying this policy.
struct ContextWitness {
    std::size_t candidate_id = 0, entry_index = 0;
    int repeat_distance = -1, spacer_distance = -1;
    std::size_t spacer_denominator = 0;
    bool repeat_reverse_complement = false, spacer_reverse_complement = false;
};
using OwnerContext = std::map<std::size_t, ContextWitness>;

inline int hamming_distance(const std::string& a, const std::string& b) {
    if (a.size() != b.size()) return -1;
    int distance = 0;
    for (std::size_t i = 0; i < a.size(); ++i) distance += a[i] != b[i];
    return distance;
}

inline ContextWitness context_witness(const read_recruitment::Candidate& focal,
        const read_recruitment::Candidate& other, std::size_t candidate_id,
        std::size_t entry_index) {
    ContextWitness result;
    if (focal.left_repeat.size() != other.left_repeat.size()) return result;
    const auto repeat_forward = hamming_distance(focal.left_repeat, other.left_repeat);
    const auto repeat_reverse = hamming_distance(focal.left_repeat,
        read_recruitment::detail::rc(other.left_repeat));
    const auto repeat_distance = std::min(repeat_forward, repeat_reverse);
    if (repeat_distance > 1) return result;
    const auto spacer_forward = read_recruitment::detail::edit_distance(focal.spacer, other.spacer);
    const auto spacer_reverse = read_recruitment::detail::edit_distance(focal.spacer,
        read_recruitment::detail::rc(other.spacer));
    const auto spacer_distance = std::min(spacer_forward, spacer_reverse);
    const auto denominator = std::max(focal.spacer.size(), other.spacer.size());
    // Inclusive .18 dissimilarity, using strand-invariant edit distance.
    if (50 * spacer_distance < 9 * denominator) return result;
    result.candidate_id = candidate_id;
    result.entry_index = entry_index;
    result.repeat_distance = repeat_distance;
    result.repeat_reverse_complement = repeat_reverse < repeat_forward;
    result.spacer_distance = static_cast<int>(spacer_distance);
    result.spacer_denominator = denominator;
    result.spacer_reverse_complement = spacer_reverse < spacer_forward;
    return result;
}

inline std::size_t supported_owners(const OwnerContext& context) {
    return std::count_if(context.begin(), context.end(),
        [](const auto& item) { return item.second.candidate_id != 0; });
}

inline std::string row_selection_status(const read_recruitment::Candidate& candidate,
        bool added, const OwnerContext& context) {
    if (added) return "novel_boundary_unresolved";
    if (!candidate.selection_eligible) return selection_status(candidate);
    if (context.empty() || supported_owners(context) != context.size())
        return "insufficient_spacer_diversity";
    return "supported";
}

// Preserve every original row and local row index. New exact observed pairs
// append in lexical order within one unambiguous graph family. This does not
// claim that a repeat family is a physical locus.
template<class Arrays, class Limits>
supported_calls::RecruitmentEvidence apply(Arrays& arrays,
        const read_recruitment::Result& result, const Limits& limits,
        const std::filesystem::path& directory,
        std::map<std::string, std::size_t>& counts) {
    using Pair = read_recruitment::Pair;
    using Key = std::pair<std::size_t, Pair>;
    std::vector<std::map<Pair, std::size_t>> original(arrays.size());
    std::map<std::string, std::set<Pair>> original_cores;
    std::size_t original_rows = 0;
    for (std::size_t a = 0; a < arrays.size(); ++a) {
        for (std::size_t e = 0; e < arrays[a].entries.size(); ++e) {
            const auto& pair = arrays[a].entries[e];
            original[a].emplace(pair, e + 1);
            original_cores[element_confidence::circular_identity(pair.first + pair.second)]
                .insert(element_confidence::pair_identity(pair));
            ++original_rows;
        }
    }
    counts["read_recruitment_seed_families"] = arrays.size();
    counts["read_recruitment_original_rows"] = original_rows;
    counts["read_recruitment_added_rows"] = 0;
    counts["read_recruitment_bound_rows"] = 0;
    counts["read_recruitment_ambiguous_new_pairs"] = 0;
    counts["read_recruitment_protected_core_pairs"] = 0;
    counts["read_recruitment_length_excluded_pairs"] = 0;
    counts["read_recruitment_original_snapshot_rows"] = original_rows;
    counts["read_recruitment_context_rows"] = 0;
    counts["read_recruitment_context_supported_owners"] = 0;
    counts["read_recruitment_policy_eligible_rows"] = 0;
    counts["read_recruitment_novel_boundary_unresolved_rows"] = 0;
    counts["read_recruitment_insufficient_spacer_diversity_rows"] = 0;
    std::map<Key, std::size_t> proofs;
    auto better = [&](std::size_t left, std::size_t right) {
        const auto& a = result.candidates[left];
        const auto& b = result.candidates[right];
        auto rank = [](const auto& c) {
            return std::make_tuple(kind_rank(c.evidence_kind), c.support.sequence_units,
                                   c.support.start_offsets);
        };
        return rank(a) > rank(b) || (rank(a) == rank(b) && left < right);
    };
    for (std::size_t id = 0; id < result.candidates.size(); ++id) {
        const auto& c = result.candidates[id];
        if (!c.accepted) continue;
        if (!c.pair_only || !c.support.accepted() || c.distinct_supported_family_spacers < 2 ||
            c.boundary_conflict || c.resource_limited || c.family_owners.empty())
            throw std::invalid_argument("Inconsistent accepted recruitment candidate");
        kind_rank(c.evidence_kind);
        for (const auto owner : c.family_owners) {
            if (owner >= arrays.size()) throw std::invalid_argument("Invalid recruitment family owner");
            for (const auto& pair : read_recruitment::oriented_pairs(c)) {
                if (!limits.pair_allowed(pair.first, pair.second)) {
                    ++counts["read_recruitment_length_excluded_pairs"]; continue;
                }
                const bool exists = original[owner].count(pair);
                if (!exists && c.family_owners.size() != 1) {
                    ++counts["read_recruitment_ambiguous_new_pairs"]; continue;
                }
                if (!exists && c.evidence_kind == "graph_junction") continue;
                const auto core = original_cores.find(element_confidence::circular_identity(pair.first + pair.second));
                if (!exists && core != original_cores.end() &&
                    !core->second.count(element_confidence::pair_identity(pair))) {
                    ++counts["read_recruitment_protected_core_pairs"]; continue;
                }
                const Key key{owner, pair};
                const auto old = proofs.find(key);
                if (old == proofs.end() || better(id, old->second)) proofs[key] = id;
            }
        }
    }
    // Context sources must have an actual proof binding to an exact ORIGINAL
    // row in this owner. Candidate.graph_proposed also recognizes RC aliases
    // and other owners and therefore cannot establish this membership.
    std::map<std::size_t, std::map<std::size_t, std::size_t>> context_sources;
    std::set<std::size_t> bound_candidates;
    for (const auto& [key, id] : proofs) {
        bound_candidates.insert(id);
        const auto found = original[key.first].find(key.second);
        if (found != original[key.first].end() && result.candidates[id].selection_eligible)
            context_sources[key.first].emplace(id, found->second);
    }
    std::map<std::size_t, OwnerContext> contexts;
    for (auto id : bound_candidates) {
        const auto& focal = result.candidates[id];
        auto& context = contexts[id];
        for (auto owner : focal.family_owners) {
            auto& witness = context[owner];
            const auto sources = context_sources.find(owner);
            if (sources == context_sources.end()) continue;
            // IDs ascend independently of row lexical order. Do not prune this
            // source set by the final policy: that would create circular gates.
            for (const auto& [other, entry_index] : sources->second) {
                witness = context_witness(focal, result.candidates[other], other + 1, entry_index);
                if (witness.candidate_id) break;
            }
        }
    }
    supported_calls::RecruitmentEvidence evidence;
    const auto candidate_path = directory / "read_recruitment_candidates.tsv";
    const auto proof_path = directory / "read_recruitment.tsv";
    const auto original_path = directory / "read_recruitment_original_rows.tsv";
    const auto context_path = directory / "read_recruitment_context.tsv";
    std::ofstream candidates(candidate_path), out(proof_path), snapshot(original_path), context_out(context_path);
    if (!candidates || !out || !snapshot || !context_out)
        throw std::runtime_error("Cannot write recruitment evidence");
    snapshot << "array_id\tentry_index\trepeat_variant\tspacer\n";
    for (std::size_t owner = 0; owner < arrays.size(); ++owner)
        for (std::size_t entry = 0; entry < arrays[owner].entries.size(); ++entry) {
            const auto& pair = arrays[owner].entries[entry];
            snapshot << "Array_" << owner + 1 << '\t' << entry + 1 << '\t'
                << pair.first << '\t' << pair.second << '\n';
        }
    context_out << "candidate_id\towner_array_id\tcontext_candidate_id\tcontext_entry_index"
        "\trepeat_distance\trepeat_reverse_complement\tspacer_distance\tspacer_denominator"
        "\tspacer_reverse_complement\tstatus\n";
    for (const auto& [id, context] : contexts) for (const auto& [owner, witness] : context) {
        context_out << id + 1 << "\tArray_" << owner + 1 << '\t' << witness.candidate_id
            << '\t' << witness.entry_index << '\t' << witness.repeat_distance
            << '\t' << witness.repeat_reverse_complement << '\t' << witness.spacer_distance
            << '\t' << witness.spacer_denominator << '\t' << witness.spacer_reverse_complement
            << '\t' << (witness.candidate_id ? "supported" : "missing") << '\n';
        ++counts["read_recruitment_context_rows"];
        counts["read_recruitment_context_supported_owners"] += witness.candidate_id != 0;
    }
    candidates << std::setprecision(17);
    candidates << "candidate_id\tfamily_owners\tleft_repeat\tspacer\tright_repeat\tevidence_kind"
        "\tgraph_proposed\tleft_complete\tright_complete\tleft_observed_bases\tright_observed_bases"
        "\tsequence_units\tstart_offsets\tunits_saturated\toffsets_saturated"
        "\tdistinct_supported_family_spacers\tboundary_conflict\tresource_limited\taccepted\treason\tcore"
        "\tpair_only\tobserved_right_repeats\tsource_candidate_ids"
        "\tselection_eligible\terror_ambiguous\terror_status\topportunity_records"
        "\topportunity_min_spacer_quality\topportunity_complete\terror_hypotheses"
        "\terror_neighbor_count\terror_neighbor_candidate_id\terror_observed_units"
        "\terror_probability_upper\terror_log_bound\terror_log_threshold"
        "\tnon_tandem_sequence_units\tnon_tandem_start_offsets\tnon_tandem_units_saturated\tnon_tandem_offsets_saturated"
        "\ttandem_ambiguous\ttandem_status\tselection_status\ttandem_explained\ttandem_model"
        "\ttandem_training_kind\ttandem_training_side\ttandem_model_orientation\ttandem_motif\ttandem_period\ttandem_spacer_template_length"
        "\ttandem_left_observed_bases\ttandem_right_observed_bases\ttandem_repeat_similarity"
        "\ttandem_left_similarity\ttandem_spacer_similarity\ttandem_right_similarity\ttandem_whole_similarity"
        "\ttandem_repeat_spacer_similarity\ttandem_repeat_max_base_fraction\ttandem_repeat_max_3mer_fraction\n";
    for (std::size_t i = 0; i < result.candidates.size(); ++i) {
        const auto& c = result.candidates[i];
        candidates << i + 1 << '\t';
        bool first = true;
        for (const auto owner : c.family_owners) {
            if (!first) candidates << ';';
            candidates << "Array_" << owner + 1; first = false;
        }
        candidates << '\t' << c.left_repeat << '\t' << c.spacer << '\t' << c.right_repeat
            << '\t' << c.evidence_kind << '\t' << c.graph_proposed << '\t' << c.left_complete
            << '\t' << c.right_complete << '\t' << c.left_observed_bases << '\t' << c.right_observed_bases
            << '\t' << c.support.sequence_units << '\t' << c.support.start_offsets
            << '\t' << c.support.units_saturated << '\t' << c.support.offsets_saturated
            << '\t' << c.distinct_supported_family_spacers << '\t' << c.boundary_conflict
            << '\t' << c.resource_limited << '\t' << c.accepted << '\t' << c.reason << '\t' << c.core
            << '\t' << c.pair_only << '\t';
        for (std::size_t j = 0; j < c.observed_right_repeats.size(); ++j) {
            if (j) candidates << ';';
            candidates << c.observed_right_repeats[j];
        }
        candidates << '\t';
        for (std::size_t j = 0; j < c.source_candidate_indices.size(); ++j) {
            if (j) candidates << ';';
            candidates << c.source_candidate_indices[j] + 1;
        }
        candidates << '\t' << c.selection_eligible << '\t' << c.error_ambiguous << '\t' << c.error_status
            << '\t' << c.opportunity_records << '\t' << c.opportunity_min_spacer_quality
            << '\t' << c.opportunity_complete << '\t' << c.error_hypotheses
            << '\t' << c.error_neighbor_count << '\t'
            << (c.error_neighbor_index == std::numeric_limits<std::size_t>::max() ? 0 : c.error_neighbor_index + 1)
            << '\t' << c.error_observed_units << '\t' << c.error_probability_upper
            << '\t' << c.error_log_bound << '\t' << c.error_log_threshold
            << '\t' << c.non_tandem_support.sequence_units << '\t' << c.non_tandem_support.start_offsets
            << '\t' << c.non_tandem_support.units_saturated << '\t' << c.non_tandem_support.offsets_saturated
            << '\t' << c.tandem_ambiguous << '\t' << c.tandem_status << '\t' << selection_status(c)
            << '\t' << c.tandem.explained << '\t' << c.tandem.model << '\t' << c.tandem.training_kind
            << '\t' << c.tandem.training_side << '\t' << c.tandem.model_orientation
            << '\t' << c.tandem.motif << '\t' << c.tandem.period
            << '\t' << c.tandem.spacer_template_length << '\t' << c.tandem.left_observed_bases
            << '\t' << c.tandem.right_observed_bases << '\t' << c.tandem.repeat_similarity
            << '\t' << c.tandem.left_similarity << '\t' << c.tandem.spacer_similarity
            << '\t' << c.tandem.right_similarity << '\t' << c.tandem.whole_similarity
            << '\t' << c.tandem.repeat_spacer_similarity << '\t' << c.tandem.repeat_max_base_fraction
            << '\t' << c.tandem.repeat_max_3mer_fraction << '\n';
    }
    out << "array_id\tentry_index\trepeat_variant\tspacer\tcandidate_id\tevidence_kind"
           "\tsequence_units\tstart_offsets\tadded_row\tselection_eligible\terror_ambiguous\terror_status"
           "\ttandem_ambiguous\ttandem_status\tselection_status\tnon_tandem_sequence_units\tnon_tandem_start_offsets"
           "\tcandidate_evidence_eligible\tcandidate_evidence_status\trow_origin\tcontext_status"
           "\tcontext_owner_count\tcontext_supported_owners\n";
    for (const auto& [key, id] : proofs) {
        const auto owner = key.first;
        const auto& pair = key.second;
        const auto& c = result.candidates[id];
        const auto found = original[owner].find(pair);
        const bool added = found == original[owner].end();
        const auto& context = contexts.at(id);
        const auto context_count = supported_owners(context);
        const auto status = row_selection_status(c, added, context);
        const bool eligible = status == "supported";
        std::size_t index;
        if (added) {
            arrays[owner].entries.push_back(pair);
            index = arrays[owner].entries.size();
            ++counts["read_recruitment_added_rows"];
        } else index = found->second;
        evidence.emplace(std::make_pair(owner, index), supported_calls::RecruitmentSupport{
            id + 1, c.support.sequence_units, c.support.start_offsets, c.evidence_kind, eligible, added, status});
        out << "Array_" << owner + 1 << '\t' << index << '\t' << pair.first << '\t' << pair.second
            << '\t' << id + 1 << '\t' << c.evidence_kind << '\t' << c.support.sequence_units
            << '\t' << c.support.start_offsets << '\t' << added << '\t' << eligible
            << '\t' << c.error_ambiguous << '\t' << c.error_status << '\t' << c.tandem_ambiguous
            << '\t' << c.tandem_status << '\t' << status << '\t' << c.non_tandem_support.sequence_units
            << '\t' << c.non_tandem_support.start_offsets << '\t' << c.selection_eligible
            << '\t' << selection_status(c) << '\t' << (added ? "recruited_observation" : "original_graph")
            << '\t' << (context_count == context.size() && !context.empty() ?
                "all_owners_supported" : "insufficient_spacer_diversity")
            << '\t' << context.size() << '\t' << context_count << '\n';
        counts["read_recruitment_policy_eligible_rows"] += eligible;
        counts["read_recruitment_novel_boundary_unresolved_rows"] += status == "novel_boundary_unresolved";
        counts["read_recruitment_insufficient_spacer_diversity_rows"] += status == "insufficient_spacer_diversity";
    }
    counts["read_recruitment_bound_rows"] = evidence.size();
    supported_calls::close_output(candidates, candidate_path);
    supported_calls::close_output(out, proof_path);
    supported_calls::close_output(snapshot, original_path);
    supported_calls::close_output(context_out, context_path);
    return evidence;
}

}} // namespace Zuck::read_recruitment_output
#endif
