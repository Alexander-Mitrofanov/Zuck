#ifndef ZUCK_SUPPORTED_CALLS_H
#define ZUCK_SUPPORTED_CALLS_H

#include <algorithm>
#include <array>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <limits>
#include <map>
#include <set>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>
#include "crispr_array/element_confidence.h"

namespace Zuck { namespace supported_calls {

// A reporting profile over existing oriented rows. It never deletes candidates,
// changes their grouping, resolves a biological boundary, or rescans reads.
constexpr std::size_t STANDARD_MINIMUM_SPACER = 26;
constexpr std::size_t STANDARD_MAXIMUM_SPACER = 50;

// A separately measured channel; it never changes literal-context counts.
// The integration layer binds each proof to its original array/entry identity.
struct RecruitmentSupport {
    std::size_t candidate_id = 0, sequence_units = 0, start_offsets = 0;
    std::string evidence_kind;
    bool selection_eligible = true, added_row = false;
    std::string selection_block_reason = "error_ambiguous";
    bool observed() const {
        return candidate_id && sequence_units >= 2 && start_offsets >= 2 && !(added_row && selection_eligible) &&
            (evidence_kind == "full_pair" || evidence_kind == "full_partial" ||
             evidence_kind == "graph_junction") && (selection_eligible ||
             selection_block_reason == "error_ambiguous" || selection_block_reason == "tandem_ambiguous" ||
             selection_block_reason == "novel_boundary_unresolved" ||
             selection_block_reason == "insufficient_spacer_diversity");
    }
    bool accepted() const { return observed() && selection_eligible; }
    std::string status() const {
        return accepted() ? "supported" : observed() ? selection_block_reason : "absent";
    }
};
using RecruitmentEvidence = std::map<std::pair<std::size_t, std::size_t>, RecruitmentSupport>;

struct Decision {
    std::size_t array_index = 0; // Zero-based index in the complete array vector.
    std::size_t entry_index = 0; // One-based index in that original array.
    std::string repeat_variant, spacer, canonical_spacer;
    bool selected = false, length_in_profile = true;
    std::string context_status, context_observability, reason;
    element_confidence::Support context;
    RecruitmentSupport recruitment;
    std::string selected_evidence;
};

struct Group {
    std::size_t array_index = 0;
    std::vector<std::size_t> entry_indices; // Original one-based row indexes.
    std::set<std::string> canonical_spacers;
};
struct RepeatHaplotypeGroup {
    std::size_t array_index = 0;
    std::string repeat_haplotype;
    std::vector<std::size_t> entry_indices;
    std::set<std::string> canonical_spacers;
};
struct CircularUnitDecision {
    std::string circular_identity, canonical_spacer, winner_spacer;
    bool winner = false, retained = false;
    std::size_t selected_occurrences = 0, literal_owners = 0,
                supported_owners = 0, sum_owner_balanced_support = 0,
                min_owner_balanced_support = 0, max_owner_balanced_support = 0;
    std::vector<std::pair<std::size_t, std::size_t>> entries;
};

struct Result {
    std::string profile;
    bool recruitment_enabled = false;
    std::vector<Decision> rows;
    std::vector<Group> groups;
    std::vector<RepeatHaplotypeGroup> repeat_haplotype_groups;
    std::vector<CircularUnitDecision> circular_unit_decisions;
    std::set<std::string> spacers; // Any selected row licenses only this inventory.
    std::map<std::string, std::size_t> counts;
};

inline void validate_profile(const std::string& profile) {
    if (profile != "standard" && profile != "context" && profile != "off")
        throw std::invalid_argument("call-profile must be standard, context, or off");
}

inline bool length_allowed(const std::string& profile, std::size_t length) {
    return profile != "standard" ||
        (length >= STANDARD_MINIMUM_SPACER && length <= STANDARD_MAXIMUM_SPACER);
}

inline void add_reason(std::string& reasons, const std::string& reason) {
    if (!reasons.empty()) reasons += ';';
    reasons += reason;
}

template <class Arrays>
void rebuild_selected_view(const Arrays& arrays, Result& result) {
    result.groups.clear();
    result.repeat_haplotype_groups.clear();
    result.spacers.clear();
    result.counts["supported_selected_rows"] = 0;
    result.counts["supported_unselected_rows"] = 0;
    result.counts["supported_family_groups"] = 0;
    result.counts["supported_singleton_groups"] = 0;
    result.counts["supported_multiple_spacer_groups"] = 0;
    result.counts["supported_repeat_haplotype_memberships"] = 0;
    if (result.recruitment_enabled)
        result.counts["supported_recruitment_selected_rows"] = 0;
    for (std::size_t array = 0; array < arrays.size(); ++array) {
        Group group;
        group.array_index = array;
        std::map<std::string, RepeatHaplotypeGroup> repeat_groups;
        for (auto& row : result.rows) {
            if (row.array_index != array) continue;
            if (!row.selected) {
                ++result.counts["supported_unselected_rows"];
                continue;
            }
            ++result.counts["supported_selected_rows"];
            if (result.recruitment_enabled)
                result.counts["supported_recruitment_selected_rows"] +=
                    row.selected_evidence == "read_recruitment";
            group.entry_indices.push_back(row.entry_index);
            group.canonical_spacers.insert(row.canonical_spacer);
            result.spacers.insert(row.canonical_spacer);
            const auto haplotype = array_database::canonical(row.repeat_variant);
            auto& repeat_group = repeat_groups[haplotype];
            repeat_group.array_index = array;
            repeat_group.repeat_haplotype = haplotype;
            repeat_group.entry_indices.push_back(row.entry_index);
            repeat_group.canonical_spacers.insert(row.canonical_spacer);
        }
        if (group.entry_indices.empty()) continue;
        for (auto& repeat_group : repeat_groups) {
            auto& candidate = repeat_group.second;
            if (candidate.canonical_spacers.size() < 2 ||
                candidate.canonical_spacers == group.canonical_spacers) continue;
            result.counts["supported_repeat_haplotype_memberships"] +=
                candidate.entry_indices.size();
            result.repeat_haplotype_groups.push_back(std::move(candidate));
        }
        ++result.counts["supported_family_groups"];
        ++result.counts[group.canonical_spacers.size() == 1 ?
            "supported_singleton_groups" : "supported_multiple_spacer_groups"];
        result.groups.push_back(std::move(group));
    }
    result.counts["supported_spacer_types"] = result.spacers.size();
    result.counts["supported_repeat_haplotype_groups"] =
        result.repeat_haplotype_groups.size();
}

// Confidence rows must have been annotated in the same flattened array/entry
// order. Confidence itself carries no sequence IDs; count checks cannot detect
// a caller supplying an unrelated, equally sized confidence vector.
template <class Arrays>
Result infer(const Arrays& arrays, const element_confidence::Result& confidence,
             const std::string& profile, const RecruitmentEvidence& recruitment = {},
             bool recruitment_enabled = false,
             bool singleton_precision_prune = false,
             bool circular_unit_dedup = false) {
    validate_profile(profile);
    Result result;
    result.profile = profile;
    result.recruitment_enabled = recruitment_enabled;
    if (!recruitment_enabled && !recruitment.empty())
        throw std::invalid_argument("Recruitment evidence supplied while disabled");
    result.counts = {{"supported_call_profile_enabled", profile != "off"},
        {"supported_input_rows", 0}, {"supported_selected_rows", 0},
        {"supported_unselected_rows", 0}, {"supported_unknown_context_rows", 0},
        {"supported_insufficient_context_rows", 0}, {"supported_out_of_profile_rows", 0},
        {"supported_family_groups", 0}, {"supported_singleton_groups", 0},
        {"supported_multiple_spacer_groups", 0}, {"supported_spacer_types", 0}};
    result.counts["supported_repeat_haplotype_groups"] = 0;
    result.counts["supported_repeat_haplotype_memberships"] = 0;
    result.counts["supported_singleton_precision_prune_enabled"] =
        singleton_precision_prune;
    result.counts["supported_singleton_precision_pruned_rows"] = 0;
    result.counts["supported_singleton_precision_pruned_types"] = 0;
    result.counts["supported_circular_unit_dedup_enabled"] = circular_unit_dedup;
    result.counts["supported_circular_unit_duplicate_units"] = 0;
    result.counts["supported_circular_unit_pruned_rows"] = 0;
    result.counts["supported_circular_unit_pruned_types"] = 0;
    result.counts["supported_circular_unit_protected_alternatives"] = 0;
    if (profile == "off") return result;
    std::size_t total = 0;
    for (const auto& array : arrays) total += array.entries.size();
    if (total != confidence.rows.size())
        throw std::invalid_argument("Supported-call confidence row count mismatch");
    for (const auto& item : recruitment)
        if (item.first.first >= arrays.size() || item.first.second == 0 ||
            item.first.second > arrays[item.first.first].entries.size())
            throw std::invalid_argument("Recruitment proof references an absent row");
    result.counts["supported_input_rows"] = total;
    if (recruitment_enabled) result.counts["supported_recruitment_selected_rows"] = 0;
    result.rows.reserve(total);
    std::size_t offset = 0;
    for (std::size_t a = 0; a < arrays.size(); ++a) {
        const auto& array = arrays[a];
        if (!array_database::dna(array.consensus))
            throw std::invalid_argument("Supported-call consensus must contain only ACGT");
        Group group;
        group.array_index = a;
        std::map<std::string, RepeatHaplotypeGroup> repeat_groups;
        for (std::size_t e = 0; e < array.entries.size(); ++e) {
            const auto& pair = array.entries[e];
            if (!array_database::dna(pair.first) || !array_database::dna(pair.second))
                throw std::invalid_argument("Supported-call pairs must contain only ACGT");
            const auto& evidence = confidence.rows[offset++];
            const bool known = evidence.context_observability == "observable";
            const bool missing = evidence.context_observability == "no_reads" ||
                evidence.context_observability == "no_quality" ||
                evidence.context_observability == "span_exceeds_max_quality_read";
            if ((!known && !missing) || known != evidence.context.observable)
                throw std::invalid_argument("Inconsistent supported-call context observability");
            Decision row;
            row.array_index = a;
            row.entry_index = e + 1;
            row.repeat_variant = pair.first;
            row.spacer = pair.second;
            row.canonical_spacer = array_database::canonical(pair.second);
            row.context = evidence.context;
            row.context_observability = evidence.context_observability;
            row.length_in_profile = length_allowed(profile, pair.second.size());
            row.context_status = !known ? "unknown" :
                evidence.context.accepted() ? "supported" : "insufficient";
            const auto recruited = recruitment.find({a, e + 1});
            if (recruited != recruitment.end()) {
                if (!recruited->second.observed())
                    throw std::invalid_argument("Invalid recruited-call support");
                row.recruitment = recruited->second;
            }
            // Original rows keep their original literal-context selection.
            // A newly appended uncertain allele cannot bypass its recruitment
            // gate merely because the downstream context pass also sees it.
            const bool literal_eligible = evidence.context.accepted() &&
                !row.recruitment.added_row;
            row.selected = row.length_in_profile &&
                (literal_eligible || row.recruitment.accepted());
            row.selected_evidence = !row.selected ? "none" : literal_eligible
                ? "literal_context" : "read_recruitment";
            if (row.recruitment.observed() && !row.recruitment.selection_eligible)
                add_reason(row.reason, "recruitment_" + row.recruitment.selection_block_reason);
            if (!row.length_in_profile) {
                add_reason(row.reason, "spacer_length_outside_profile");
                ++result.counts["supported_out_of_profile_rows"];
            }
            if (!known) {
                add_reason(row.reason, evidence.context_observability);
                ++result.counts["supported_unknown_context_rows"];
            } else if (!evidence.context.accepted()) {
                if (evidence.context.sequence_units < 2)
                    add_reason(row.reason, "insufficient_context_units");
                if (evidence.context.start_offsets < 2)
                    add_reason(row.reason, "insufficient_context_offsets");
                ++result.counts["supported_insufficient_context_rows"];
            }
            if (row.selected) {
                row.reason = row.selected_evidence == "read_recruitment"
                    ? "passes_read_recruitment" : "passes_profile";
                if (row.selected_evidence == "read_recruitment")
                    ++result.counts["supported_recruitment_selected_rows"];
                ++result.counts["supported_selected_rows"];
                group.entry_indices.push_back(row.entry_index);
                group.canonical_spacers.insert(row.canonical_spacer);
                result.spacers.insert(row.canonical_spacer);
                const auto haplotype = array_database::canonical(pair.first);
                auto& repeat_group = repeat_groups[haplotype];
                repeat_group.array_index = a;
                repeat_group.repeat_haplotype = haplotype;
                repeat_group.entry_indices.push_back(row.entry_index);
                repeat_group.canonical_spacers.insert(row.canonical_spacer);
            } else {
                ++result.counts["supported_unselected_rows"];
            }
            result.rows.push_back(std::move(row));
        }
        if (!group.entry_indices.empty()) {
            for (auto& repeat_group : repeat_groups) {
                auto& candidate = repeat_group.second;
                if (candidate.canonical_spacers.size() < 2 ||
                    candidate.canonical_spacers == group.canonical_spacers) continue;
                result.counts["supported_repeat_haplotype_memberships"] +=
                    candidate.entry_indices.size();
                result.repeat_haplotype_groups.push_back(std::move(candidate));
            }
            ++result.counts["supported_family_groups"];
            ++result.counts[group.canonical_spacers.size() == 1 ?
                "supported_singleton_groups" : "supported_multiple_spacer_groups"];
            result.groups.push_back(std::move(group));
        }
    }
    result.counts["supported_spacer_types"] = result.spacers.size();
    result.counts["supported_repeat_haplotype_groups"] =
        result.repeat_haplotype_groups.size();
    if (singleton_precision_prune) {
        std::map<std::size_t, std::size_t> family_spacer_types;
        for (const auto& group : result.groups)
            family_spacer_types[group.array_index] = group.canonical_spacers.size();
        std::map<std::string, std::vector<std::size_t>> occurrences;
        for (std::size_t index = 0; index < result.rows.size(); ++index)
            if (result.rows[index].selected)
                occurrences[result.rows[index].canonical_spacer].push_back(index);
        std::set<std::string> pruned;
        for (const auto& item : occurrences) {
            bool eligible = true;
            for (auto index : item.second) {
                const auto& row = result.rows[index];
                eligible &= family_spacer_types[row.array_index] == 1 &&
                    row.selected_evidence == "literal_context" &&
                    row.context.sequence_units <= 2 && row.context.start_offsets <= 2;
            }
            if (eligible) pruned.insert(item.first);
        }
        for (auto& row : result.rows) if (row.selected && pruned.count(row.canonical_spacer)) {
            row.selected = false;
            row.selected_evidence = "none";
            row.reason = "singleton_minimal_literal_context_pruned";
            ++result.counts["supported_singleton_precision_pruned_rows"];
        }
        result.counts["supported_singleton_precision_pruned_types"] = pruned.size();

    }
    if (circular_unit_dedup) {
        struct OwnerEvidence { bool literal = false; std::size_t balanced = 0; };
        struct TypeEvidence {
            std::vector<std::size_t> rows;
            std::map<std::size_t, OwnerEvidence> owners;
        };
        std::map<std::string, std::map<std::string, TypeEvidence>> units;
        for (std::size_t index = 0; index < result.rows.size(); ++index) {
            const auto& row = result.rows[index];
            if (!row.selected) continue;
            auto& evidence = units[element_confidence::circular_identity(
                row.repeat_variant + row.spacer)][row.canonical_spacer];
            evidence.rows.push_back(index);
            auto& owner = evidence.owners[row.array_index];
            const bool literal = row.selected_evidence == "literal_context";
            const auto balanced = literal
                ? std::min(row.context.sequence_units, row.context.start_offsets)
                : std::min(row.recruitment.sequence_units,
                           row.recruitment.start_offsets);
            if (literal && !owner.literal) {
                owner.literal = true;
                owner.balanced = balanced;
            } else if (literal == owner.literal) {
                owner.balanced = std::max(owner.balanced, balanced);
            }
        }
        const auto better = [](const auto& left, const auto& right) {
            const auto& a = left.second;
            const auto& b = right.second;
            const auto summarize = [](const TypeEvidence& value) {
                std::size_t literal = 0, sum = 0,
                            minimum = std::numeric_limits<std::size_t>::max(), maximum = 0;
                for (const auto& owner : value.owners) {
                    literal += owner.second.literal;
                    sum += owner.second.balanced;
                    minimum = std::min(minimum, owner.second.balanced);
                    maximum = std::max(maximum, owner.second.balanced);
                }
                if (value.owners.empty()) minimum = 0;
                return std::array<std::size_t, 6>{literal, value.owners.size(), sum,
                                                  minimum, maximum, value.rows.size()};
            };
            const auto ar = summarize(a), br = summarize(b);
            if (ar != br) return ar > br;
            return left.first < right.first;
        };
        std::set<std::string> winners, selected_types;
        for (const auto& unit : units) {
            auto winner = unit.second.begin();
            for (auto candidate = std::next(winner); candidate != unit.second.end(); ++candidate)
                if (better(*candidate, *winner)) winner = candidate;
            winners.insert(winner->first);
            for (const auto& item : unit.second) selected_types.insert(item.first);
            if (unit.second.size() > 1)
                ++result.counts["supported_circular_unit_duplicate_units"];
        }
        std::set<std::string> pruned;
        std::set_difference(selected_types.begin(), selected_types.end(),
                            winners.begin(), winners.end(),
                            std::inserter(pruned, pruned.end()));
        for (const auto& unit : units) {
            if (unit.second.size() < 2) continue;
            auto winner = unit.second.begin();
            for (auto candidate = std::next(winner); candidate != unit.second.end(); ++candidate)
                if (better(*candidate, *winner)) winner = candidate;
            for (const auto& item : unit.second) {
                CircularUnitDecision decision;
                decision.circular_identity = unit.first;
                decision.canonical_spacer = item.first;
                decision.winner_spacer = winner->first;
                decision.winner = item.first == winner->first;
                decision.retained = !pruned.count(item.first);
                decision.selected_occurrences = item.second.rows.size();
                decision.supported_owners = item.second.owners.size();
                decision.min_owner_balanced_support =
                    std::numeric_limits<std::size_t>::max();
                for (const auto& owner : item.second.owners) {
                    decision.literal_owners += owner.second.literal;
                    decision.sum_owner_balanced_support += owner.second.balanced;
                    decision.min_owner_balanced_support = std::min(
                        decision.min_owner_balanced_support, owner.second.balanced);
                    decision.max_owner_balanced_support = std::max(
                        decision.max_owner_balanced_support, owner.second.balanced);
                }
                if (item.second.owners.empty())
                    decision.min_owner_balanced_support = 0;
                for (const auto index : item.second.rows)
                    decision.entries.push_back({result.rows[index].array_index,
                                                result.rows[index].entry_index});
                if (!decision.winner && decision.retained)
                    ++result.counts["supported_circular_unit_protected_alternatives"];
                result.circular_unit_decisions.push_back(std::move(decision));
            }
        }
        for (auto& row : result.rows) if (row.selected && pruned.count(row.canonical_spacer)) {
            row.selected = false;
            row.selected_evidence = "none";
            row.reason = "circular_unit_alternative_partition_pruned";
            ++result.counts["supported_circular_unit_pruned_rows"];
        }
        result.counts["supported_circular_unit_pruned_types"] = pruned.size();
    }
    if (singleton_precision_prune || circular_unit_dedup)
        rebuild_selected_view(arrays, result);
    return result;
}

inline void close_output(std::ofstream& out, const std::filesystem::path& path) {
    out.close();
    if (!out) throw std::runtime_error("Error writing supported calls: " + path.string());
}

// The enclosing result-bundle transaction owns all standard output files and
// optional precision-demotion ledgers, including
// stale-file removal when the profile is disabled. Array and row IDs are never
// renumbered. Singleton families remain explicit sequence groups, not proven loci.
template <class Arrays>
void write_tables(const std::string& directory, const Arrays& arrays, const Result& result) {
    validate_profile(result.profile);
    if (result.profile == "off") return;
    std::size_t row_number = 0;
    std::vector<Group> groups;
    std::set<std::string> spacers;
    std::set<std::pair<std::size_t, std::size_t>> selected_rows;
    for (std::size_t a = 0; a < arrays.size(); ++a) {
        Group group;
        group.array_index = a;
        for (std::size_t e = 0; e < arrays[a].entries.size(); ++e) {
            if (row_number >= result.rows.size())
                throw std::invalid_argument("Supported-call output row count mismatch");
            const auto& row = result.rows[row_number++];
            if (row.array_index != a || row.entry_index != e + 1 ||
                row.repeat_variant != arrays[a].entries[e].first ||
                row.spacer != arrays[a].entries[e].second ||
                row.canonical_spacer != array_database::canonical(row.spacer))
                throw std::invalid_argument("Supported-call output row identity mismatch");
            if (row.selected) {
                selected_rows.insert({a, row.entry_index});
                group.entry_indices.push_back(row.entry_index);
                group.canonical_spacers.insert(row.canonical_spacer);
                spacers.insert(row.canonical_spacer);
            }
        }
        if (!group.entry_indices.empty()) groups.push_back(std::move(group));
    }
    if (row_number != result.rows.size() || spacers != result.spacers || groups.size() != result.groups.size())
        throw std::invalid_argument("Supported-call output inventory mismatch");
    for (std::size_t i = 0; i < groups.size(); ++i)
        if (groups[i].array_index != result.groups[i].array_index ||
            groups[i].entry_indices != result.groups[i].entry_indices ||
            groups[i].canonical_spacers != result.groups[i].canonical_spacers)
            throw std::invalid_argument("Supported-call output group membership mismatch");
    for (const auto& group : result.repeat_haplotype_groups) {
        if (group.array_index >= arrays.size() || group.canonical_spacers.size() < 2)
            throw std::invalid_argument("Invalid repeat-haplotype group");
        std::set<std::size_t> unique_entries;
        std::set<std::string> observed_spacers;
        for (const auto entry : group.entry_indices) {
            if (!entry || entry > arrays[group.array_index].entries.size() ||
                !unique_entries.insert(entry).second ||
                !selected_rows.count({group.array_index, entry}) ||
                array_database::canonical(arrays[group.array_index].entries[entry - 1].first) !=
                    group.repeat_haplotype)
                throw std::invalid_argument("Repeat-haplotype membership differs from emitted row");
            observed_spacers.insert(
                array_database::canonical(arrays[group.array_index].entries[entry - 1].second));
        }
        if (observed_spacers != group.canonical_spacers)
            throw std::invalid_argument("Repeat-haplotype spacer inventory differs");
    }

    const auto root = std::filesystem::path(directory);
    const auto decisions_path = root / "call_decisions.tsv";
    const auto entries_path = root / "supported_entries.tsv";
    const auto arrays_path = root / "supported_arrays.tsv";
    const auto fasta_path = root / "supported_spacers.fasta";
    const auto haplotypes_path = root / "repeat_haplotype_groups.tsv";
    std::ofstream decisions(decisions_path), entries(entries_path), array_out(arrays_path),
        fasta(fasta_path), haplotypes(haplotypes_path);
    if (!decisions || !entries || !array_out || !fasta || !haplotypes)
        throw std::runtime_error("Cannot write supported-call output files: " + directory);
    decisions << "array_id\tentry_index\trepeat_variant\tspacer\tcanonical_spacer\tprofile\tselected"
        "\tcontext_status\tcontext_observability\tcontext_sequence_units\tcontext_start_offsets"
        "\tcontext_units_saturated\tcontext_offsets_saturated\tlength_in_profile\treason";
    if (result.recruitment_enabled)
        decisions << "\tselected_evidence\trecruitment_status\trecruitment_kind\trecruitment_candidate_id"
                     "\trecruitment_sequence_units\trecruitment_start_offsets";
    decisions << '\n';
    entries << "array_id\tentry_index\trepeat_variant\tspacer\tcanonical_spacer\n";
    array_out << "array_id\trepeat_consensus\tentry_count\tdistinct_canonical_spacer_count\tgroup_status\tprofile\n";
    haplotypes << "group_id\tparent_array_id\trepeat_haplotype\tentry_index\n";
    for (const auto& row : result.rows) {
        decisions << "Array_" << row.array_index + 1 << '\t' << row.entry_index << '\t'
            << row.repeat_variant << '\t' << row.spacer << '\t' << row.canonical_spacer << '\t'
            << result.profile << '\t' << row.selected << '\t' << row.context_status << '\t'
            << row.context_observability << '\t' << row.context.sequence_units << '\t'
            << row.context.start_offsets << '\t' << row.context.units_saturated << '\t'
            << row.context.offsets_saturated << '\t' << row.length_in_profile << '\t' << row.reason;
        if (result.recruitment_enabled)
            decisions << '\t' << row.selected_evidence << '\t'
                << row.recruitment.status() << '\t'
                << row.recruitment.evidence_kind << '\t' << row.recruitment.candidate_id << '\t'
                << row.recruitment.sequence_units << '\t' << row.recruitment.start_offsets;
        decisions << '\n';
        if (row.selected)
            entries << "Array_" << row.array_index + 1 << '\t' << row.entry_index << '\t'
                << row.repeat_variant << '\t' << row.spacer << '\t' << row.canonical_spacer << '\n';
    }
    for (const auto& group : result.groups)
        array_out << "Array_" << group.array_index + 1 << '\t' << arrays[group.array_index].consensus << '\t'
            << group.entry_indices.size() << '\t' << group.canonical_spacers.size() << '\t'
            << (group.canonical_spacers.size() == 1 ? "singleton" : "multiple_spacers") << '\t'
            << result.profile << '\n';
    for (std::size_t i = 0; i < result.repeat_haplotype_groups.size(); ++i) {
        const auto& group = result.repeat_haplotype_groups[i];
        for (const auto entry : group.entry_indices) {
            haplotypes << "RepeatHaplotype_" << i + 1 << "\tArray_" << group.array_index + 1
                << '\t' << group.repeat_haplotype << '\t' << entry << '\n';
        }
    }
    std::size_t number = 0;
    for (const auto& spacer : result.spacers)
        fasta << ">Zuck_supported_" << ++number << '\n' << spacer << '\n';
    close_output(decisions, decisions_path);
    close_output(entries, entries_path);
    close_output(array_out, arrays_path);
    close_output(fasta, fasta_path);
    close_output(haplotypes, haplotypes_path);
    if (result.counts.at("supported_singleton_precision_prune_enabled")) {
        const auto pruned_path = root / "singleton_precision_pruned.tsv";
        std::ofstream pruned(pruned_path);
        if (!pruned) throw std::runtime_error(
            "Cannot write singleton precision ledger: " + pruned_path.string());
        pruned << "array_id\tentry_index\trepeat_variant\tspacer\tcanonical_spacer"
                  "\tcontext_sequence_units\tcontext_start_offsets\treason\n";
        for (const auto& row : result.rows)
            if (row.reason == "singleton_minimal_literal_context_pruned")
                pruned << "Array_" << row.array_index + 1 << '\t' << row.entry_index << '\t'
                    << row.repeat_variant << '\t' << row.spacer << '\t' << row.canonical_spacer
                    << '\t' << row.context.sequence_units << '\t' << row.context.start_offsets
                    << '\t' << row.reason << '\n';
        close_output(pruned, pruned_path);
    }
    if (result.counts.at("supported_circular_unit_dedup_enabled")) {
        const auto ledger_path = root / "circular_unit_dedup.tsv";
        std::ofstream ledger(ledger_path);
        if (!ledger) throw std::runtime_error(
            "Cannot write circular-unit dedup ledger: " + ledger_path.string());
        ledger << "circular_identity\tcanonical_spacer\twinner_spacer\twinner\tretained"
                  "\tselected_occurrences\tliteral_owners\tsupported_owners"
                  "\tsum_owner_balanced_support\tmin_owner_balanced_support"
                  "\tmax_owner_balanced_support"
                  "\tarray_entries\treason\n";
        for (const auto& item : result.circular_unit_decisions) {
            ledger << item.circular_identity << '\t' << item.canonical_spacer << '\t'
                   << item.winner_spacer << '\t' << item.winner << '\t' << item.retained
                   << '\t' << item.selected_occurrences << '\t' << item.literal_owners
                   << '\t' << item.supported_owners << '\t'
                   << item.sum_owner_balanced_support << '\t'
                   << item.min_owner_balanced_support << '\t'
                   << item.max_owner_balanced_support << '\t';
            for (std::size_t i = 0; i < item.entries.size(); ++i) {
                if (i) ledger << ',';
                ledger << "Array_" << item.entries[i].first + 1 << ':' << item.entries[i].second;
            }
            ledger << '\t' << (item.winner ? "evidence_rank_winner" : item.retained
                ? "protected_winner_in_other_unit" : "alternative_partition_pruned") << '\n';
        }
        close_output(ledger, ledger_path);
    }
}

}} // namespace Zuck::supported_calls
#endif
