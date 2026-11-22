#ifndef ZUCK_READ_FIRST_OUTPUT_H
#define ZUCK_READ_FIRST_OUTPUT_H

#include "core/settings.h"
#include "crispr_array/read_first_overlap.h"
#include "crispr_array/read_first_variants.h"
#include "crispr_array/read_first_family_recut.h"
#include "crispr_array/supported_calls.h"
#include <fstream>

namespace Zuck { namespace read_first {

// Require two substantially different, quality-qualified spacer interiors before
// using a read family to add any types. Repeated copies of a genomic tandem
// unit can otherwise supply many distinct templates for one homologous tract.
inline std::set<std::string> diverse_families(const std::vector<Call>& strict) {
    std::map<std::string, std::vector<const std::string*>> first_interiors;
    std::set<std::string> established;
    for (const auto& call : strict) {
        if (established.count(call.repeat)) continue;
        auto& earlier = first_interiors[call.repeat];
        for (const auto* other : earlier) {
            if (4 * read_recruitment::detail::edit_distance(call.spacer, *other) >=
                std::max(call.spacer.size(), other->size())) {
                established.insert(call.repeat);
                break;
            }
        }
        // Calls are in stable family/spacer order. A bounded set of witnesses
        // limits work in a large population of near-identical variants.
        if (!established.count(call.repeat) && earlier.size() < 32) earlier.push_back(&call.spacer);
    }
    return established;
}

// Verify both junctions and the whole observed unit before comparing endpoint evidence.
inline bool terminal_recut(const Call& observed, std::string repeat, std::string spacer) {
    if (repeat.size() <= observed.repeat.size() || repeat.size() - observed.repeat.size() > 2) return false;
    const auto extra = repeat.size() - observed.repeat.size();
    if (spacer.size() + extra != observed.spacer.size()) return false;
    for (unsigned strand = 0; strand < 2; ++strand) {
        for (std::size_t left = 0; left <= extra; ++left) {
            const auto right = extra - left;
            if (repeat.compare(left, observed.repeat.size(), observed.repeat) == 0 &&
                repeat.compare(0, left, observed.spacer, observed.spacer.size() - left, left) == 0 &&
                repeat.compare(left + observed.repeat.size(), right, observed.spacer, 0, right) == 0 &&
                observed.spacer.compare(right, spacer.size(), spacer) == 0) return true;
        }
        repeat = read_linked::reverse_complement(repeat);
        spacer = read_linked::reverse_complement(spacer);
    }
    return false;
}

inline void write_inventory(const std::filesystem::path& directory,
                            const supported_calls::Result& graph_calls,
                            const Settings& settings, std::map<std::string, std::size_t>& counts) {
    const auto paths = settings.effective_grouping_reads();
    const auto families = discover(paths, counts);
    auto observations = recruit(families, paths, counts);
    const auto established = diverse_families(observations);
    counts["read_first_diverse_families"] = established.size();
    const auto within_read = observations.size();
    auto complementary = recruit_overlap(families, paths, observations, counts);
    std::vector<VariantDetail> variant_evidence;
    auto variants = recruit_variants(families, paths, observations, counts, &variant_evidence);
    observations.insert(observations.end(), complementary.begin(), complementary.end());
    const auto through_complementary = observations.size();
    std::set<std::string> existing_spacers, existing_units;
    for (const auto& call : observations) {
        existing_spacers.insert(array_database::canonical(call.spacer));
        existing_units.insert(element_confidence::circular_identity(call.repeat + call.spacer));
    }
    for (const auto& call : variants) {
        if (!existing_spacers.count(array_database::canonical(call.spacer)) &&
            !existing_units.count(element_confidence::circular_identity(call.repeat + call.spacer)))
            observations.push_back(call);
    }
    std::ofstream allele_table(directory / "read_first_variant_evidence.tsv");
    allele_table << "family_repeat\tspacer\tproof\ttemplates\tfull_pairs\texact_templates\tshort_templates\tflank_substitutions\n";
    for (const auto& detail : variant_evidence) {
        allele_table << detail.repeat << '\t' << detail.spacer << '\t' << detail.proof << '\t'
                    << detail.templates << '\t' << detail.full_pairs << '\t' << detail.exact_templates
                    << '\t' << detail.short_templates << '\t';
        bool first_signature = true;
        for (const auto& signature : detail.substitution_signatures) {
            if (!first_signature) allele_table << ',';
            allele_table << signature;
            first_signature = false;
        }
        allele_table << '\n';
    }
    std::ofstream entries(directory / "read_first_entries.tsv");
    std::ofstream fasta(directory / "read_first_spacers.fasta");
    std::ofstream models(directory / "read_first_families.tsv");
    std::ofstream evidence(directory / "read_first_observations.tsv");
    std::ofstream boundaries(directory / "read_first_graph_boundary.tsv");
    entries << "origin\tfamily_repeat\tspacer\tcanonical_spacer\ttemplates\tfull_pairs\tsource_group\tentry_index\n";
    models << "family_repeat\ttemplates\tcontexts\n";
    std::map<std::string, const Family*> by_repeat;
    for (const auto& family : families) {
        models << family.repeat << '\t' << family.templates << '\t' << family.contexts << '\n';
        by_repeat.emplace(family.repeat, &family);
    }
    evidence << "family_repeat\tspacer\tcanonical_spacer\ttemplates\tfull_pairs\tdecision\tgeometry\n";
    boundaries << "graph_entry_id\tgraph_repeat\tgraph_spacer\tcircular_unit\tread_repeat\tread_spacer\tdecision\tother_contexts\n";
    std::map<std::string, std::vector<const Call*>> by_unit;
    for (const auto& call : observations) {
        if (established.count(call.repeat) && settings.dna_sequence_settings.pair_allowed(call.repeat, call.spacer) &&
            supported_calls::length_allowed(graph_calls.profile, call.spacer.size()))
            by_unit[element_confidence::circular_identity(call.repeat + call.spacer)].push_back(&call);
    }
    family_recut::Index family_index;
    std::set<std::string> conflicting_units;
    std::ofstream family_log;
    if (settings.read_first_family_recut) {
        family_log.open(directory / "read_first_family_recut.tsv");
        family_log << "graph_entry_id\tgraph_repeat\tgraph_spacer\tread_repeat\tread_spacer\ttemplates\tfull_pairs"
                     "\tdecision\tother_contexts\tcompeting_repeat\traw_q20_full_templates\traw_support_status\tfocal_support_route\n";
        std::map<std::string,std::pair<std::string,std::string>> original;
        for (const auto& row : graph_calls.rows) {
            if (!row.selected) continue;
            const auto unit = element_confidence::circular_identity(row.repeat_variant + row.spacer);
            const auto model = std::min(std::make_pair(row.repeat_variant,row.spacer),
                                        std::make_pair(read_linked::reverse_complement(row.repeat_variant),
                                                       read_linked::reverse_complement(row.spacer)));
            if (original.emplace(unit, model).first->second != model) conflicting_units.insert(unit);
        }
        // Both candidate branches must explain an existing short-model graph unit.
        // Families observed only in complementary fragments cannot certify a longer cut.
        for (const auto& entry : original) {
            const auto it = by_unit.find(entry.first);
            if (conflicting_units.count(entry.first) || it == by_unit.end() || it->second.size() != 1) continue;
            const auto* call = it->second.front();
            if (static_cast<std::size_t>(call - observations.data()) < within_read)
                family_index.register_call(entry.second.first,entry.second.second,*call,*by_repeat.at(call->repeat));
        }
        family_index.check_focal_templates(paths);
        family_index.check_transitions(paths);
        counts["read_first_family_recut_rows"] = 0;
    }
    std::set<std::string> inventory, graph_units, retained;
    std::set<const Call*> promoted;
    for (const auto& row : graph_calls.rows) {
        if (!row.selected) continue;
        const auto unit = element_confidence::circular_identity(row.repeat_variant + row.spacer);
        graph_units.insert(unit);
        const auto alternatives = by_unit.find(unit);
        const Call* replacement = nullptr;
        const Call* candidate = nullptr;
        family_recut::Decision family_decision;
        std::size_t contexts = 0;
        if (alternatives != by_unit.end() && alternatives->second.size() == 1) {
            const auto& call = *alternatives->second.front();
            candidate = &call;
            const bool within = static_cast<std::size_t>(&call - observations.data()) < within_read;
            if (within && call.templates >= 2 && terminal_recut(call, row.repeat_variant, row.spacer)) {
                contexts = shorter_cut_contexts(*by_repeat.at(call.repeat), row.repeat_variant, call.spacer);
                if (contexts >= 2) replacement = &call;
            }
            if (settings.read_first_family_recut) {
                if (replacement) family_decision = {"prior_shorter_recut",contexts,{}};
                else if (conflicting_units.count(unit)) family_decision.reason = "conflicting_graph_models";
                else if (!within) family_decision.reason = "requires_within_read";
                else {
                    family_decision = family_index.decide(row.repeat_variant,row.spacer,call);
                    if (family_decision.other_contexts >= 2) {
                        replacement = &call;
                        contexts = family_decision.other_contexts;
                        ++counts["read_first_family_recut_rows"];
                    }
                }
            }
        } else if (alternatives != by_unit.end()) {
            family_decision.reason = "ambiguous_read_models";
        }
        if (settings.read_first_family_recut) {
            family_log << "Array_" << row.array_index+1 << ':' << row.entry_index << '\t'
                       << row.repeat_variant << '\t' << row.spacer << '\t';
            if (candidate) family_log << candidate->repeat << '\t' << candidate->spacer << '\t'
                                      << candidate->templates << '\t' << candidate->full_pairs;
            else family_log << "\t\t\t";
            family_log << '\t' << family_decision.reason << '\t' << family_decision.other_contexts
                       << '\t' << family_decision.competing_repeat;
            family_index.log_focal(family_log,candidate);
            family_log << '\n';
        }
        boundaries << "Array_" << row.array_index + 1 << ':' << row.entry_index << '\t'
                   << row.repeat_variant << '\t' << row.spacer << '\t' << unit << '\t';
        if (replacement) {
            promoted.insert(replacement);
            ++counts["read_first_graph_replaced_rows"];
            boundaries << replacement->repeat << '\t' << replacement->spacer << "\treplaced\t" << contexts << '\n';
        } else {
            retained.insert(row.canonical_spacer);
            entries << "graph\t" << row.repeat_variant << '\t' << row.spacer << '\t'
                    << row.canonical_spacer << "\t\t\tArray_" << row.array_index + 1
                    << '\t' << row.entry_index << '\n';
            boundaries << "\t\tretained\t0\n";
        }
    }
    inventory = retained;
    for (std::size_t index = 0; index < observations.size(); ++index) {
        const auto& call = observations[index];
        const auto canonical = array_database::canonical(call.spacer);
        const auto unit = element_confidence::circular_identity(call.repeat + call.spacer);
        const char* decision = "selected";
        if (!settings.dna_sequence_settings.pair_allowed(call.repeat, call.spacer) ||
            !supported_calls::length_allowed(graph_calls.profile, call.spacer.size())) decision = "length";
        else if (!established.count(call.repeat)) decision = "low_family_diversity";
        else if (promoted.count(&call)) decision = "selected_recut";
        else if (retained.count(canonical)) decision = "graph_spacer";
        else if (graph_units.count(unit)) decision = "graph_boundary";
        else if (by_unit.at(unit).size() > 1) decision = "ambiguous_boundary";
        if (std::string(decision).compare(0, 8, "selected") == 0) {
            inventory.insert(canonical);
            entries << "read_first\t" << call.repeat << '\t' << call.spacer << '\t' << canonical
                    << '\t' << call.templates << '\t' << call.full_pairs << '\t' << call.repeat << "\t\n";
        }
        evidence << call.repeat << '\t' << call.spacer << '\t' << canonical << '\t'
                 << call.templates << '\t' << call.full_pairs << '\t' << decision << '\t'
                 << (index < within_read ? "one_read" : index < through_complementary ?
                     "complementary_junctions" : "variant_junctions") << '\n';
    }
    std::size_t id = 0;
    for (const auto& spacer : inventory) fasta << ">spacer_" << ++id << '\n' << spacer << '\n';
    for (auto* stream : {&entries, &fasta, &models, &evidence, &boundaries, &allele_table}) {
        stream->close();
        if (!*stream) throw std::runtime_error("Failed to write read-first result tables");
    }
    if (settings.read_first_family_recut) {
        family_log.close();
        if (!family_log) throw std::runtime_error("Failed to write read-first family recut evidence");
    }
    counts.try_emplace("read_first_graph_replaced_rows", 0);
    counts["read_first_added_spacer_types"] = std::count_if(inventory.begin(), inventory.end(),
        [&](const auto& type) { return !graph_calls.spacers.count(type); });
    counts["read_first_removed_spacer_types"] = std::count_if(graph_calls.spacers.begin(), graph_calls.spacers.end(),
        [&](const auto& type) { return !inventory.count(type); });
    counts["read_first_combined_spacer_types"] = inventory.size();
    std::cout << "  Read-first discovery: " << families.size() << " repeat families, "
              << counts["read_first_added_spacer_types"] << " additional spacer types, "
              << counts["read_first_graph_replaced_rows"] << " graph boundaries revised\n";
}

}} // namespace Zuck::read_first
#endif
