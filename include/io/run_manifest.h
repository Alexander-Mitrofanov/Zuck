#ifndef ZUCK_RUN_MANIFEST_H
#define ZUCK_RUN_MANIFEST_H

#include <map>
#include "core/build_identity.h"
#include "core/settings.h"
#include "io/result_bundle.h"

namespace result_bundle {

inline void write_manifest(const fs::path& directory, const std::string& input,
                           const std::string& input_checksum, const Settings& settings,
                           const std::map<std::string, std::size_t>& counts,
                           const std::map<std::string, double>& thresholds,
                           double seconds) {
    const fs::path path = directory / "run_manifest.json";
    std::ofstream out(path);
    if (!out) throw std::runtime_error("Cannot open run manifest: " + path.string());
    out << std::setprecision(17) << std::boolalpha;
    const std::string mode = settings.branch_aware_boundaries ? "branch-aware" :
                            settings.two_sided_boundaries ? "two-sided" : "legacy";
    out << "{\n  \"schema_version\": 1,\n  \"tool\": \"Zuck\",\n"
        << "  \"tool_version\": " << json_string(ZUCK_VERSION) << ",\n"
        << "  \"source_sha256\": " << json_string(ZUCK_SOURCE_SHA256) << ",\n"
        << "  \"compiler\": " << json_string(__VERSION__) << ",\n"
        << "  \"checksum_algorithm\": \"xxh64-seed0\",\n"
        << "  \"input_mode\": " << json_string(!settings.cycles_input.empty() ? "cycles" :
            !settings.graph_input.empty() ? "graph" : "reads") << ",\n"
        << "  \"cycle_input\": {\"path\": " << json_string(fs::absolute(input).string())
        << ", \"bytes\": " << fs::file_size(input)
        << ", \"checksum\": " << json_string(input_checksum) << "},\n"
        << "  \"upstream_inputs\": {\"reads\": " << json_string(settings.input_files)
        << ", \"graph\": " << json_string(settings.graph_input)
        << ", \"checksummed\": false},\n"
        << "  \"settings\": {\"boundary_mode\": " << json_string(mode)
        << ", \"short_spacer_filter\": " << settings.short_spacer_filter
        << ", \"short_spacer_filter_applied\": "
        << (settings.two_sided_boundaries && settings.short_spacer_filter)
        << ", \"terminal_variant_rescue\": " << settings.terminal_variant_rescue
        << ", \"terminal_variant_rescue_applied\": "
        << (settings.two_sided_boundaries && settings.terminal_variant_rescue)
        << ", \"read_linked_groups\": " << settings.read_linked_groups
        << ", \"array_database\": " << json_string(settings.array_database)
        << ", \"database_read_validation\": " << settings.database_read_validation
        << ", \"database_novel_spacers\": " << settings.database_novel_spacers
        << ", \"database_groups\": " << settings.database_groups
        << ", \"database_include_possible\": " << settings.database_include_possible
        << ", \"confidence_scores\": " << settings.confidence_scores
        << ", \"call_profile\": " << json_string(settings.call_profile)
        << ", \"family_repeat_context\": " << settings.family_repeat_context
        << (settings.read_first_family_recut ? ", \"read_first_family_recut\": true" : "")
        << (settings.read_first_family_admission ? ", \"read_first_family_admission\": true" : "")
        << (settings.read_first_shallow_rescue ? ", \"read_first_shallow_rescue\": true" : "")
        << (settings.read_first_sparse_junctions ? ", \"read_first_sparse_junctions\": true" : "")
        << ", \"read_recruitment\": " << settings.read_recruitment
        << ", \"local_boundary_graph\": " << settings.local_boundary_graph
        << ", \"motif_cliff_boundary\": " << settings.motif_cliff_boundary
        << ", \"singleton_precision_prune\": " << settings.singleton_precision_prune
        << ", \"circular_unit_dedup\": " << settings.circular_unit_dedup
        << ", \"directed_primary_groups\": " << settings.directed_primary_groups
        << ", \"cross_variant_scaffolds\": " << settings.cross_variant_scaffolds
        << ", \"repeat_mixture_refinement\": " << settings.repeat_mixture_refinement
        << ", \"candidate_evidence\": " << settings.candidate_evidence
        << ", \"boundary_radius\": " << settings.boundary_radius
        << ", \"boundary_policy\": " << json_string(settings.boundary_policy)
        << ", \"read_linked_groups_applied\": "
        << (counts.count("read_linkage_applied") && counts.at("read_linkage_applied") != 0)
        << ", \"coalesce_strand_equivalent_groups\": " << settings.two_sided_boundaries
        << ", \"experimental_alternative_views\": " << settings.experimental_alternative_views
        << ", \"threads\": " << settings.threads << ", \"ram_gib\": " << settings.ram
        << ", \"retain_graph\": " << settings.retain_graph
        << ", \"retain_cycles\": " << settings.retain_cycles << "},\n"
        << "  \"cycle_search_settings\": {\"applied_in_this_run\": " << settings.cycles_input.empty()
        << ", \"minimum_length\": " << settings.cycle_finder_settings.cycle_min_length
        << ", \"maximum_length\": " << settings.cycle_finder_settings.cycle_max_length
        << ", \"multiplicity_threshold\": " << settings.cycle_finder_settings.threshold_multiplicity
        << ", \"low_abundance\": " << settings.cycle_finder_settings.low_abundance
        << ", \"singleton_edges_removed\": " << !settings.cycle_finder_settings.low_abundance
        << ", \"traversal_partitions\": " << settings.cycle_finder_settings.traversal_partitions
        << ", \"partition_policy\": \"sorted_bucket_contiguous_balanced_v1\"},\n";
    if (settings.read_first)
        out << "  \"read_first\": {\"enabled\": true, \"inventory\": \"read_first_spacers.fasta\", "
               "\"evidence\": \"read_first_entries.tsv\"},\n";
    if (settings.read_first_family_recut)
        out << "  \"read_first_family_recut_evidence\": {\"ledger\": \"read_first_family_recut.tsv\", "
               "\"direct_model\": \"literal_whole_RSR_one_mate_either_strand\", \"minimum_base_quality\": 20, "
               "\"direct_minimum_distinct_templates\": 2, \"template_identity\": \"paired_sequence_fingerprint\", "
               "\"truncated_counts_are_lower_bounds\": true, \"truncated_geometry_abstains\": true},\n";
    if (settings.adaptive_coverage) {
        const auto& support = settings.coverage_support;
        const auto quartile = [&](unsigned value) {
            return support.repeated_edges ? std::to_string(value) : std::string("null");
        };
        out << "  \"coverage_support\": {\"policy\": \"repeated_edge_quartiles_v1\""
            << ", \"measured\": " << support.measured << ", \"low_support\": " << support.low_support
            << ", \"real_edges\": " << support.real_edges << ", \"repeated_edges\": " << support.repeated_edges
            << ", \"singleton_edges\": " << support.real_edges - support.repeated_edges
            << ", \"q25\": " << quartile(support.q25) << ", \"q50\": " << quartile(support.q50)
            << ", \"q75\": " << quartile(support.q75) << ", \"last_bin\": 255"
            << ", \"threshold_explicit\": " << settings.cycle_finder_settings.threshold_explicit
            << ", \"low_abundance_explicit\": " << settings.cycle_finder_settings.low_abundance_explicit
            << "},\n";
    }
    const auto& lengths = settings.dna_sequence_settings;
    out << "  \"output_length_constraints\": {\"applied\": " << lengths.active()
        << ", \"zero_means\": \"no_additional_bound\", \"scope\": \"rows_and_group_consensus_before_sidecars\""
        << ", \"spacer_minimum\": " << lengths.spacer_min_length
        << ", \"spacer_maximum\": " << lengths.spacer_max_length
        << ", \"repeat_minimum\": " << lengths.repeat_min_length
        << ", \"repeat_maximum\": " << lengths.repeat_max_length << "},\n";
    out << "  \"candidate_evidence_interpretation\": {\"enabled\": " << settings.candidate_evidence
        << ", \"table\": \"candidate_alleles.tsv\", \"unit\": \"exact_oriented_pre_collapse_pair\""
        << ", \"adds_primary_calls\": false, \"selection_uses_read_scores\": false"
        << ", \"traversal_counts_are_read_depth\": false, \"calibrated_probability\": false"
        << ", \"read_support_count_cap\": 32, \"database_only_candidates_included\": false},\n";
    out << "  \"local_boundary_graph_interpretation\": {\"enabled\": " << settings.local_boundary_graph
        << ", \"experimental\": true, \"table\": \"local_boundary_graph.tsv\""
        << ", \"scope\": \"family_aggregate_existing_primary_rows_only\", \"adds_primary_calls\": false"
        << ", \"extension_range\": [5,16], \"minimum_distinct_supported_family_rows\": 3"
        << ", \"owner_orientation\": \"candidate_repeat_contains_original_repeat_same_strand\""
        << ", \"boundary_cliff\": \"both_one_base_children_at_most_half_parent_row_support\""
        << ", \"application\": \"algebraic_content_preserving_recut_of_all_eligible_family_rows\""
        << ", \"cross_haplotype_transition_policy\": \"zero_supported_cross_transitions\""
        << ", \"minimum_sequence_units_per_family_candidate\": 1"
        << ", \"minimum_start_offsets_per_family_candidate\": 1"
        << ", \"reference_sequences_used\": false, \"calibrated_probability\": false},\n";
    out << "  \"motif_cliff_boundary_interpretation\": {\"enabled\": "
        << settings.motif_cliff_boundary
        << ", \"experimental\": true, \"table\": \"motif_cliff_boundary.tsv\""
        << ", \"scope\": \"exact_oriented_repeat_haplotypes_over_original_rows_after_exact_graph_recut\""
        << ", \"skip_exact_graph_owners\": true, \"adds_primary_calls\": false"
        << ", \"extension_range\": [5,16], \"minimum_distinct_original_spacers\": 3"
        << ", \"maximum_terminal_variants\": 2, \"minimum_top_variant_fraction\": 0.8"
        << ", \"maximum_child_support_ratio\": 0.7, \"both_sides_policy\": \"abstain\""
        << ", \"application\": \"circular_content_preserving_replacement_only\""
        << ", \"reference_sequences_used\": false, \"calibrated_probability\": false},\n";
    out << "  \"supported_call_interpretation\": {\"enabled\": " << (settings.call_profile != "off")
        << ", \"profile\": " << json_string(settings.call_profile)
        << ", \"decision_table\": \"call_decisions.tsv\", \"selected_table\": \"supported_entries.tsv\""
        << ", \"sequence_inventory\": \"supported_spacers.fasta\", \"family_table\": \"supported_arrays.tsv\""
        << ", \"primary_rows_preserved\": true, \"unit\": \"oriented_primary_or_supplemental_row\""
        << ", \"selection\": " << json_string(settings.read_recruitment
            ? "literal_context_or_verified_graph_seeded_read_recruitment"
            : "observable_context_with_at_least_two_sequence_units_and_start_offsets")
        << ", \"spacer_minimum\": " << (settings.call_profile == "standard" ? 26 : 0)
        << ", \"spacer_maximum\": " << (settings.call_profile == "standard" ? 50 : 0)
        << ", \"zero_length_bound_means\": \"no_profile_length_bound\""
        << ", \"read_span\": " << json_string(settings.family_repeat_context
            ? "whole_literal_repeat_spacer_observed_same_family_repeat_prefix16_or_whole_reverse_complement"
            : "whole_literal_repeat_spacer_next_repeat_prefix16_or_whole_reverse_complement")
        << ", \"next_repeat_allele_policy\": " << json_string(settings.family_repeat_context
            ? "same_output_family_same_length_at_most_four_substitutions_after_orientation"
            : "same_as_upstream_repeat")
        << ", \"minimum_quality\": 10, \"maximum_expected_errors\": 1"
        << ", \"independent_molecules\": false, \"calibrated_probability\": false"
        << ", \"original_memberships_retained\": true, \"singletons_retained_and_labeled\": true"
        << ", \"singleton_precision_prune\": {\"enabled\": "
        << settings.singleton_precision_prune
        << ",\"scope\":\"standard_selected_view_only\",\"inclusive_rows_retained\":true"
        << ",\"ledger\":\"singleton_precision_pruned.tsv\",\"canonical_type_rule\":\"all_selected_occurrences_are_literal_context_in_singleton_selected_families_with_sequence_units_at_most_two_and_start_offsets_at_most_two\"}"
        << ", \"circular_unit_dedup\": {\"enabled\": " << settings.circular_unit_dedup
        << ",\"scope\":\"standard_selected_view_only\",\"inclusive_rows_retained\":true"
        << ",\"ledger\":\"circular_unit_dedup.tsv\""
        << ",\"identity\":\"lexicographically_minimum_rotation_over_repeat_plus_spacer_and_reverse_complement\""
        << ",\"rank\":\"literal_owner_count_then_supported_owner_count_then_sum_min_max_balanced_owner_support_then_selected_row_count_then_lexicographic_spacer\""
        << ",\"cross_unit_type_conflict\":\"protect_globally\",\"reference_sequences_used\":false}"
        << ", \"canonical_inventory\": \"any_selected_row_does_not_validate_other_pairings\""
        << ", \"unknown_context\": [\"no_reads\",\"no_quality\",\"span_exceeds_max_quality_read\"]"
        << ", \"unknown_is_biological_rejection\": false"
        << ", \"locus_assignment\": \"unresolved\", \"spacer_order\": \"unspecified\"},\n";
    out << "  \"repeat_haplotype_group_interpretation\": {\"enabled\": "
        << (settings.call_profile != "off")
        << ", \"table\": \"repeat_haplotype_groups.tsv\""
        << ", \"scope\": \"selected_rows_within_original_family\""
        << ", \"partition\": \"exact_reverse_complement_canonical_observed_repeat_allele\""
        << ", \"minimum_distinct_spacers\": 2, \"content_calls_changed\": false"
        << ", \"original_families_preserved\": true, \"overlapping_groups\": false"
        << ", \"reference_sequences_used\": false},\n";
    out << "  \"directed_primary_group_interpretation\": {\"enabled\": "
        << settings.directed_primary_groups
        << ", \"membership_table\": \"directed_haplotype_paths.tsv\""
        << ", \"edge_table\": \"directed_haplotype_edges.tsv\""
        << ", \"scope\": \"final_selected_original_row_occurrences_after_boundary_and_pruning\""
        << ", \"repeat_partition\": \"exact_reverse_complement_canonical_observed_repeat_haplotype_within_parent\""
        << ", \"edge_direction\": \"ordered_left_spacer_to_right_spacer_across_exact_repeat\""
        << ", \"minimum_sequence_units\": 2, \"minimum_start_offsets\": 2"
        << ", \"edge_selection\": \"unique_reciprocal_best_by_min_support_then_total_support\""
        << ", \"reverse_direction_policy\": \"reject_any_observed_reverse_support\""
        << ", \"cycle_policy\": \"reject_all_cycle_edges\""
        << ", \"two_vertex_paths_allowed\": true, \"palindromic_repeat_policy\": \"residual_singletons\""
        << ", \"complete_selected_occurrence_partition\": true, \"row_disjoint\": true"
        << ", \"residual_policy\": \"one_singleton_group_per_unassigned_owner_haplotype_canonical_type_with_complete_row_fiber\""
        << ", \"recruitment_added_rows_may_create_edges\": false"
        << ", \"content_calls_changed\": false, \"reference_sequences_used\": false},\n";
    out << "  \"cross_variant_scaffold_interpretation\": {\"enabled\": "
        << settings.cross_variant_scaffolds
        << ", \"membership_table\": \"cross_variant_scaffolds.tsv\""
        << ", \"base_membership_table\": \"directed_haplotype_paths.tsv\""
        << ", \"read_link_membership_table\": \"read_linked_scaffolds.tsv\""
        << ", \"read_link_edge_table\": \"read_linked_scaffold_edges.tsv\""
        << ", \"read_link_component_table\": \"read_linked_scaffold_components.tsv\""
        << ", \"read_link_component_edge_table\": \"read_linked_scaffold_component_edges.tsv\""
        << ", \"scope\": \"final_selected_original_row_occurrences_after_boundary_and_pruning\""
        << ", \"merge_policy\": \"union_directed_components_sharing_a_native_path_or_deterministic_maximum_evidence_degree_two_forest_group_then_fallback_small_parent_partitions_to_content_family\""
        << ", \"cross_repeat_haplotype_scaffolds_allowed\": true"
        << ", \"minimum_sequence_units\": 2, \"minimum_start_offsets\": 2"
        << ", \"resolved_component_edge_order\": \"descending_min_support_then_total_support_then_max_support_then_sequence_units_then_start_offsets_then_lexical_endpoints\""
        << ", \"resolved_component_edge_acceptance\": \"accept_if_both_endpoint_degrees_below_two_and_edge_does_not_create_cycle\""
        << ", \"resolved_component_edge_ledger\": \"all_edges_recorded_with_selected_flag_and_selection_status\""
        << ", \"read_link_group_validation\": \"membership_nodes_equal_selected_edge_endpoints_complete_row_fibers_connected_degree_at_most_two_edges_equal_nodes_minus_one\""
        << ", \"parent_partition_gate\": \"preserve_read_partition_only_when_final_pre_gate_group_count_is_at_least_four_else_union_complete_content_family\""
        << ", \"parent_partition_gate_minimum_groups\": 4"
        << ", \"complete_selected_occurrence_partition\": true, \"row_disjoint\": true"
        << ", \"base_groups_split\": false, \"unlinked_base_groups_preserved\": true"
        << ", \"recruitment_added_rows_may_join_groups\": false"
        << ", \"content_calls_changed\": false, \"reference_sequences_used\": false},\n";
    out << "  \"read_recruitment_interpretation\": {\"enabled\": " << settings.read_recruitment
        << ", \"schema_version\": 4, \"evidence_table\": \"read_recruitment.tsv\""
        << ", \"candidate_table\": \"read_recruitment_candidates.tsv\""
        << ", \"seed_origin\": \"final_de_novo_graph_families\""
        << ", \"candidate_origin\": \"observed_complete_spacers_between_verified_repeat_copies\""
        << ", \"evidence_kinds\": [\"full_pair\",\"full_partial\",\"graph_junction\"]"
        << ", \"minimum_sequence_units\": 2, \"minimum_start_offsets\": 2"
        << ", \"minimum_family_spacers\": 2, \"support_count_cap\": 32"
        << ", \"seed_length\": 11, \"minimum_partial_repeat\": 16"
        << ", \"maximum_repeat_substitutions\": 1"
        << ", \"minimum_quality\": 10, \"maximum_expected_errors\": 1"
        << ", \"spacer_interior\": \"literal_observed_alleles_never_endpoint_collapsed\""
        << ", \"read_unit\": \"distinct_full_read_or_pair_canonical_under_reverse_complement_and_mate_swap\""
        << ", \"unit_identity\": \"two_independently_seeded_64_bit_hashes\""
        << ", \"offset_coordinate\": \"oriented_spacer_start_minus_16\""
        << ", \"pooling\": \"set_union_per_exact_oriented_repeat_spacer_pair\""
        << ", \"raw_triplets\": \"provenance_only\", \"aggregate_reverse_pair\": \"never_implied\""
        << ", \"maximum_candidates\": 65536, \"resource_overflow\": \"fail_closed_preserve_baseline\""
        << ", \"baseline_rows_preserved\": true, \"existing_core_recut_policy\": \"evidence_only\""
        << ", \"novel_ambiguous_family_policy\": \"evidence_only\""
        << ", \"legacy_context_measurements_preserved\": true"
        << ", \"tandem_ambiguity_model\": {\"version\": 1"
        << ", \"scope\": \"physical_raw_triplet_shared_by_both_oriented_pairs\""
        << ", \"geometry\": \"all_verified_observed_flank_bases\""
        << ", \"training\": \"left_full_then_right_full_reverse_complement_else_graph_proposed_left\""
        << ", \"model_order\": [\"short_period\",\"repeat_copy\"]"
        << ", \"motif_learning\": \"training_repeat_only_spacer_held_out\""
        << ", \"period_maximum\": \"floor_min_full_training_repeat_and_spacer_lengths_divided_by_two\""
        << ", \"similarity_threshold\": 0.82, \"threshold_comparison\": \"strictly_greater\""
        << ", \"spacer_distance\": \"levenshtein_to_cyclic_template_with_phase_continuation\""
        << ", \"flank_distance\": \"hamming_on_observed_bases\""
        << ", \"spacer_template_maximum\": \"ceil_50_times_spacer_length_divided_by_41\""
        << ", \"fit_coordinates\": \"chosen_training_orientation\""
        << ", \"repeat_diagnostics\": \"canonical_raw_left_repeat_all_3mer_windows\""
        << ", \"repeat_complexity_veto\": false"
        << ", \"non_tandem_support\": \"set_union_of_units_and_oriented_offsets_from_unexplained_raw_sources\""
        << ", \"minimum_non_tandem_sequence_units\": 2, \"minimum_non_tandem_start_offsets\": 2"
        << ", \"ambiguous_alleles\": \"retained_inclusive_without_recruitment_promotion\""
        << ", \"original_literal_selection_preserved\": true, \"calibrated_biological_probability\": false}"
        << ", \"error_ambiguity_model\": {\"version\": 2, \"alpha\": 0.05"
        << ", \"candidate_scope\": \"pre_gate_accepted_exact_pairs\""
        << ", \"opportunities\": \"all_quality_qualified_suffix16_repeat_spacer_length_arbitrary16_windows_both_orientations\""
        << ", \"quality\": \"minimum_spacer_phred_across_all_opportunity_windows\""
        << ", \"null_bound\": \"min_1_choose_N_K_times_p_power_K\""
        << ", \"successes\": \"unioned_non_tandem_sequence_units\""
        << ", \"tandem_ambiguous_status\": \"not_evaluated_tandem_ambiguous\""
        << ", \"multiplicity\": \"three_times_sum_pre_gate_accepted_pair_spacer_lengths\""
        << ", \"neighbor\": \"same_exact_oriented_repeat_same_length_hamming_one\""
        << ", \"independent_calibrated_errors_assumed\": true, \"calibrated_biological_probability\": false"
        << ", \"ambiguous_alleles\": \"retained_inclusive_without_recruitment_promotion\""
        << ", \"original_literal_selection_preserved\": true, \"new_ambiguous_literal_bypass\": false}"
        << ", \"row_promotion_model\": {\"version\":1,\"scope\":\"pre_recruitment_exact_array_oriented_repeat_spacer_rows\",\"original_rows_table\":\"read_recruitment_original_rows.tsv\",\"context_table\":\"read_recruitment_context.tsv\",\"candidate_evidence\":\"schema3_tandem_and_error_eligibility_unchanged\",\"context_universe\":\"evidence_eligible_candidates_with_actual_original_row_proof_binding_before_row_policy\",\"owner_scope\":\"every_declared_candidate_owner\",\"repeat_neighbor\":\"equal_length_minimum_hamming_to_either_strand_at_most_one\",\"focal_spacer_distance\":\"minimum_levenshtein_to_either_strand_divided_by_maximum_length\",\"minimum_spacer_distance\":0.18,\"distance_comparison\":\"greater_than_or_equal\",\"integer_distance_test\":\"50_times_distance_at_least_9_times_maximum_length\",\"deterministic_context\":\"lowest_candidate_id_then_lowest_original_entry_index\",\"reverse_complement_ties\":\"forward\",\"recursive_context_pruning\":false,\"novel_pairs\":\"retained_inclusive_boundary_unresolved\",\"original_literal_selection_preserved\":true,\"new_novel_literal_bypass\":false,\"read_linkage\":\"inferred_before_recruitment\",\"supplemental_views\":\"recruitment_added_rows_excluded\",\"family_context\":\"recruitment_added_rows_cannot_supply_repeat_alleles\",\"status_precedence\":[\"novel_boundary_unresolved\",\"tandem_ambiguous\",\"error_ambiguous\",\"insufficient_spacer_diversity\",\"supported\"],\"calibrated_biological_probability\":false}"
        << ", \"independent_molecules\": false, \"calibrated_probability\": false},\n";
    out << "  \"boundary_candidate_interpretation\": {\"applied\": " << (counts.count("boundary_seed_pairs") != 0)
        << ", \"independent_junction_shifts\": true, \"reference_recall_used_for_inference\": false"
        << ", \"selection_is_additive\": true, \"same_source_array_anchors\": true"
        << ", \"known_focal_spacer_and_core_excluded_from_corroboration\": true"
        << ", \"read_support_count_cap\": 32, \"calibrated_probability\": false},\n";
    out << "  \"repeat_mixture_interpretation\": {\"enabled\": " << settings.repeat_mixture_refinement
        << ", \"experimental\": true, \"maximum_partition_levels\": 2"
        << ", \"additional_signature_bases\": 5, \"minimum_between_branch_edits\": 3"
        << ", \"minimum_distinct_cores_per_branch\": 3, \"maximum_branches\": 4"
        << ", \"minimum_core_coverage_fraction\": 0.8, \"minimum_repeat_extension\": 5"
        << ", \"reference_sequences_used\": false"
        << ", \"conserved_spacer_termini_can_be_ambiguous\": true},\n";
    out << "  \"grouping_read_inputs\": [";
    const auto grouping_inputs = settings.effective_grouping_reads();
    for (std::size_t i = 0; i < grouping_inputs.size(); ++i) {
        if (i) out << ',';
        out << json_string(fs::absolute(grouping_inputs[i]).string());
    }
    out << "],\n";
    out << "  \"array_database_input\": ";
    if (settings.array_database.empty()) out << "null";
    else out << "{\"path\": " << json_string(fs::absolute(settings.array_database).string())
             << ", \"bytes\": " << fs::file_size(settings.array_database)
             << ", \"checksum\": " << json_string(file_checksum(settings.array_database)) << "}";
    out << ",\n  \"database_interpretation\": {"
        << "\"mode\": \"additive_exact_cycle_repeat_spacer_support\","
        << "\"groups\": \"exact_repeat_families_not_verified_loci\","
        << "\"de_novo_groups_and_linkage_preserved\": true,"
        << "\"source_confidence_is_calibrated_probability\": false,"
        << "\"duplicates_accumulate_votes\": false,"
        << "\"requires_same_array_spacer_matches\": true,"
        << "\"known_spacer_support\": \"exact_same_source_array\","
        << "\"novel_spacer_support\": \"trusted_repeat_with_other_known_same_array_spacers\","
        << "\"read_span\": \"whole_literal_repeat_spacer_next_repeat_prefix16_or_whole_reverse_complement\","
        << "\"read_unit\": \"distinct_full_read_or_ordered_pair_sequences\","
        << "\"read_units_are_independent_molecules\": false,"
        << "\"tied_boundary_interpretations\": \"abstain\","
        << "\"score\": \"category_weight * confidence * (2*min(distinct_spacer_matches,4)/4 + min(distinct_repeat_matches,2)/2)\"},\n";
    out << "  \"read_linkage_interpretation\": {"
        << "\"unit\": " << json_string(grouping_inputs.size() == 2 ? "distinct_exact_sequence_pair" : "distinct_exact_read_sequence")
        << ", \"duplicate_definition\": \"identical full sequences in the same input positions\""
        << ", \"independent_molecules\": false"
        << ", \"view_kind\": \"supported_unbranched_content_subset\""
        << ", \"primary_rows_preserved\": true"
        << ", \"locus_assignment\": \"unresolved\", \"biological_spacer_order\": \"unspecified\"},\n";
    // Search settings cannot describe how an externally supplied cycle file was made.
    out << "  \"postprocessing_thresholds\": {";
    bool first = true;
    for (const auto& [name, value] : thresholds) {
        if (!first) out << ',';
        out << "\n    " << json_string(name) << ": " << value;
        first = false;
    }
    out << "\n  },\n  \"counts\": {";
    first = true;
    for (const auto& [name, count] : counts) {
        if (!first) out << ',';
        out << "\n    " << json_string(name) << ": " << count;
        first = false;
    }
    out << "\n  },\n  \"postprocessing_seconds\": " << seconds << ",\n"
        << "  \"interpretation\": {\"unit\": \"sequence_group\", "
        << "\"locus_assignment\": \"unresolved\", \"spacer_order\": \"unspecified\", "
        << "\"row_orientation\": \"preserved within emitted pairs; may differ within a group\", "
        << "\"terminal_repeat\": \"unobserved\", \"support_counts_are_read_depth\": false},\n"
        << "  \"scoring\": ";
    if (!settings.confidence_scores) out << "null";
    else out << "{\"schema_version\":1,\"table\":\"element_confidence.tsv\","
             << "\"kind\":\"heuristic_evidence_indices\",\"calibrated_probability\":false,"
             << "\"primary_rows_preserved\":true,\"independent_molecules\":false,"
             << "\"read_quality_minimum\":10,\"maximum_expected_errors\":1,"
             << "\"support_count_cap\":32,\"unknown_read_support\":\"NA\","
             << "\"read_index\":\"min(sequence_units,4)/4 * min(start_offsets,2)/2\","
             << "\"database_index\":\"category_weight * confidence * min(min(other_spacers,other_cores),3)/3\","
             << "\"boundary_index\":\"database_index / catalogue_supported_decompositions\","
             << "\"database_exclusion\":\"focal_spacer_and_circular_core\","
             << "\"boundary_endpoints_separately_scored\":false,"
             << "\"boundary_scope\":\"emitted_competing_circular_core_splits\"}";
    out << ",\n  \"files\": [";
    std::vector<fs::path> paths;
    for (const auto& entry : fs::directory_iterator(directory))
        if (entry.path() != path) paths.push_back(entry.path());
    std::sort(paths.begin(), paths.end());
    first = true;
    for (const auto& file : paths) {
        if (!first) out << ',';
        out << "\n    {\"path\": " << json_string(file.filename().string())
            << ", \"bytes\": " << fs::file_size(file)
            << ", \"checksum\": " << json_string(file_checksum(file)) << '}';
        first = false;
    }
    out << "\n  ]\n}\n";
    close_checked(out, path);
}
} // namespace result_bundle
#endif
