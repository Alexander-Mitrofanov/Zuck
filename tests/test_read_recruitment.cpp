#include "crispr_array/read_recruitment.h"
#include <filesystem>
#include <fstream>
#include <random>
#include <sstream>
#include <stdexcept>

namespace rr = Zuck::read_recruitment;
namespace rl = Zuck::read_linked;
struct Array { std::string consensus; std::vector<rr::Pair> entries; };
static void require(bool condition, const char* reason) {
    if (!condition) throw std::runtime_error(reason);
}
static std::set<std::string> selected(const rr::Result& result) {
    std::set<std::string> spacers;
    for (const auto& candidate : result.candidates)
        if (candidate.accepted) spacers.insert(rl::canonical(candidate.spacer));
    return spacers;
}
static std::set<std::string> promoted(const rr::Result& result) {
    std::set<std::string> spacers;
    for (const auto& candidate : result.candidates)
        if (candidate.selection_eligible) spacers.insert(rl::canonical(candidate.spacer));
    return spacers;
}
static std::string inventory(const rr::Result& result) {
    std::ostringstream out;
    for (const auto& candidate : result.candidates) {
        out << candidate.evidence_kind << ' ' << candidate.left_repeat << ' ' << candidate.spacer << ' '
            << candidate.right_repeat << ' ' << candidate.support.sequence_units << ' '
            << candidate.support.start_offsets << ' ' << candidate.accepted << ' ' << candidate.reason << '\n';
        out << candidate.non_tandem_support.sequence_units << ' ' << candidate.non_tandem_support.start_offsets << ' '
            << candidate.tandem_ambiguous << ' ' << candidate.tandem_status << ' ' << candidate.tandem.model << ' '
            << candidate.tandem.training_kind << ' ' << candidate.tandem.motif << ' ' << candidate.tandem.period << ' '
            << candidate.tandem.spacer_template_length << ' ' << candidate.tandem.whole_similarity << '\n';
        for (auto owner : candidate.family_owners) out << owner << ',';
    }
    return out.str();
}
int main() {
    const auto root = std::filesystem::temp_directory_path() /
        ("zuck-read-recruitment-" + std::to_string(std::random_device{}()));
    std::filesystem::create_directory(root);
    try {
        std::mt19937 random(916113);
        auto dna = [&](std::size_t length) {
            std::string sequence;
            while (length--) sequence += "ACGT"[random() % 4];
            return sequence;
        };
        const auto repeat = dna(30);
        const std::vector<std::string> spacers{dna(32), dna(33), dna(31), dna(32), dna(34)};
        const Array family{repeat, {{repeat, spacers[0]}, {repeat, spacers[1]}}};
        const std::vector<Array> arrays{family};
        const auto original = family.entries;
        auto fastq = [&](const std::string& name, const std::vector<std::string>& reads, int q = 40) {
            const auto path = root / name;
            std::ofstream out(path);
            for (std::size_t i = 0; i < reads.size(); ++i)
                out << '@' << i << '\n' << reads[i] << "\n+\n" << std::string(reads[i].size(), char(q + 33)) << '\n';
            return path.string();
        };
        auto full = [&](const std::string& spacer, const std::string& next) { return repeat + spacer + next; };
        auto witness = [&](std::vector<std::string>& reads, const std::string& sequence, bool terminal = false) {
            reads.push_back("A" + sequence + (terminal ? "" : "C"));
            reads.push_back("GGG" + sequence + (terminal ? "" : "TTT"));
        };
        auto run = [&](const std::vector<std::string>& reads, const std::vector<Array>& seeds = std::vector<Array>{},
                       const rr::Config& config = rr::Config{}) {
            return rr::infer(seeds.empty() ? arrays : seeds, {fastq("reads.fq", reads)}, config);
        };
        std::string variant = repeat;
        variant[12] = variant[12] == 'A' ? 'C' : 'A';
        std::vector<std::string> reads;
        witness(reads, full(spacers[2], variant));
        witness(reads, full(spacers[3], repeat));
        const auto result = run(reads);
        require(selected(result) == std::set<std::string>{rl::canonical(spacers[2]), rl::canonical(spacers[3])},
                "Full repeat pairs must recover exact graph-missing spacers across one repeat SNP");
        require(promoted(result) == selected(result), "Distinct full-repeat contexts failed promotion");
        require(family.entries == original, "Read recruitment mutated graph rows");
        bool saw_variant = false;
        for (const auto& candidate : result.candidates) {
            require(candidate.accepted == (candidate.pair_only && candidate.support.accepted() && !candidate.boundary_conflict &&
                !candidate.resource_limited && candidate.distinct_supported_family_spacers >= 2), "Acceptance contract mismatch");
            require(candidate.selection_eligible == (candidate.accepted && !candidate.error_ambiguous && !candidate.tandem_ambiguous),
                    "Selection contract does not distinguish observation from ambiguity");
            if (rl::canonical(candidate.spacer) != rl::canonical(spacers[2])) continue;
            for (const auto& pair : rr::oriented_pairs(candidate))
                saw_variant |= pair == rr::Pair{rl::reverse_complement(variant), rl::reverse_complement(spacers[2])};
        }
        require(saw_variant, "Reverse observed adjacency failed to swap flanking repeats");
        auto reversed_reads = reads;
        for (auto& read : reversed_reads) read = rl::reverse_complement(read);
        require(inventory(run(reversed_reads)) == inventory(result), "Read reverse complement changed the evidence graph");
        std::reverse(reads.begin(), reads.end());
        require(inventory(run(reads)) == inventory(result), "Read arrival order changed candidates or support");
        auto duplicate_reads = reads;
        duplicate_reads.insert(duplicate_reads.end(), reversed_reads.begin(), reversed_reads.end());
        require(inventory(run(duplicate_reads)) == inventory(result), "Read/RC duplicates inflated support");

        reads.clear();
        witness(reads, full(spacers[2], repeat.substr(0, 16)), true);
        witness(reads, full(spacers[3], repeat.substr(0, 19)), true);
        const auto partial = run(reads);
        require(selected(partial).size() == 2, "Terminal next-repeat prefixes must recover complete novel spacers");
        require(promoted(partial) == selected(partial), "Distinct partial-repeat contexts failed promotion");
        for (const auto& candidate : partial.candidates) if (candidate.accepted) {
            require(candidate.evidence_kind == "full_partial" && !candidate.right_complete &&
                    candidate.right_repeat.size() >= 16 && rr::oriented_pairs(candidate).size() == 1,
                    "Partial repeat evidence was promoted to a full invented repeat");
        }
        reversed_reads = reads;
        for (auto& read : reversed_reads) read = rl::reverse_complement(read);
        require(inventory(run(reversed_reads)) == inventory(partial), "Clipped-repeat evidence depends on read strand");
        auto contradictory = repeat.substr(0, 20);
        contradictory[17] = contradictory[17] == 'A' ? 'C' : 'A';
        contradictory[18] = contradictory[18] == 'A' ? 'C' : 'A';
        reads.clear();
        witness(reads, full(spacers[2], contradictory), true);
        witness(reads, full(spacers[3], contradictory), true);
        require(selected(run(reads)).empty(), "Partial evidence ignored contradictory available repeat bases");
        reads.clear();
        witness(reads, repeat + spacers[2].substr(0, 20), true);
        witness(reads, repeat + spacers[3].substr(0, 20), true);
        require(selected(run(reads)).empty(), "Incomplete spacer interior was invented");

        reads.clear();
        for (int s : {0, 1}) {
            reads.push_back(repeat.substr(14) + spacers[s] + repeat.substr(0, 16));
            reads.push_back(repeat.substr(13) + spacers[s] + repeat.substr(0, 17));
        }
        const auto junction = run(reads);
        require(selected(junction).size() == 2, "Exact graph proposals need support from two clipped repeat junctions");
        require(promoted(junction) == selected(junction), "Distinct graph-junction contexts failed promotion");
        for (const auto& candidate : junction.candidates) if (candidate.accepted)
            require(candidate.graph_proposed && candidate.evidence_kind == "graph_junction" &&
                    candidate.support.start_offsets == 2, "Clipped junction offsets or evidence kind were lost");
        auto unproposed = arrays;
        unproposed[0].entries = {{repeat, spacers[2]}, {repeat, spacers[3]}};
        require(selected(run(reads, unproposed)).empty(), "Two partial repeats invented a graph-missing spacer");

        reads.clear();
        witness(reads, full(spacers[2], repeat));
        witness(reads, full(spacers[3], repeat));
        require(selected(rr::infer(arrays, {fastq("low.fq", reads, 9)})).empty(), "Low-quality evidence accepted");
        require(selected(rr::infer(arrays, {fastq("ee.fq", reads, 10)})).empty(), "Expected-error limit ignored");
        std::vector<std::string> same_offsets;
        for (int s : {2, 3}) {
            same_offsets.push_back("A" + full(spacers[s], repeat) + "C");
            same_offsets.push_back("G" + full(spacers[s], repeat) + "T");
        }
        require(selected(run(same_offsets)).empty(), "Different sequences at one offset manufactured support");
        std::vector<std::string> only_one;
        witness(only_one, full(spacers[2], repeat));
        require(selected(run(only_one)).empty(), "Homogeneous tandem sequence passed family diversity");
        std::vector<std::string> mate1{reads[0], rl::reverse_complement(reads[1]), reads[2], rl::reverse_complement(reads[3])};
        std::vector<std::string> mate2{reads[1], rl::reverse_complement(reads[0]), reads[3], rl::reverse_complement(reads[2])};
        const auto paired = rr::infer(arrays, {fastq("r1.fq", mate1), fastq("r2.fq", mate2)});
        require(selected(paired).empty(), "Swapped reverse-complement read pairs counted as independent units");
        for (const auto& candidate : paired.candidates) require(candidate.support.sequence_units == 1, "Pair unit deduplication failed");

        // Pool independent witnesses for the exact emitted pair across evidence
        // classes and next-repeat alleles. Each direction retains its own proof.
        std::vector<std::string> mixed;
        for (int s : {2, 3}) {
            mixed.push_back("A" + full(spacers[s], repeat) + "C");
            mixed.push_back("GGG" + full(spacers[s], repeat.substr(0, 16)));
        }
        const auto mixed_result = run(mixed);
        require(selected(mixed_result).size() == 2, "Full and partial witnesses for one exact pair were not united");
        require(promoted(mixed_result) == selected(mixed_result), "Mixed non-tandem evidence was not promoted");
        for (const auto& candidate : mixed_result.candidates) if (candidate.accepted)
            require(candidate.pair_only && candidate.evidence_kind == "full_partial" &&
                    candidate.source_candidate_indices.size() >= 2 && rr::oriented_pairs(candidate).size() == 1,
                    "Pooled evidence lost provenance or licensed an unverified reverse pair");
        std::vector<std::string> different_right;
        for (int s : {2, 3}) {
            different_right.push_back("A" + full(spacers[s], repeat) + "C");
            different_right.push_back("GGG" + full(spacers[s], variant) + "TTT");
        }
        const auto pooled_right = run(different_right);
        require(selected(pooled_right).size() == 2, "Observed next-repeat variants split valid upstream support");
        bool pooled_upstream = false;
        const rr::Pair rare_reverse{rl::reverse_complement(variant), rl::reverse_complement(spacers[2])};
        for (const auto& candidate : pooled_right.candidates) {
            const auto pairs = rr::oriented_pairs(candidate);
            if (candidate.pair_only && candidate.left_repeat == repeat && candidate.spacer == spacers[2]) {
                require(candidate.support.sequence_units == 2 && candidate.accepted && candidate.observed_right_repeats.size() == 2,
                        "Exact upstream pair failed to combine next-repeat alleles");
                pooled_upstream = true;
            }
            if (std::find(pairs.begin(), pairs.end(), rare_reverse) != pairs.end())
                require(!candidate.accepted && candidate.support.sequence_units == 1, "Rare reverse repeat borrowed unrelated upstream proof");
        }
        require(pooled_upstream, "Missing directed aggregate");
        std::vector<std::string> one_reverse_offset;
        for (int s : {2, 3}) {
            one_reverse_offset.push_back("A" + full(spacers[s], repeat) + "C");
            one_reverse_offset.push_back("GGG" + full(spacers[s], repeat) + "T");
        }
        const auto directed = run(one_reverse_offset);
        require(selected(directed).size() == 2, "Valid forward placement lost support");
        for (const auto& candidate : directed.candidates) {
            if (!candidate.accepted) continue;
            require(candidate.pair_only && candidate.left_repeat == repeat,
                    "Raw canonical-triplet offsets licensed an unsupported reverse placement");
        }
        std::vector<std::string> mixed_one_offset;
        for (int s : {2, 3}) {
            mixed_one_offset.push_back("A" + full(spacers[s], repeat) + "C");
            mixed_one_offset.push_back("A" + full(spacers[s], repeat.substr(0, 16)));
        }
        const auto shared_offset = run(mixed_one_offset);
        require(selected(shared_offset).empty(), "Different evidence coordinate origins manufactured start offsets");
        for (const auto& candidate : shared_offset.candidates)
            if (candidate.pair_only && candidate.left_repeat == repeat)
                require(candidate.support.sequence_units == 2 && candidate.support.start_offsets == 1,
                        "Pooled offsets are not normalized to the same spacer anchor");
        const auto pooled_mates = rr::infer(arrays, {
            fastq("pool1.fq", {mixed[0], mixed[2]}), fastq("pool2.fq", {mixed[1], mixed[3]})});
        require(selected(pooled_mates).empty(), "Two evidence classes from one paired unit manufactured independent support");
        for (const auto& candidate : pooled_mates.candidates)
            require(candidate.support.sequence_units == 1, "Pair aggregation added counts instead of unioning unit identities");

        auto ambiguous = arrays;
        ambiguous.push_back(family);
        const auto ambiguity = run(reads, ambiguous);
        require(!ambiguity.candidates.empty(), "Ambiguous family evidence disappeared");
        for (const auto& candidate : ambiguity.candidates)
            require(candidate.family_owners == std::set<std::size_t>{0, 1}, "An arbitrary seed family won ownership");

        // Same circular DNA, a competing cut. Both interpretations stay in the
        // ledger; read-span support cannot decide the biological boundary.
        const auto core = repeat + spacers[2];
        const auto shifted_repeat = core.substr(0, repeat.size() + 2);
        const auto shifted_spacer = core.substr(repeat.size() + 2);
        auto protected_arrays = arrays;
        protected_arrays[0].entries.push_back({shifted_repeat, shifted_spacer});
        const auto conflict = run(reads, protected_arrays);
        bool conflict_seen = false;
        for (const auto& candidate : conflict.candidates)
            if (rl::canonical(candidate.spacer) == rl::canonical(spacers[2])) {
                conflict_seen |= candidate.boundary_conflict;
                require(!candidate.accepted, "A supported read silently recut an existing graph core");
            }
        require(conflict_seen, "Alternative circular cuts were not recorded");
        const auto body = dna(27);
        const auto prefix = std::string("ATTGC");
        const auto complete_spacer = prefix + body;
        const std::vector<Array> unequal_seeds{{repeat, {
            {repeat, complete_spacer}, {repeat, spacers[0]}, {repeat + prefix, spacers[1]}}}};
        std::vector<std::string> unequal_reads;
        witness(unequal_reads, full(complete_spacer, repeat));
        witness(unequal_reads, full(spacers[0], repeat));
        const auto unequal = run(unequal_reads, unequal_seeds);
        bool tainted_pool_seen = false;
        for (const auto& candidate : unequal.candidates)
            if (candidate.pair_only && rl::canonical(candidate.spacer) == rl::canonical(body)) {
                require(candidate.boundary_conflict && !candidate.accepted,
                        "Pair aggregation laundered a competing cut visible only on the opposite repeat");
                tainted_pool_seen = true;
            }
        require(tainted_pool_seen, "Unequal-flank boundary fixture did not exercise pooled conflict inheritance");

        // Distinct exact interior alleles with identical terminal anchors survive.
        auto minority = spacers[2]; minority[16] = minority[16] == 'A' ? 'C' : 'A';
        std::vector<std::string> alleles;
        witness(alleles, full(spacers[2], repeat));
        witness(alleles, full(minority, repeat));
        require(selected(run(alleles)).size() == 2, "Shared endpoints collapsed a true internal spacer allele");

        // The error screen conditions on independent calibrated substitution
        // errors. It never erases an exact observed allele or changes support.
        auto quality_panel = [&](std::size_t major_count, std::size_t minor_count, int q) {
            std::vector<std::string> panel;
            for (std::size_t i = 1; i <= major_count; ++i)
                panel.push_back(std::string(i, 'A') + full(spacers[2], repeat) + std::string(i, 'C'));
            for (std::size_t i = 1; i <= minor_count; ++i)
                panel.push_back(std::string(i, 'G') + full(minority, repeat) + std::string(i, 'T'));
            return std::make_pair(panel, rr::infer(arrays, {fastq("quality.fq", panel, q)}));
        };
        const auto high_quality = quality_panel(8, 2, 40);
        require(promoted(high_quality.second).size() == 2, "High-quality independently supported minority allele was collapsed");
        const auto error_plausible = quality_panel(8, 2, 20);
        require(selected(error_plausible.second).size() == 2 &&
                promoted(error_plausible.second) == std::set<std::string>{rl::canonical(spacers[2])},
                "Plausible error must remain inclusive while only automatic promotion is withheld");
        for (const auto& candidate : error_plausible.second.candidates) if (candidate.accepted) {
            require(candidate.opportunity_records == 10 && candidate.opportunity_min_spacer_quality == 20 &&
                    candidate.opportunity_complete && candidate.error_hypotheses == 12 * spacers[2].size() &&
                    candidate.error_neighbor_count == 1, "Opportunity or comparison-universe accounting is incorrect");
            if (rl::canonical(candidate.spacer) == rl::canonical(minority)) {
                require(candidate.error_ambiguous && candidate.error_status == "error_ambiguous" &&
                        std::abs(candidate.error_log_bound - std::log(45.0 * 0.01 * 0.01)) < 1e-10,
                        "Independent-error union bound was not computed in log space");
            }
        }
        require(promoted(quality_panel(32, 8, 20).second).size() == 2,
                "Strong low-quality minority evidence must be distinguished from sparse coincident errors");
        auto extra_opportunities = high_quality.first;
        extra_opportunities.push_back(repeat + std::string(spacers[2].size(), 'N') + std::string(16, 'T'));
        const auto arbitrary = rr::infer(arrays, {fastq("arbitrary.fq", extra_opportunities, 40)});
        for (const auto& candidate : arbitrary.candidates) if (candidate.accepted && candidate.left_repeat == repeat)
            require(candidate.opportunity_records == 11, "Opportunity denominator excluded arbitrary spacer/next bases");
        auto duplicated_source = error_plausible.first;
        for (int copy = 0; copy < 8; ++copy)
            duplicated_source.insert(duplicated_source.end(), error_plausible.first.begin(), error_plausible.first.begin() + 8);
        const auto duplicate_opportunities = rr::infer(arrays, {fastq("duplicate-source.fq", duplicated_source, 20)});
        for (const auto& candidate : duplicate_opportunities.candidates)
            if (candidate.accepted && rl::canonical(candidate.spacer) == rl::canonical(minority))
                require(candidate.error_ambiguous && candidate.opportunity_records == 74 &&
                        candidate.support.sequence_units == 2, "Duplicate source opportunities strengthened an uncertain allele");
        rr::Config opportunity_limit;
        opportunity_limit.maximum_opportunity_matches_per_read = 1;
        const auto incomplete = rr::infer(arrays, {fastq("limited-opportunities.fq", high_quality.first)}, opportunity_limit);
        require(selected(incomplete).size() == 2 && promoted(incomplete).empty(),
                "Incomplete error opportunities erased observations or promoted incomplete claims");
        for (const auto& candidate : incomplete.candidates) if (candidate.accepted)
            require(candidate.error_status == "opportunity_unresolved" && !candidate.opportunity_complete,
                    "Incomplete statistical evidence was not explicit");
        require(incomplete.counts.at("read_recruitment_opportunity_read_records") == high_quality.first.size(),
                "Opportunity overflow skipped input validation");
        require(std::isfinite(rr::independent_error_log_bound(1000000, 32, 1e-9)),
                "Small-tail error bound lost numerical stability");

        // A similar allele belonging to a different exact repeat is not a
        // competing source, and sub-threshold-quality reads cannot create one.
        const auto other_repeat = dna(30);
        std::vector<Array> different_parent = arrays;
        different_parent.push_back({other_repeat, {{other_repeat, spacers[2]}, {other_repeat, spacers[0]}}});
        std::vector<std::string> separate;
        witness(separate, full(minority, repeat));
        witness(separate, full(spacers[3], repeat));
        witness(separate, other_repeat + spacers[2] + other_repeat);
        witness(separate, other_repeat + spacers[0] + other_repeat);
        const auto separate_families = run(separate, different_parent);
        for (const auto& candidate : separate_families.candidates)
            if (candidate.accepted && candidate.left_repeat == repeat && candidate.spacer == minority)
                require(candidate.error_neighbor_count == 0 && candidate.selection_eligible,
                        "Other-repeat allele was used as an error source");
        const auto mixed_quality_path = root / "mixed-quality.fq";
        {
            std::ofstream out(mixed_quality_path);
            std::vector<std::string> good;
            witness(good, full(minority, repeat)); witness(good, full(spacers[3], repeat));
            std::vector<std::string> bad;
            witness(bad, full(spacers[2], repeat));
            std::size_t id = 0;
            for (const auto& sequence : good)
                out << '@' << id++ << '\n' << sequence << "\n+\n" << std::string(sequence.size(), 'I') << '\n';
            for (const auto& sequence : bad)
                out << '@' << id++ << '\n' << sequence << "\n+\n" << std::string(sequence.size(), '*') << '\n';
        }
        const auto low_parent = rr::infer(arrays, {mixed_quality_path.string()});
        for (const auto& candidate : low_parent.candidates)
            if (candidate.accepted && candidate.left_repeat == repeat && candidate.spacer == minority)
                require(candidate.error_neighbor_count == 0 && candidate.selection_eligible,
                        "Low-quality competitor licensed an error veto");

        // V3 retains observations but needs independent evidence that lacks a
        // repeat-trained tandem explanation before automatic promotion.
        auto raw_candidate = [](const std::string& left, const std::string& spacer,
                                const std::string& right, std::size_t left_bases = 0,
                                std::size_t right_bases = 0) {
            rr::Candidate candidate;
            candidate.left_repeat = left; candidate.spacer = spacer; candidate.right_repeat = right;
            candidate.left_observed_bases = left_bases ? left_bases : left.size();
            candidate.right_observed_bases = right_bases ? right_bases : right.size();
            candidate.right_complete = right.size() >= left.size();
            candidate.evidence_kind = candidate.left_observed_bases < left.size() ? "graph_junction" :
                candidate.right_complete ? "full_pair" : "full_partial";
            return candidate;
        };
        const std::string motif9 = "ACGTTGCAA";
        const auto periodic_repeat = motif9 + motif9 + motif9;
        const auto periodic_spacer = periodic_repeat + motif9;
        const auto periodic_long = periodic_spacer + motif9;
        const auto periodic_full = rr::tandem_evidence(raw_candidate(periodic_repeat, periodic_spacer, periodic_repeat));
        const auto periodic_partial = rr::tandem_evidence(raw_candidate(periodic_repeat, periodic_spacer,
            periodic_repeat.substr(0, 16)));
        const auto periodic_graph = rr::tandem_evidence(raw_candidate(periodic_repeat, periodic_spacer,
            periodic_repeat, 16, 16));
        require(periodic_full.explained && periodic_partial.explained && periodic_graph.explained &&
                periodic_full.period == 9 && periodic_partial.period == 9 && periodic_graph.period == 9,
                "Clipping a period-nine path bypassed the tandem model universe");
        require(periodic_graph.training_kind == "graph_proposed_repeat" &&
                periodic_partial.training_kind == "observed_full_repeat",
                "Proposed unobserved repeat bases were described as read observations");
        // Extra verified flank bases must be retained even when the first 16
        // alone support a periodic explanation. The complete graph repeat is a
        // proposal; neither clipped repeat here is physically full.
        const std::string motif8 = "ACGTTCGA";
        auto extended_repeat = motif8 + motif8 + motif8 + motif8 + motif8;
        for (std::size_t i = 16; i < 20; ++i)
            extended_repeat[i] = extended_repeat[i] == 'A' ? 'C' : 'A';
        const auto extended_spacer = motif8 + motif8 + motif8 + motif8;
        require(rr::tandem_evidence(raw_candidate(extended_repeat, extended_spacer, extended_repeat, 16, 16)).explained &&
                !rr::tandem_evidence(raw_candidate(extended_repeat, extended_spacer, extended_repeat, 16, 20)).explained,
                "Verified bases beyond the minimum anchor were discarded by the model");
        const std::vector<Array> extended_seeds{{extended_repeat,
            {{extended_repeat, extended_spacer}, {extended_repeat, spacers[0]}}}};
        const auto actual_geometry = run({extended_repeat.substr(24) + extended_spacer + extended_repeat.substr(0, 20),
            extended_repeat.substr(23) + extended_spacer + extended_repeat.substr(0, 21)}, extended_seeds);
        bool saw_actual_geometry = false;
        for (const auto& candidate : actual_geometry.candidates)
            if (!candidate.pair_only && candidate.evidence_kind == "graph_junction" && candidate.left_repeat == extended_repeat &&
                candidate.left_observed_bases == 16 && candidate.right_observed_bases == 20) {
                saw_actual_geometry = true;
                require(!candidate.tandem.explained, "Graph source lost its contradictory verified tail");
            }
        require(saw_actual_geometry, "Graph witnesses were not keyed by actual verified flank geometry");
        auto shifted_majority = periodic_repeat;
        for (const auto i : {0, 1, 9, 10})
            shifted_majority[i] = shifted_majority[i] == 'A' ? 'C' : 'A';
        const auto right_trained = rr::tandem_evidence(raw_candidate(shifted_majority, periodic_spacer, periodic_repeat));
        require(right_trained.explained && right_trained.training_side == "right" &&
                right_trained.model_orientation == "reverse" && right_trained.training_kind == "observed_full_repeat" &&
                right_trained.left_observed_bases == 27 && right_trained.right_observed_bases == 27,
                "A valid right-trained physical explanation was discarded");
        const auto right_only_full = rr::tandem_evidence(raw_candidate(periodic_repeat, periodic_spacer, periodic_repeat, 20, 27));
        require(right_only_full.explained && right_only_full.training_side == "right" &&
                right_only_full.training_kind == "observed_full_repeat" &&
                right_only_full.left_observed_bases == 27 && right_only_full.right_observed_bases == 20,
                "A graph source did not preserve its physically full right repeat and swapped spans");
        auto deleted_spacer = periodic_spacer;
        deleted_spacer.erase(13, 1);
        const auto deletion_model = rr::tandem_evidence(raw_candidate(periodic_repeat, deleted_spacer, periodic_repeat));
        require(deletion_model.explained && deletion_model.spacer_template_length != deleted_spacer.size(),
                "One spacer indel hid the next-repeat periodic phase");
        require(!rr::tandem_evidence(raw_candidate(periodic_repeat, spacers[2], periodic_repeat)).explained,
                "A periodic repeat alone vetoed a diverse spacer");
        require(rr::tandem_evidence(raw_candidate(repeat, repeat, repeat)).explained,
                "A full repeat-sized spacer escaped the three-copy explanation");
        auto repeat_like = repeat;
        repeat_like[repeat_like.size() / 2] = repeat_like[repeat_like.size() / 2] == 'A' ? 'C' : 'A';
        require(rr::tandem_evidence(raw_candidate(repeat, repeat_like, repeat)).explained,
                "A one-substitution repeat-sized spacer escaped the tandem explanation");
        const std::string mono_repeat(50, 'A');
        require(!rr::tandem_evidence(raw_candidate(mono_repeat, std::string(41, 'A') + std::string(9, 'C'),
                    mono_repeat)).explained &&
                rr::tandem_evidence(raw_candidate(mono_repeat, std::string(42, 'A') + std::string(8, 'C'),
                    mono_repeat)).explained,
                "The tandem similarity threshold is not strictly greater than .82");
        const std::vector<Array> periodic_seeds{{periodic_repeat,
            {{periodic_repeat, periodic_spacer}, {periodic_repeat, periodic_long}}}};
        std::vector<std::string> tandem_reads;
        for (const auto& spacer : {periodic_spacer, periodic_long}) {
            witness(tandem_reads, periodic_repeat + spacer + periodic_repeat);
            witness(tandem_reads, periodic_repeat + spacer + periodic_repeat.substr(0, 16), true);
            tandem_reads.push_back(periodic_repeat.substr(periodic_repeat.size() - 16) + spacer + periodic_repeat.substr(0, 16));
            tandem_reads.push_back(periodic_repeat.substr(periodic_repeat.size() - 17) + spacer + periodic_repeat.substr(0, 17));
        }
        const auto tandem_result = run(tandem_reads, periodic_seeds);
        require(!selected(tandem_result).empty() && promoted(tandem_result).empty(),
                "Supported tandem evidence was erased or automatically promoted");
        for (const auto& candidate : tandem_result.candidates) if (candidate.accepted)
            require(candidate.tandem_ambiguous && !candidate.non_tandem_support.accepted() &&
                    candidate.error_status == "not_evaluated_tandem_ambiguous",
                    "Tandem evidence was credited to the allele-error rejection");
        auto tandem_reversed = tandem_reads;
        for (auto& read : tandem_reversed) read = rl::reverse_complement(read);
        const auto tandem_rc = run(tandem_reversed, periodic_seeds);
        require(inventory(tandem_rc) == inventory(tandem_result) && promoted(tandem_rc).empty(),
                "Whole-read reverse complement changed tandem eligibility");

        // A supported contradictory next repeat can distinguish an exact pair
        // from the tandem explanation. Count only those independent witnesses,
        // never all units minus explained units or the sum across source rows.
        auto contradictory_repeat = periodic_repeat;
        for (std::size_t i = 0; i < contradictory_repeat.size(); i += 4)
            contradictory_repeat[i] = contradictory_repeat[i] == 'A' ? 'C' : 'A';
        const std::vector<Array> contradictory_seeds{{periodic_repeat,
            {{periodic_repeat, spacers[0]}, {contradictory_repeat, spacers[1]}}}};
        std::vector<std::string> mixed_tandem;
        witness(mixed_tandem, periodic_repeat + periodic_spacer + periodic_repeat);
        witness(mixed_tandem, periodic_repeat + spacers[3] + periodic_repeat);
        mixed_tandem.push_back("A" + periodic_repeat + periodic_spacer + contradictory_repeat + "C");
        const auto one_non_tandem = run(mixed_tandem, contradictory_seeds);
        bool saw_one = false;
        for (const auto& candidate : one_non_tandem.candidates)
            if (candidate.accepted && candidate.left_repeat == periodic_repeat && candidate.spacer == periodic_spacer) {
                saw_one = true;
                require(candidate.non_tandem_support.sequence_units == 1 && candidate.tandem_ambiguous &&
                        candidate.error_observed_units == 1 && !candidate.selection_eligible,
                        "Tandem-compatible units strengthened non-tandem support or error K");
            }
        require(saw_one, "Missing inclusive mixed tandem/non-tandem pair");
        mixed_tandem.push_back("GGG" + periodic_repeat + periodic_spacer + contradictory_repeat + "TTT");
        const auto two_non_tandem = run(mixed_tandem, contradictory_seeds);
        bool saw_two = false;
        for (const auto& candidate : two_non_tandem.candidates)
            if (candidate.accepted && candidate.left_repeat == periodic_repeat && candidate.spacer == periodic_spacer) {
                saw_two = true;
                require(candidate.non_tandem_support.sequence_units == 2 && candidate.non_tandem_support.accepted() &&
                        candidate.error_observed_units == 2 && candidate.selection_eligible,
                        "Independent contradictory repeat contexts did not restore pair eligibility");
            }
        require(saw_two, "Missing independently non-tandem pair");

        std::vector<std::string> many;
        for (int s : {2, 3}) for (std::size_t i = 1; i <= 40; ++i)
            many.push_back(std::string(i, 'A') + full(spacers[s], repeat) + std::string(i, 'C'));
        const auto saturated = run(many);
        require(selected(saturated).size() == 2, "Saturation changed accepted inventory");
        for (const auto& candidate : saturated.candidates)
            require(candidate.support.sequence_units == rr::SUPPORT_CAP && candidate.support.start_offsets == rr::SUPPORT_CAP &&
                candidate.support.units_saturated && candidate.support.offsets_saturated, "Evidence memory was not capped explicitly");
        std::shuffle(many.begin(), many.end(), random);
        require(inventory(run(many)) == inventory(saturated), "Saturation depended on first-arriving reads");
        rr::Config bounded;
        bounded.maximum_candidates_per_family = 1;
        const auto limited = run(reads, {}, bounded);
        require(limited.candidates.empty() && limited.counts.at("read_recruitment_resource_limited_families") == 1,
                "Family overload retained arbitrary first candidates");
        std::reverse(reads.begin(), reads.end());
        require(inventory(run(reads, {}, bounded)) == inventory(limited), "Family overload depends on arrival order");
        bounded = rr::Config{}; bounded.maximum_total_candidates = 1;
        const auto global = run(reads, {}, bounded);
        require(global.candidates.empty() && global.counts.at("read_recruitment_global_resource_limit") == 1,
                "Global evidence budget retained first-arriving successes");
        std::reverse(reads.begin(), reads.end());
        require(inventory(run(reads, {}, bounded)) == inventory(global), "Global overload depends on arrival order");
        bounded = rr::Config{}; bounded.maximum_seed_placements_per_read = 1;
        require(selected(run(reads, {}, bounded)).empty(), "Read overload partially credited earlier hits");
        bounded = rr::Config{}; bounded.maximum_read_length = 20;
        require(run(reads, {}, bounded).counts.at("read_recruitment_resource_limited_units") == reads.size(),
                "Overlong reads bypassed resource bound");

        const auto plain = fastq("good.fq", reads);
        { std::ifstream in(plain); std::string text((std::istreambuf_iterator<char>(in)), {});
          gzFile out = gzopen((root / "good.fq.gz").c_str(), "wb");
          require(out != nullptr, "gzip open failed"); gzwrite(out, text.data(), static_cast<unsigned>(text.size())); gzclose(out); }
        require(inventory(rr::infer(arrays, {(root / "good.fq.gz").string()})) == inventory(run(reads)), "Gzip altered recruitment");
        auto throws = [&](const std::vector<std::string>& files) {
            bool failed = false;
            try { rr::infer(arrays, files); } catch (const std::exception&) { failed = true; }
            require(failed, "Malformed or unpaired reads were accepted");
        };
        { std::ofstream out(root / "bad.fq"); out << "@x\nACGT\n+\nII\n"; }
        throws({(root / "bad.fq").string()});
        throws({plain, fastq("short.fq", {reads[0]})});
        { std::ofstream out(root / "wrong.fq"); out << "@wrong\nACGT\n+\nIIII\n"; }
        throws({plain, (root / "wrong.fq").string()});
        { std::ofstream out(root / "reads.fa"); for (std::size_t i = 0; i < reads.size(); ++i) out << '>' << i << '\n' << reads[i] << '\n'; }
        const auto fasta = rr::infer(arrays, {(root / "reads.fa").string()});
        require(fasta.candidates.empty() && fasta.counts.at("read_recruitment_records_without_quality") == reads.size(),
                "FASTA created invented quality support");
    } catch (...) { std::filesystem::remove_all(root); throw; }
    std::filesystem::remove_all(root);
}
