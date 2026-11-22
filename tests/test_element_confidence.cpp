#include "crispr_array/element_confidence.h"
#include <chrono>
#include <filesystem>
#include <iostream>

namespace ec = Zuck::element_confidence;
namespace db = Zuck::array_database;
void require(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}
const std::string repeat = "GTTCACTGCCGTACAGGCAGCTTAGAAA";
const std::string spacer = "ACGTGGCATGACTGATCGTACGATCGTACGAA";
const std::string other = "AGGTCTCGATAGCATGCCGATATCCGTAGTCA";
std::string span(const std::string& r, const std::string& s) { return r + s + r.substr(0, 16); }
db::Record record(std::string id, std::string r, std::vector<std::string> spacers) {
    db::Record result;
    result.id = std::move(id); result.category = "Bona-fide"; result.bona_fide = true;
    result.confidence = 0.9; result.repeats.insert(db::canonical(r));
    for (const auto& s : spacers) result.spacers.insert(db::canonical(s));
    return result;
}
int main() {
    const auto root = std::filesystem::temp_directory_path() /
        ("zuck-confidence-" + std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()));
    std::filesystem::create_directory(root);
    const auto path = root / "reads.fastq";
    const std::vector<std::pair<std::string, std::string>> pairs{{repeat, spacer}};
    const auto context = span(repeat, spacer);
    auto write = [&](const std::vector<std::string>& reads, char quality = 'I') {
        std::ofstream out(path);
        for (std::size_t i = 0; i < reads.size(); ++i)
            out << '@' << i << '\n' << reads[i] << "\n+\n" << std::string(reads[i].size(), quality) << '\n';
    };
    auto run = [&](const std::vector<std::pair<std::string, std::string>>& p,
                   const std::vector<db::Record>& records = std::vector<db::Record>{}) {
        return ec::annotate(p, {path.string()}, records);
    };
    write({context + "AA", "CCC" + context});
    auto result = run(pairs);
    require(result.rows[0].context.accepted(), "Two quality spans at different offsets must qualify");
    require(result.rows[0].context_score == 0.5 && result.rows[0].repeat_score == 0.5 &&
            result.rows[0].spacer_score == 0.5, "Read-presence indices must be explicit saturated counts");
    require(std::isnan(result.rows[0].boundary_score) && result.rows[0].boundary_status == "unresolved",
            "A unique emitted interpretation must not imply a certain boundary");
    write({db::reverse_complement(context) + "AA", "CCC" + db::reverse_complement(context)});
    require(run(pairs).rows[0].context.accepted(), "Whole-span reverse complement must support presence");
    const auto independently_reversed = std::vector<std::pair<std::string, std::string>>{
        {db::reverse_complement(repeat), db::reverse_complement(spacer)}};
    write({context + "AA", "CCC" + context});
    require(!run(independently_reversed).rows[0].context.accepted(),
            "Independent component reverse complements are not whole-span adjacency");
    write({context, context, context});
    require(run(pairs).rows[0].context.sequence_units == 1, "Duplicate read strings cannot inflate unit counts");
    write({context + "AA", context + "CC"});
    require(!run(pairs).rows[0].context.accepted(), "Different sequence units need different offsets");
    write({context + context, context + context});
    require(!run(pairs).rows[0].context.accepted(), "Multiple positions in one unique read are not independent units");
    write({context + "AA", "CCC" + context}, '!');
    require(run(pairs).rows[0].context_score == 0, "Low-quality exact spans cannot support the sequence");
    write({context + "AA", "CCC" + context}, '+');
    require(run(pairs).rows[0].context_score == 0, "Q10 bases still need the total expected-error bound");
    const std::string threshold_repeat = repeat + std::string(30 - repeat.size(), 'C');
    for (std::size_t n : {53, 54, 55}) {
        const std::string threshold_spacer = spacer + std::string(n - spacer.size(), 'G');
        const auto threshold_span = span(threshold_repeat, threshold_spacer);
        const std::vector<std::pair<std::string, std::string>> threshold_pairs{{threshold_repeat, threshold_spacer}};
        write({threshold_span + "AA", std::string(50, 'T') + threshold_span}, '5');
        double direct_errors = 0;
        for (std::size_t i = 0; i < threshold_span.size(); ++i)
            direct_errors += std::pow(10.0, -20.0 / 10.0);
        const bool expected = direct_errors <= 1.0;
        require(Zuck::database_reads::validate(threshold_pairs, {path.string()}).support[0].accepted() == expected,
                "Database quality spans must preserve the direct-summation threshold");
        require(run(threshold_pairs).rows[0].context.accepted() == expected,
                "Confidence decisions below, at and above 1 must agree with independent direct summation");
    }
    write({repeat + spacer + "AA", "CCC" + repeat + spacer});
    require(run(pairs).rows[0].spacer.accepted() && !run(pairs).rows[0].context.accepted(),
            "Spacer presence must not be confused with next-repeat adjacency");
    write({repeat, "CCC" + repeat});
    require(run(pairs).rows[0].context_observability == "span_exceeds_max_quality_read" &&
            std::isnan(run(pairs).rows[0].context_score), "Too-short reads produce unknown context evidence");
    { std::ofstream out(path); out << ">a\n" << context << "\n>b\nCCC" << context << '\n'; }
    require(run(pairs).rows[0].context_observability == "no_quality" &&
            std::isnan(run(pairs).rows[0].repeat_score), "FASTA is quality-unassessed, not negative support");
    require(ec::annotate(pairs, {}).rows[0].context_observability == "no_reads" &&
            std::isnan(ec::annotate(pairs, {}).rows[0].spacer_score), "No reads produce missing scores");

    // Real repeat copies can differ: requiring R1-S-prefix(R1) rejects an
    // observed R1-S-prefix(R2). Alternate anchors must come from a compatible
    // observed row in this family, not a different family or an arbitrary SNP.
    struct FamilyArray { std::vector<std::pair<std::string, std::string>> entries; };
    std::string next_repeat = repeat;
    next_repeat[5] = next_repeat[5] == 'A' ? 'C' : 'A';
    const std::vector<std::pair<std::string, std::string>> family_pairs{
        {repeat, spacer}, {next_repeat, other}};
    const auto family_prefixes = ec::family_repeat_prefixes(std::vector<FamilyArray>{{family_pairs}});
    const auto mixed_span = repeat + spacer + next_repeat.substr(0, 16);
    write({mixed_span + "AA", "CCC" + mixed_span});
    require(!run(family_pairs).rows[0].context.accepted(), "Baseline must isolate next-repeat allele loss");
    auto family_result = ec::annotate(family_pairs, {path.string()}, {}, family_prefixes);
    require(family_result.rows[0].context.accepted() &&
            family_result.rows[0].context.sequence_units == 2,
            "Compatible observed next-repeat allele must support literal context");
    const auto restricted_prefixes = ec::family_repeat_prefixes(
        std::vector<FamilyArray>{{family_pairs}}, {{0, 2}});
    const auto original_prefixes = ec::family_repeat_prefixes(
        std::vector<FamilyArray>{{{{repeat, spacer}}}});
    require(restricted_prefixes.size() == family_pairs.size() &&
            restricted_prefixes[0] == original_prefixes[0] &&
            !ec::annotate(family_pairs, {path.string()}, {}, restricted_prefixes).rows[0].context.accepted(),
            "A recruited repeat allele must not lend a new prefix to an original row");
    write({context + "AA", "CCC" + context, mixed_span + "CC", "GGG" + mixed_span});
    const auto original_annotation = ec::annotate(pairs, {path.string()}, {}, original_prefixes);
    const auto restricted_annotation = ec::annotate(family_pairs, {path.string()}, {}, restricted_prefixes);
    require(original_annotation.rows[0].context.sequence_units == restricted_annotation.rows[0].context.sequence_units &&
            original_annotation.rows[0].context.start_offsets == restricted_annotation.rows[0].context.start_offsets &&
            original_annotation.rows[0].context.accepted() == restricted_annotation.rows[0].context.accepted(),
            "Extra evidence-only patterns must preserve original literal support exactly");
    bool invalid_exclusion = false;
    try { ec::family_repeat_prefixes(std::vector<FamilyArray>{{family_pairs}}, {{0, 3}}); }
    catch (const std::invalid_argument&) { invalid_exclusion = true; }
    require(invalid_exclusion, "Out-of-range family context exclusions must fail closed");
    write({mixed_span + "AA", "CCC" + mixed_span});
    const auto split_prefixes = ec::family_repeat_prefixes(std::vector<FamilyArray>{
        {{{repeat, spacer}}}, {{{next_repeat, other}}}});
    require(!ec::annotate(family_pairs, {path.string()}, {}, split_prefixes).rows[0].context.accepted(),
            "Different families must not lend next-repeat alleles");
    auto divergent = next_repeat;
    for (std::size_t p : {0, 1, 2, 3, 4}) divergent[p] = divergent[p] == 'A' ? 'C' : 'A';
    const auto divergent_prefixes = ec::family_repeat_prefixes(std::vector<FamilyArray>{
        {{{repeat, spacer}, {divergent, other}}}});
    require(divergent_prefixes[0].size() == 1, "Distant repeats within a content group cannot lend anchors");
    const auto reverse_allele_prefixes = ec::family_repeat_prefixes(std::vector<FamilyArray>{
        {{{repeat, spacer}, {db::reverse_complement(next_repeat), other}}}});
    require(reverse_allele_prefixes[0] == family_prefixes[0],
            "Reverse-strand family rows must contribute the correctly oriented repeat prefix");
    write({db::reverse_complement(mixed_span) + "AA", "CCC" + db::reverse_complement(mixed_span)});
    require(ec::annotate(family_pairs, {path.string()}, {}, family_prefixes).rows[0].context.accepted(),
            "Whole mixed-context reverse complement must qualify");
    write({mixed_span + "AA", "CCC" + mixed_span}, '!');
    require(!ec::annotate(family_pairs, {path.string()}, {}, family_prefixes).rows[0].context.accepted(),
            "Family context does not bypass quality filtering");
    write({mixed_span, mixed_span});
    require(ec::annotate(family_pairs, {path.string()}, {}, family_prefixes).rows[0].context.sequence_units == 1,
            "Family anchors do not multiply duplicate read evidence");
    write({mixed_span + "AA", context + "CC"});
    require(!ec::annotate(family_pairs, {path.string()}, {}, family_prefixes).rows[0].context.accepted(),
            "Support across different anchors still requires distinct offsets");
    write({mixed_span + "AA", "CCC" + context});
    require(ec::annotate(family_pairs, {path.string()}, {}, family_prefixes).rows[0].context.accepted(),
            "Distinct reads across observed next-repeat alleles may corroborate one spacer context");
    write({repeat + spacer + "AAAAAAAAAAAAAAAAAA", "CCC" + repeat + spacer + "AAAAAAAAAAAAAAAA"});
    require(!ec::annotate(family_pairs, {path.string()}, {}, family_prefixes).rows[0].context.accepted(),
            "Unobserved next-repeat anchors must not qualify");
    require(ec::annotate(family_pairs, {}, {}, family_prefixes).rows[0].context_observability == "no_reads",
            "Family mode preserves explicit evidence missingness");

    std::vector<std::string> many;
    for (std::size_t i = 0; i < 40; ++i) many.push_back(std::string(i, 'C') + context + "GG");
    write(many); const auto capped = run(pairs);
    require(capped.rows[0].context.sequence_units == ec::SUPPORT_CAP &&
            capped.rows[0].context.start_offsets == ec::SUPPORT_CAP &&
            capped.rows[0].context.units_saturated && capped.rows[0].context.offsets_saturated,
            "Memory-bounded counts must expose inclusive saturation");
    std::reverse(many.begin(), many.end()); write(many);
    const auto reversed = run(pairs);
    require(reversed.rows[0].context_score == capped.rows[0].context_score &&
            reversed.rows[0].context.sequence_units == capped.rows[0].context.sequence_units,
            "Count caps and scores must be read-order invariant");

    const auto shorter = repeat.substr(0, repeat.size() - 1);
    const auto shifted = repeat.substr(repeat.size() - 1) + spacer;
    const auto shifted_other = repeat.substr(repeat.size() - 1) + other;
    const std::vector<std::pair<std::string, std::string>> alternatives = {
        {repeat, spacer}, {shorter, shifted}, {repeat, other}, {shorter, shifted_other}};
    require(span(shorter, shifted) == context, "Fixture must have literally identical spans under shifted cuts");
    require(ec::pair_identity({repeat, spacer}) != ec::pair_identity({shorter, shifted}),
            "Different cut positions remain different identities");
    require(ec::circular_identity(repeat + spacer) == ec::circular_identity(shorter + shifted),
            "Boundary alternatives must compete within the same circular core");
    require(ec::circular_identity(repeat + spacer) ==
            ec::circular_identity(db::reverse_complement(repeat + spacer)),
            "Competition cores must be strand invariant");
    write({context + "AA", "CCC" + context, span(repeat, other) + "AA", "CCC" + span(repeat, other)});
    const auto unscored = run(alternatives);
    require(unscored.rows[0].context_score == unscored.rows[1].context_score &&
            unscored.rows[0].boundary_alternatives == 2 &&
            unscored.rows[0].boundary_status == "ambiguous_unscored",
            "Indistinguishable read spans must expose cut ambiguity");
    auto source = record("source", repeat, {spacer, other});
    const auto preferred = run(alternatives, {source});
    require(preferred.rows[0].database_other_spacers == 1 && preferred.rows[0].database_other_cores == 1 &&
            std::abs(preferred.rows[0].database_score - 0.3) < 1e-12,
            "Same-array evidence excludes the query's own spacer and circular core");
    require(preferred.rows[0].boundary_status == "catalogue_preferred" &&
            preferred.rows[1].boundary_status == "alternative_to_catalogue_cut",
            "Only independent catalogue corroboration may provide boundary preference");
    auto copied = source; copied.id = "source_copy";
    const auto duplicate = run(alternatives, {source, copied});
    require(duplicate.rows[0].database_score == preferred.rows[0].database_score,
            "Copied source arrays cannot accumulate confidence");
    const auto two_sources = run(alternatives, {source, record("alternative", shorter, {shifted, shifted_other})});
    require(two_sources.rows[0].boundary_status == "catalogue_ambiguous" &&
            two_sources.rows[0].boundary_database_alternatives == 2 &&
            two_sources.rows[0].boundary_score < preferred.rows[0].boundary_score,
            "Two supported indistinguishable cuts must reduce specificity and remain ambiguous");
    require(run(alternatives, {record("self_only", repeat, {spacer})}).rows[0].database_score == 0,
            "One source spacer cannot vote for its own cut");
    write({context + "AA", "CCC" + context});
    require(run(alternatives, {source}).rows[0].database_score == 0,
            "Unobserved source spacers cannot provide leave-one-out votes");
    const auto shifted_only = std::vector<std::pair<std::string, std::string>>{{repeat, spacer}, {shorter, shifted}};
    auto mixed_source = record("mixed", repeat, {spacer, shifted});
    mixed_source.repeats.insert(db::canonical(shorter));
    require(run(shifted_only, {mixed_source}).rows[0].database_score == 0,
            "Alternate splits of one circular core cannot corroborate each other");

    { std::ofstream out(path); out << "@broken\nACGT\n+\nII\n"; }
    bool malformed = false;
    try { run(pairs); } catch (const std::runtime_error&) { malformed = true; }
    require(malformed, "Malformed FASTQ must fail confidence scoring");
    struct Array { std::vector<std::pair<std::string, std::string>> entries; };
    ec::write_tables(root.string(), std::vector<Array>{{pairs}}, ec::annotate(pairs, {}));
    std::ifstream table(root / "element_confidence.tsv");
    std::string header, row; std::getline(table, header); std::getline(table, row);
    require(header.find("boundary_status") != std::string::npos && row.find("\tNA\t") != std::string::npos,
            "Unknown evidence must serialize explicitly as NA");
    std::filesystem::remove_all(root);
    std::cout << "Element confidence presence, quality, duplication, caps, observability, leave-one-out, and boundary controls passed\n";
}
