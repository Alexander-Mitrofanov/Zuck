#include "crispr_array/supported_calls.h"
#include "crispr_array/read_recruitment_output.h"
#include <chrono>
#include <filesystem>
#include <iostream>
#include <sstream>

namespace sc = Zuck::supported_calls;
namespace ec = Zuck::element_confidence;
namespace db = Zuck::array_database;

struct Array {
    std::string consensus;
    std::vector<std::pair<std::string, std::string>> entries;
};
void require(bool ok, const char* message) {
    if (!ok) throw std::runtime_error(message);
}
template <class Fn> void rejects(Fn function, const char* message) {
    bool rejected = false;
    try { function(); } catch (const std::invalid_argument&) { rejected = true; }
    require(rejected, message);
}
std::string sequence(std::size_t length) {
    std::string result;
    while (result.size() < length) result += "ACGT";
    return result.substr(0, length);
}
ec::Row evidence(std::size_t units = 2, std::size_t offsets = 2,
                 const std::string& observability = "observable") {
    ec::Row row;
    row.context.sequence_units = units;
    row.context.start_offsets = offsets;
    row.context.observable = observability == "observable";
    row.context_observability = observability;
    return row;
}
std::string contents(const std::filesystem::path& path) {
    std::ifstream in(path);
    require(static_cast<bool>(in), "Expected output file");
    return {std::istreambuf_iterator<char>(in), std::istreambuf_iterator<char>()};
}
std::map<std::string, std::string> files(const std::filesystem::path& directory) {
    std::map<std::string, std::string> result;
    for (const auto& file : std::filesystem::directory_iterator(directory))
        result.emplace(file.path().filename().string(), contents(file.path()));
    return result;
}
struct Temporary {
    std::filesystem::path path = std::filesystem::temp_directory_path() /
        ("zuck-supported-" + std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()));
    Temporary() { std::filesystem::create_directory(path); }
    ~Temporary() { std::error_code ignored; std::filesystem::remove_all(path, ignored); }
};

void recruitment_graph_policy() {
    namespace rr = Zuck::read_recruitment;
    namespace ro = Zuck::read_recruitment_output;
    struct Limits { bool pair_allowed(const std::string&, const std::string&) const { return true; } } limits;
    const std::string repeat = "GTTCACTGCCGTACAGGCAGCTTAGAAA";
    const std::string a(30, 'A'), c(30, 'C'), novel = std::string(15, 'A') + std::string(15, 'C');
    auto candidate = [&](std::string r, std::string s, std::set<std::size_t> owners) {
        rr::Candidate x;
        x.left_repeat = x.right_repeat = std::move(r); x.spacer = std::move(s);
        x.family_owners = std::move(owners); x.pair_only = x.accepted = x.selection_eligible = true;
        x.support = x.non_tandem_support = {2, 2, false, false};
        x.distinct_supported_family_spacers = 2; x.evidence_kind = "full_pair";
        x.tandem_status = "independent_non_tandem_support"; x.error_status = "no_supported_neighbor";
        return x;
    };
    auto apply = [&](std::vector<Array>& arrays, rr::Result result, Temporary& out,
                     std::map<std::string, std::size_t>& counts) {
        return ro::apply(arrays, result, limits, out.path, counts);
    };

    std::vector<Array> arrays{{repeat, {{repeat, a}, {repeat, c}, {repeat, a}}}};
    rr::Result result;
    result.candidates = {candidate(repeat, a, {0}), candidate(repeat, c, {0}),
                         candidate(repeat, novel, {0}), candidate(db::reverse_complement(repeat),
                         db::reverse_complement(a), {0})};
    result.candidates[3].graph_proposed = true; // RC alias is NOT an original exact row.
    Temporary out;
    std::map<std::string, std::size_t> counts;
    const auto proof = apply(arrays, result, out, counts);
    require(proof.at({0, 1}).accepted() && proof.at({0, 2}).accepted(),
            "Original graph rows with fixed eligible diverse context must remain promotable");
    require(counts.at("read_recruitment_original_snapshot_rows") == 3 &&
            counts.at("read_recruitment_added_rows") == 2 &&
            counts.at("read_recruitment_novel_boundary_unresolved_rows") == 2,
            "Novel and RC-alias pairs must remain inclusive without original membership");
    for (const auto& item : proof) if (item.second.added_row)
        require(!item.second.accepted() && item.second.status() == "novel_boundary_unresolved",
                "Novel graph-compatible alleles cannot become selected through the new channel");
    const auto snapshot = contents(out.path / "read_recruitment_original_rows.tsv");
    require(snapshot.find("Array_1\t3\t" + repeat + "\t" + a + "\n") != std::string::npos &&
            snapshot.find(novel) == std::string::npos,
            "Original snapshot must preserve duplicate row identities and exclude appended observations");
    const auto contexts = contents(out.path / "read_recruitment_context.tsv");
    require(contexts.find("1\tArray_1\t2\t2\t0\t0\t30\t30\t0\tsupported\n") != std::string::npos &&
            contexts.find("2\tArray_1\t1\t1\t") != std::string::npos,
            "Context must bind deterministic candidate and first exact original row IDs");
    ec::Result all_literal;
    for (std::size_t i = 0; i < arrays[0].entries.size(); ++i) all_literal.rows.push_back(evidence());
    const auto calls = sc::infer(arrays, all_literal, "standard", proof, true);
    require(calls.rows[0].selected && calls.rows[1].selected && calls.rows[2].selected,
            "Literal original rows including duplicates must remain selected");
    for (std::size_t i = 3; i < calls.rows.size(); ++i)
        require(!calls.rows[i].selected && calls.rows[i].context_status == "supported",
                "Novel rows may not bypass the policy through downstream literal context");
    auto forged = proof;
    for (auto& item : forged) if (item.second.added_row) item.second.selection_eligible = true;
    rejects([&] { sc::infer(arrays, all_literal, "standard", forged, true); },
            "An appended row marked promotable must fail closed before selection");

    // An attractive novel spacer cannot supply the original-row context prior.
    arrays = {{repeat, {{repeat, a}}}};
    result.candidates = {candidate(repeat, a, {0}), candidate(repeat, c, {0})};
    Temporary novel_context;
    counts.clear();
    const auto no_context = apply(arrays, result, novel_context, counts);
    require(no_context.at({0, 1}).status() == "insufficient_spacer_diversity" &&
            !no_context.at({0, 1}).accepted(), "Novel context must not license an original focal row");
    ec::Result original_literal; original_literal.rows = {evidence(), evidence()};
    const auto kept = sc::infer(arrays, original_literal, "standard", no_context, true);
    require(kept.rows[0].selected && !kept.rows[1].selected,
            "Missing policy diversity must preserve original literal calls only");

    // Both declared owners must supply context, even when only one row is emitted.
    arrays = {{repeat, {{repeat, a}, {repeat, c}}}, {repeat, {{repeat, a}}}};
    result.candidates = {candidate(repeat, a, {0, 1}), candidate(repeat, c, {0})};
    Temporary owners;
    counts.clear();
    const auto missing_owner = apply(arrays, result, owners, counts);
    require(!missing_owner.at({0, 1}).accepted() && !missing_owner.at({1, 1}).accepted() &&
            contents(owners.path / "read_recruitment_context.tsv").find(
                "1\tArray_2\t0\t0\t-1\t0\t-1\t0\t0\tmissing\n") != std::string::npos,
            "A supported owner cannot erase another declared owner's missing context");

    // Context is fixed before final gating: its own second owner may be missing.
    arrays = {{repeat, {{repeat, a}, {repeat, c}}}, {repeat, {{repeat, c}}}};
    result.candidates = {candidate(repeat, a, {0}), candidate(repeat, c, {0, 1})};
    Temporary fixed;
    counts.clear();
    const auto fixed_context = apply(arrays, result, fixed, counts);
    require(fixed_context.at({0, 1}).accepted() && !fixed_context.at({0, 2}).accepted(),
            "Iterative pruning must not remove an evidence-eligible original context source");
    result.candidates[1].selection_eligible = false;
    result.candidates[1].error_ambiguous = true;
    arrays = {{repeat, {{repeat, a}, {repeat, c}}}, {repeat, {{repeat, c}}}};
    Temporary ineligible;
    counts.clear();
    require(!apply(arrays, result, ineligible, counts).at({0, 1}).accepted(),
            "V3 evidence-ineligible contexts cannot license diversity");

    // Boundary .18 is inclusive; reverse orientations are compared directly.
    auto focal = candidate(repeat, std::string(50, 'A'), {0});
    auto boundary = candidate(repeat, std::string(41, 'A') + std::string(9, 'C'), {0});
    require(ro::context_witness(focal, boundary, 2, 1).candidate_id == 2,
            "Exact .18 spacer dissimilarity must pass");
    boundary.spacer = std::string(42, 'A') + std::string(8, 'C');
    require(!ro::context_witness(focal, boundary, 2, 1).candidate_id,
            "Subthreshold spacer dissimilarity must not pass");
    boundary.left_repeat = db::reverse_complement(repeat);
    boundary.spacer = std::string(50, 'C');
    require(ro::context_witness(focal, boundary, 2, 1).repeat_reverse_complement,
            "Repeat neighbor restriction must allow RC with explicit witness orientation");
    boundary.spacer = std::string(50, 'T');
    require(!ro::context_witness(focal, boundary, 2, 1).candidate_id,
            "A reverse-complement spacer cannot manufacture diversity");
    boundary.spacer = std::string(50, 'C'); boundary.left_repeat += 'A';
    require(!ro::context_witness(focal, boundary, 2, 1).candidate_id,
            "Repeat context must have exactly equal length");
    boundary.left_repeat = repeat; boundary.left_repeat[0] = 'A'; boundary.left_repeat[1] = 'A';
    require(!ro::context_witness(focal, boundary, 2, 1).candidate_id,
            "Two-repeat-mismatch context must not lend support");

    arrays = {{repeat, {{repeat, a}}}};
    Temporary empty;
    counts.clear();
    require(apply(arrays, {}, empty, counts).empty() &&
            counts.at("read_recruitment_original_snapshot_rows") == 1 &&
            counts.at("read_recruitment_context_rows") == 0,
            "Empty recruitment must still export the full original snapshot and explicit empty context table");
}

int main() {
    recruitment_graph_policy();
    const std::string repeat = "GTTCACTGCCGTACAGGCAGCTTAGAAA";
    const std::string alternative_repeat = "ATTCACTGCCGTACAGGCAGCTTAGAAA";
    std::vector<Array> arrays{{repeat, {{repeat, sequence(25)}, {repeat, sequence(26)},
                                      {repeat, sequence(50)}, {repeat, sequence(51)}}},
                              {alternative_repeat, {{alternative_repeat, sequence(26)},
                                  {db::reverse_complement(repeat), db::reverse_complement(sequence(26))}}}};
    ec::Result confidence;
    confidence.rows = {evidence(), evidence(), evidence(), evidence(), evidence(1, 2), evidence()};
    const auto selected = sc::infer(arrays, confidence, "standard");
    require(selected.rows.size() == 6 && selected.rows[1].selected && selected.rows[2].selected,
            "Inclusive profile boundaries 26 and 50 must be accepted");
    require(!selected.rows[0].selected && !selected.rows[3].selected &&
            selected.rows[0].reason == "spacer_length_outside_profile" &&
            selected.rows[3].context_status == "supported",
            "Lengths 25 and 51 are out-of-profile, not unsupported sequence evidence");
    require(!selected.rows[4].selected && selected.rows[5].selected && selected.spacers.size() == 2,
            "Any-row inventory support must not license another repeat pairing");
    require(selected.rows[4].reason == "insufficient_context_units" &&
            selected.rows[4].array_index == 1 && selected.rows[4].entry_index == 1,
            "Decisions must explain eligibility and retain original coordinates");
    require(selected.groups.size() == 2 && selected.groups[1].entry_indices == std::vector<std::size_t>{2} &&
            selected.groups[1].canonical_spacers.size() == 1 &&
            selected.counts.at("supported_singleton_groups") == 1,
            "Singleton groups survive and cannot borrow rows from the same spacer elsewhere");
    require(selected.counts.at("supported_selected_rows") == 3 &&
            selected.counts.at("supported_unselected_rows") == 3 &&
            selected.counts.at("supported_spacer_types") == 2,
            "Counts distinguish oriented rows from canonical spacer types");

    std::vector<Array> prune_arrays{
        {repeat, {{repeat, sequence(30)}, {repeat, sequence(30)}}},
        {alternative_repeat, {{alternative_repeat, sequence(31)}}},
        {repeat, {{repeat, sequence(32)}, {repeat, sequence(33)}}},
    };
    ec::Result prune_confidence;
    prune_confidence.rows = {evidence(2, 2), evidence(2, 2), evidence(3, 2),
                             evidence(2, 2), evidence(2, 2)};
    const auto pruned = sc::infer(
        prune_arrays, prune_confidence, "standard", {}, false, true);
    require(!pruned.rows[0].selected && !pruned.rows[1].selected &&
            pruned.rows[0].reason == "singleton_minimal_literal_context_pruned" &&
            pruned.counts.at("supported_singleton_precision_pruned_rows") == 2 &&
            pruned.counts.at("supported_singleton_precision_pruned_types") == 1,
            "Every weak literal occurrence of a canonical singleton type must be demoted together");
    require(pruned.rows[2].selected && pruned.rows[3].selected && pruned.rows[4].selected &&
            pruned.groups.size() == 2 && pruned.spacers.size() == 3,
            "Strong-in-one-dimension singletons and weak rows in multi-spacer families must be preserved");
    Temporary pruned_output;
    sc::write_tables(pruned_output.path.string(), prune_arrays, pruned);
    require(contents(pruned_output.path / "singleton_precision_pruned.tsv").find(
                "singleton_minimal_literal_context_pruned") != std::string::npos &&
            files(pruned_output.path).size() == 6,
            "An enabled precision prune must retain an explicit demotion ledger");

    const std::string circular_spacer = "TGCATGACCTAGGCTAACCGTTAGGATCCA";
    const std::string circular_unit = repeat + circular_spacer;
    const std::string shifted_unit = circular_unit.substr(1) + circular_unit.substr(0, 1);
    const std::string shifted_repeat = shifted_unit.substr(0, repeat.size());
    const std::string shifted_spacer = shifted_unit.substr(repeat.size());
    require(ec::circular_identity(circular_unit) == ec::circular_identity(shifted_unit) &&
            db::canonical(circular_spacer) != db::canonical(shifted_spacer),
            "Fixture must encode two spacer partitions of one exact circular unit");
    std::vector<Array> circular_arrays{
        {repeat, {{repeat, circular_spacer}}},
        {shifted_repeat, {{shifted_repeat, shifted_spacer}}},
    };
    ec::Result circular_confidence;
    circular_confidence.rows = {evidence(3, 3), evidence(6, 5)};
    const auto circular = sc::infer(
        circular_arrays, circular_confidence, "standard", {}, false, false, true);
    require(!circular.rows[0].selected && circular.rows[1].selected &&
            circular.rows[0].reason == "circular_unit_alternative_partition_pruned" &&
            circular.spacers.size() == 1 && circular.groups.size() == 1 &&
            circular.counts.at("supported_circular_unit_duplicate_units") == 1 &&
            circular.counts.at("supported_circular_unit_pruned_types") == 1,
            "One strongest supported partition must represent an exact circular unit");
    std::reverse(circular_arrays.begin(), circular_arrays.end());
    std::reverse(circular_confidence.rows.begin(), circular_confidence.rows.end());
    const auto circular_reordered = sc::infer(
        circular_arrays, circular_confidence, "standard", {}, false, false, true);
    require(circular_reordered.spacers == circular.spacers,
            "Circular-unit selection must be invariant to array order");

    ec::Result tied_confidence;
    tied_confidence.rows = {evidence(2, 2), evidence(2, 2)};
    const auto tied = sc::infer(
        std::vector<Array>{{repeat, {{repeat, circular_spacer}}},
                           {shifted_repeat, {{shifted_repeat, shifted_spacer}}}},
        tied_confidence, "standard", {}, false, false, true);
    require(tied.spacers.size() == 1 &&
            tied.spacers.count(std::min(db::canonical(circular_spacer),
                                        db::canonical(shifted_spacer))) == 1,
            "Evidence ties must use the lexicographically smallest canonical spacer");

    const auto singleton_before_circular = sc::infer(
        std::vector<Array>{{repeat, {{repeat, circular_spacer}}},
                           {shifted_repeat, {{shifted_repeat, shifted_spacer}}}},
        tied_confidence, "standard", {}, false, true, true);
    require(singleton_before_circular.spacers.empty() &&
            singleton_before_circular.counts.at(
                "supported_singleton_precision_pruned_types") == 2 &&
            singleton_before_circular.counts.at(
                "supported_circular_unit_duplicate_units") == 0 &&
            singleton_before_circular.circular_unit_decisions.empty(),
            "Singleton pruning must precede circular-unit canonicalization");

    ec::Result literal_confidence;
    literal_confidence.rows = {evidence(2, 2), evidence(1, 1)};
    sc::RecruitmentEvidence circular_recruitment{
        {{1, 1}, {1, 10, 10, "full_pair", true, false}}};
    const auto literal_wins = sc::infer(
        std::vector<Array>{{repeat, {{repeat, circular_spacer}}},
                           {shifted_repeat, {{shifted_repeat, shifted_spacer}}}},
        literal_confidence, "standard", circular_recruitment, true, false, true);
    require(literal_wins.rows[0].selected && !literal_wins.rows[1].selected,
            "Literal owner evidence must outrank recruitment-only support");

    const std::string other_repeat = "AACCGGTTAACCGGTTAACCGGTTAACC";
    std::vector<Array> protected_arrays{
        {repeat, {{repeat, circular_spacer}}},
        {shifted_repeat, {{shifted_repeat, shifted_spacer}}},
        {other_repeat, {{other_repeat, circular_spacer}}},
    };
    ec::Result protected_confidence;
    protected_confidence.rows = {evidence(2, 2), evidence(6, 6), evidence(2, 2)};
    const auto protected_type = sc::infer(
        protected_arrays, protected_confidence, "standard", {}, false, false, true);
    require(protected_type.rows[0].selected && protected_type.rows[1].selected &&
            protected_type.rows[2].selected &&
            protected_type.counts.at("supported_circular_unit_pruned_types") == 0 &&
            protected_type.counts.at("supported_circular_unit_protected_alternatives") == 1,
            "A canonical type supported under another unit must fail closed and remain selected");

    Temporary circular_output;
    sc::write_tables(circular_output.path.string(),
                     std::vector<Array>{{repeat, {{repeat, circular_spacer}}},
                                        {shifted_repeat, {{shifted_repeat, shifted_spacer}}}},
                     circular);
    require(contents(circular_output.path / "circular_unit_dedup.tsv").find(
                "alternative_partition_pruned") != std::string::npos &&
            files(circular_output.path).size() == 6,
            "Circular-unit decisions must be retained in a complete audit ledger");
    const auto unrestricted = sc::infer(arrays, confidence, "context");
    require(unrestricted.rows[0].selected && unrestricted.rows[3].selected &&
            unrestricted.counts.at("supported_selected_rows") == 5,
            "The context profile imposes no reporting-length restriction");

    sc::RecruitmentEvidence recruited{
        {{0, 2}, {1, 2, 2, "full_partial", false, false}},
        {{1, 1}, {2, 2, 2, "graph_junction", true, false}}};
    auto recovered = sc::infer(arrays, confidence, "standard", recruited, true);
    require(recovered.rows[1].selected && recovered.rows[1].selected_evidence == "literal_context" &&
            recovered.rows[4].selected && recovered.rows[4].selected_evidence == "read_recruitment",
            "Ambiguous recruitment preserves original literal calls; eligible junctions rescue weak context");
    recruited.at({0, 2}).added_row = true;
    recovered = sc::infer(arrays, confidence, "standard", recruited, true);
    require(!recovered.rows[1].selected && recovered.rows[1].context_status == "supported" &&
            recovered.rows[1].reason == "recruitment_error_ambiguous",
            "New ambiguous rows remain inclusive and cannot bypass their gate through literal context");
    recruited.at({0, 2}).selection_block_reason = "tandem_ambiguous";
    recovered = sc::infer(arrays, confidence, "standard", recruited, true);
    require(!recovered.rows[1].selected && recovered.rows[1].context_status == "supported" &&
            recovered.rows[1].reason == "recruitment_tandem_ambiguous" &&
            recovered.rows[1].recruitment.status() == "tandem_ambiguous",
            "Tandem ambiguity retains its reason and cannot bypass the new-row gate");
    recruited.at({1, 1}).selection_eligible = false;
    require(!sc::infer(arrays, confidence, "standard", recruited, true).rows[4].selected,
            "Ambiguous recruitment cannot promote a previously unselected graph row");
    recruited.emplace(std::make_pair(5, 99), sc::RecruitmentSupport{3, 2, 2, "full_pair"});
    rejects([&] { sc::infer(arrays, confidence, "standard", recruited, true); },
            "Unbound recruitment row identities must fail closed");

    std::vector<Array> duplicated{{repeat, {{repeat, sequence(26)}}}};
    ec::Result duplicate_confidence;
    duplicate_confidence.rows = {evidence()};
    const auto one = sc::infer(duplicated, duplicate_confidence, "standard");
    duplicated[0].entries.push_back({repeat, sequence(26)});
    duplicated.push_back({db::reverse_complement(repeat),
                         {{db::reverse_complement(repeat), db::reverse_complement(sequence(26))}}});
    duplicate_confidence.rows = {evidence(), evidence(), evidence()};
    const auto three = sc::infer(duplicated, duplicate_confidence, "standard");
    require(one.spacers == three.spacers && three.counts.at("supported_selected_rows") == 3,
            "Duplicate and reverse-complement rows cannot inflate canonical inventory");

    std::vector<Array> unknown{{repeat, {{repeat, sequence(26)}, {repeat, sequence(27)},
                                       {repeat, sequence(28)}}}};
    ec::Result missing;
    missing.rows = {evidence(0, 0, "no_reads"), evidence(0, 0, "no_quality"),
                    evidence(0, 0, "span_exceeds_max_quality_read")};
    const auto undecided = sc::infer(unknown, missing, "standard");
    require(undecided.spacers.empty() && undecided.groups.empty() &&
            undecided.counts.at("supported_unknown_context_rows") == 3 &&
            undecided.counts.at("supported_insufficient_context_rows") == 0,
            "Missing evidence is unknown, not negative support");
    for (std::size_t i = 0; i < missing.rows.size(); ++i)
        require(undecided.rows[i].context_status == "unknown" &&
                undecided.rows[i].reason == missing.rows[i].context_observability,
                "Unknown evidence retains its specific reason");
    missing.rows = {evidence(2, 1), evidence(1, 1), evidence(0, 0)};
    const auto insufficient = sc::infer(unknown, missing, "standard");
    require(insufficient.rows[0].reason == "insufficient_context_offsets" &&
            insufficient.rows[1].reason == "insufficient_context_units;insufficient_context_offsets",
            "Both independent support dimensions are required with deterministic reasons");
    unknown[0].entries[0].second = sequence(25);
    missing.rows[0] = evidence(0, 0, "no_reads");
    const auto compound = sc::infer(unknown, missing, "standard");
    require(compound.rows[0].context_status == "unknown" &&
            compound.rows[0].reason == "spacer_length_outside_profile;no_reads",
            "Length exclusions must not erase evidence missingness");

    rejects([&] { sc::infer(arrays, confidence, "best"); }, "Unknown profiles must fail");
    auto short_confidence = confidence;
    short_confidence.rows.pop_back();
    rejects([&] { sc::infer(arrays, short_confidence, "standard"); }, "Too few evidence rows must fail");
    auto long_confidence = confidence;
    long_confidence.rows.push_back(evidence());
    rejects([&] { sc::infer(arrays, long_confidence, "standard"); }, "Extra evidence rows must fail");
    auto inconsistent = confidence;
    inconsistent.rows[0].context.observable = false;
    rejects([&] { sc::infer(arrays, inconsistent, "standard"); }, "Contradictory observability must fail");
    inconsistent = confidence;
    inconsistent.rows[0].context_observability = "unexpected";
    rejects([&] { sc::infer(arrays, inconsistent, "standard"); }, "Unknown evidence schema must fail");
    auto invalid = arrays;
    invalid[0].entries[0].second = "ACGT\twrong";
    rejects([&] { sc::infer(invalid, confidence, "standard"); }, "Invalid sequences must fail before TSV writing");
    require(sc::infer(arrays, short_confidence, "off").rows.empty(),
            "Disabled selection does not require confidence evidence");

    Temporary output;
    sc::write_tables(output.path.string(), arrays, selected);
    const auto initial_files = files(output.path);
    require(initial_files.size() == 5, "Five native supported output artifacts expected");
    require(initial_files.at("supported_entries.tsv").find("Array_2\t2\t") != std::string::npos &&
            initial_files.at("supported_entries.tsv").find("Array_2\t1\t") == std::string::npos,
            "Selected entries must retain original row identifiers without licensing other pairs");
    require(initial_files.at("supported_arrays.tsv").find("\tsingleton\tstandard\n") != std::string::npos,
            "Singleton status must be visible in native group output");
    require(initial_files.at("repeat_haplotype_groups.tsv") ==
            "group_id\tparent_array_id\trepeat_haplotype\tentry_index\n",
            "Families without a proper repeat-haplotype subset emit only the explicit header");
    std::ostringstream expected_fasta;
    std::size_t id = 0;
    for (const auto& s : selected.spacers) expected_fasta << ">Zuck_supported_" << ++id << '\n' << s << '\n';
    require(initial_files.at("supported_spacers.fasta") == expected_fasta.str(),
            "FASTA must emit sorted distinct canonical sequences with deterministic identifiers");
    sc::write_tables(output.path.string(), arrays, sc::infer(arrays, confidence, "standard"));
    require(initial_files == files(output.path), "Writing the same view twice must be byte deterministic");

    std::vector<Array> haplotype_arrays{{repeat, {
        {repeat, sequence(26)}, {repeat, sequence(27)},
        {alternative_repeat, sequence(28)}, {alternative_repeat, sequence(29)}}}};
    ec::Result haplotype_confidence;
    haplotype_confidence.rows = {evidence(), evidence(), evidence(), evidence()};
    const auto haplotype_calls = sc::infer(
        haplotype_arrays, haplotype_confidence, "standard");
    require(haplotype_calls.repeat_haplotype_groups.size() == 2 &&
            haplotype_calls.counts.at("supported_repeat_haplotype_memberships") == 4,
            "Exact observed repeat alleles with two spacers must create additive subdivisions");
    Temporary haplotype_output;
    sc::write_tables(haplotype_output.path.string(), haplotype_arrays, haplotype_calls);
    const auto haplotype_table = contents(
        haplotype_output.path / "repeat_haplotype_groups.tsv");
    require(haplotype_table.find("RepeatHaplotype_1\tArray_1\t") != std::string::npos &&
            haplotype_table.find("RepeatHaplotype_2\tArray_1\t") != std::string::npos,
            "Repeat-haplotype memberships must be exported with stable parent and row IDs");
    auto changed_arrays = arrays;
    std::swap(changed_arrays[0].entries[1], changed_arrays[0].entries[2]);
    rejects([&] { sc::write_tables(output.path.string(), changed_arrays, selected); },
            "Writer must reject stale/reordered row identities");
    require(initial_files == files(output.path), "Validation must precede opening output files");
    auto changed_result = selected;
    changed_result.groups[0].entry_indices.pop_back();
    rejects([&] { sc::write_tables(output.path.string(), arrays, changed_result); },
            "Writer must reject inconsistent group membership");

    Temporary empty_output, disabled_output;
    const std::vector<Array> empty_arrays;
    const ec::Result empty_confidence;
    const auto empty = sc::infer(empty_arrays, empty_confidence, "standard");
    sc::write_tables(empty_output.path.string(), empty_arrays, empty);
    require(contents(empty_output.path / "supported_spacers.fasta").empty() &&
            files(empty_output.path).size() == 5,
            "Empty inference writes explicit empty supported inventory and table headers");
    sc::write_tables(disabled_output.path.string(), arrays, sc::infer(arrays, {}, "off"));
    require(files(disabled_output.path).empty(), "Disabled selection must not emit supported outputs");
    std::cout << "Supported profiles, missingness, row identity, canonical inventory, groups and deterministic output passed\n";
}
