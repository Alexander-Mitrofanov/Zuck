#include "crispr_array/array_database.h"
#include <chrono>
#include <filesystem>
#include <iostream>

namespace rd = Zuck::database_reads;
namespace db = Zuck::array_database;
static void require(bool value, const char* message) {
    if (!value) throw std::runtime_error(message);
}
const std::string r = "GTTCACTGCCGTACAGGCAGCTTAGAAA";
const std::string a = "ACGTGGCATGACTGATCGTACGATCGTACGAA";
const std::string b = "AGGTCTCGATAGCATGCCGATATCCGTAGTCA";
const std::string c = "AATGCCGATGTCGTACTAGCTAGTGCACTTGA";

int main() {
    // Overlapping suffix patterns must all be reported, while ambiguous bases
    // interrupt matches. Check coordinates against an independent literal scan.
    Zuck::read_evidence::ExactMatcher matcher;
    const std::vector<std::string> patterns{"AC", "CAC", "ACAC", "C", "GAC"};
    for (std::size_t id = 0; id < patterns.size(); ++id) matcher.add(patterns[id], id);
    matcher.build();
    const std::string sequence = "ACACNACGACCAC";
    std::set<std::pair<std::size_t, std::size_t>> observed, expected;
    matcher.scan(sequence, [&](std::size_t id, std::size_t end) { observed.emplace(id, end); });
    for (std::size_t id = 0; id < patterns.size(); ++id)
        for (std::size_t start = 0; start + patterns[id].size() <= sequence.size(); ++start)
            if (sequence.compare(start, patterns[id].size(), patterns[id]) == 0)
                expected.emplace(id, start + patterns[id].size());
    require(observed == expected, "Shared matcher must preserve overlapping suffix matches and coordinates");
    // Large low-quality prefixes must not change the decision for a later
    // high-quality span through floating-point cancellation in prefix sums.
    const std::string long_quality = std::string(2000000, '!') + std::string(10, '+');
    Zuck::read_evidence::QualitySpans long_spans(long_quality, 10);
    double direct_errors = 0;
    for (int i = 0; i < 10; ++i) direct_errors += std::pow(10.0, -10.0 / 10.0);
    require(long_spans.passes(2000000, 2000010, 1.0) == (direct_errors <= 1.0),
            "Quality threshold must agree with direct summation after long error-heavy prefixes");
    const auto root = std::filesystem::temp_directory_path() /
        ("zuck-db-reads-" + std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()));
    std::filesystem::create_directory(root);
    const auto path = root / "reads.fastq";
    const auto span = r + a + r.substr(0, 16);
    auto write = [&](const std::vector<std::string>& reads, char quality = 'I') {
        std::ofstream out(path);
        for (std::size_t i = 0; i < reads.size(); ++i)
            out << '@' << i << '\n' << reads[i] << "\n+\n" << std::string(reads[i].size(), quality) << '\n';
    };
    const std::vector<std::pair<std::string, std::string>> pairs = {
        {r, a}, {db::reverse_complement(r), db::reverse_complement(a)}};
    write({span + "AA", "CCC" + span});
    auto result = rd::validate(pairs, {path.string()});
    require(result.support[0].accepted(), "Two exact high-quality spans at distinct starts must pass");
    require(!result.support[1].accepted(), "Do not reverse-complement repeat and spacer independently");
    write({db::reverse_complement(span) + "AA", "CCC" + db::reverse_complement(span)});
    require(rd::validate(pairs, {path.string()}).support[0].accepted(), "Whole-span reverse complement must pass");
    write({span, span, span});
    require(!rd::validate(pairs, {path.string()}).support[0].accepted(), "Duplicated reads cannot supply two units");
    write({span + "AA", span + "CC"});
    require(!rd::validate(pairs, {path.string()}).support[0].accepted(), "Distinct units need distinct read start offsets");
    write({span + span, span + span});
    require(!rd::validate(pairs, {path.string()}).support[0].accepted(), "Repeated pattern positions in one sequence cannot supply two units");
    write({span + "AA", "CCC" + span}, '!');
    require(!rd::validate(pairs, {path.string()}).support[0].accepted(), "Low-quality spans cannot pass");
    write({span + "AA", "CCC" + span}, '+');
    require(!rd::validate(pairs, {path.string()}).support[0].accepted(), "Q10 alone must not evade expected-error bound");
    write({r + a + "AA", "CCC" + r + a});
    require(!rd::validate(pairs, {path.string()}).support[0].accepted(), "A literal pair without next-repeat context is insufficient");
    { std::ofstream out(path); out << ">a\n" << span << "\n>b\nCCC" << span << '\n'; }
    require(!rd::validate(pairs, {path.string()}).support[0].accepted(), "FASTA provides no quality support");

    const auto catalog = root / "db.tsv";
    { std::ofstream out(catalog); out << "array_id\taccession\tcategory\tconfidence\trepeat_count\tspacer_count\trepeats\tspacers\n"
        << "source\tGCA_1.1\tBona-fide\t0.9\t3\t2\t" << r << '\t' << a << ',' << b << '\n'; }
    db::Database database(catalog.string());
    const auto bspan = r + b + r.substr(0, 16), cspan = r + c + r.substr(0, 16);
    write({span + "AA", "CCC" + span});
    require(database.infer({r+a, r+b}, {path.string()}, true).groups.empty(),
            "One read-validated known spacer cannot borrow an unsupported spacer vote");
    write({span + "AA", "CCC" + span, bspan + "AA", "CCC" + bspan, cspan + "AA", "CCC" + cspan});
    const auto known = database.infer({r+a, r+b, r+c}, {path.string()}, true);
    require(known.groups.size() == 1 && known.groups[0].entries.size() == 2,
            "Known-only default must not emit a novel spacer");
    const auto novel = database.infer({r+a, r+b, r+c}, {path.string()}, true, true);
    require(novel.groups.size() == 1 && novel.groups[0].entries.size() == 3,
            "Explicit novel arm needs known family anchors and its own read span");
    bool found = false;
    for (const auto& e : novel.groups[0].evidence) if (e.spacer == c) {
        found = true;
        require(e.novel_spacer && e.read_validation_applied && e.distinct_spacer_matches == 2,
                "Novel spacer must be labelled and cannot add its own known-match vote");
    }
    require(found, "Novel fixture must expose its distinct observed sequence");
    bool invalid = false;
    try { database.infer({r+a, r+b, r+c}, {}, false, true); }
    catch (const std::runtime_error&) { invalid = true; }
    require(invalid, "Novel candidates require read support");
    std::filesystem::remove_all(root);
    std::cout << "Database linear-span, quality, duplication, strand, same-array, and novel controls passed\n";
}
