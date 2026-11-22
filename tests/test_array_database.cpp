#include "crispr_array/array_database.h"
#include <chrono>
#include <filesystem>
#include <iostream>

namespace db = Zuck::array_database;
static void require(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}
const std::string repeat = "GTTCACTGCCGTACAGGCAGCTTAGAAA";
const std::string s1 = "ACGTGGCATGACTGATCGTACGATCGTACGAA";
const std::string s2 = "AGGTCTCGATAGCATGCCGATATCCGTAGTCA";
const std::string s3 = "AATGCCGATGTCGTACTAGCTAGTGCACTTGA";
const std::string header = "array_id\taccession\tcategory\tconfidence\trepeat_count\tspacer_count\trepeats\tspacers\n";

int main() {
    const auto root = std::filesystem::temp_directory_path() /
        ("zuck-array-db-" + std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()));
    std::filesystem::create_directory(root);
    const auto path = root / "db.tsv";
    auto row = [&](std::string id, std::string category, std::string score,
                   std::string r, std::string spacers) {
        return id + "\tGCA_123.1\t" + category + "\t" + score + "\t4\t3\t" + r + "\t" + spacers + "\n";
    };
    auto infer = [&](const std::string& rows, const std::vector<std::string>& cores) {
        { std::ofstream out(path); out << header << rows; }
        return db::Database(path.string()).infer(cores);
    };
    const auto all = s1 + ',' + s2 + ',' + s3;
    const auto good = row("a", "Bona-fide", "0.9", repeat, all);
    const std::vector<std::string> cores = {repeat + s1, repeat + s2};
    auto result = infer(good, cores);
    require(result.groups.size() == 1 && result.groups[0].entries.size() == 2,
            "Two same-array, exact observed spacers must support Bona-fide rescue");
    require(result.groups[0].evidence[0].distinct_spacer_matches == 2,
            "Count distinct spacer sequences");
    require(result.groups[0].evidence[0].distinct_repeat_matches == 1,
            "Repeated copies of one repeat are one allele");
    require(std::abs(result.groups[0].evidence[0].score - 1.35) < 1e-12,
            "Confidence and bounded evidence counts must affect score");
    for (const auto& pair : result.groups[0].entries)
        require(pair.second != s3, "A database-only spacer must never be emitted");
    require(infer(good, {repeat + s1, repeat + s1, repeat + s1}).groups.empty(),
            "Traversal multiplicity cannot substitute for spacer diversity");
    require(infer(row("a", "Bona-fide", "0.74", repeat, all), cores).groups.empty(),
            "Low numeric confidence cannot pass the Bona-fide gate");
    require(infer(row("a", "Possible", "0.6", repeat, all), cores).groups.empty(),
            "Possible arrays require stronger same-array diversity");
    require(infer(row("a", "Possible", "0.6", repeat, all),
                  {repeat + s1, repeat + s2, repeat + s3}).groups.size() == 1,
            "Three distinct matches can support Possible rescue");
    require(infer(row("a", "Possible", "0.49", repeat, all),
                  {repeat + s1, repeat + s2, repeat + s3}).groups.empty(),
            "Low numeric confidence cannot pass Possible gate");
    require(infer(row("bad", "Bona-fide", "0.4", repeat, all) +
                  row("valid", "Possible", "0.6", repeat, all),
                  {repeat + s1, repeat + s2, repeat + s3}).groups.size() == 1,
            "Ineligible category records cannot suppress eligible duplicate evidence");
    const auto duplicate = infer(good + row("b", "Bona-fide", "0.9", repeat, all), cores);
    require(duplicate.counts.at("database_distinct_array_signatures") == 1 &&
            duplicate.groups[0].evidence[0].score == result.groups[0].evidence[0].score,
            "Database copies cannot inflate evidence");
    require(infer(row("a", "Bona-fide", "0.9", repeat, s1) +
                  row("b", "Bona-fide", "0.9", repeat, s2), cores).groups.empty(),
            "Global spacer hits from separate arrays do not establish same-array support");
    const auto alternative = row("b", "Bona-fide", "0.9", repeat + "A",
        s1.substr(1) + ',' + s2.substr(1) + ',' + s3.substr(1));
    const auto tied = infer(good + alternative, cores);
    require(tied.groups.empty() && tied.counts.at("database_ambiguous_cores") == 2,
            "Equal evidence for different repeat/spacer boundaries must abstain");
    const auto weaker = row("b", "Bona-fide", "0.8", repeat + "A",
        s1.substr(1) + ',' + s2.substr(1) + ',' + s3.substr(1));
    const auto ranked = infer(good + weaker, cores);
    require(ranked.groups.size() == 1 && ranked.groups[0].consensus == db::canonical(repeat),
            "A stronger complete evidence record can resolve competing boundaries");
    std::vector<std::string> rotated;
    for (const auto& c : cores) rotated.push_back(c.substr(13) + c.substr(0, 13));
    require(infer(good, rotated).groups[0].entries == result.groups[0].entries,
            "Cycle serialization rotation cannot change inferred boundaries");
    std::vector<std::string> reverse;
    for (const auto& c : cores) reverse.push_back(db::reverse_complement(c));
    const auto reversed = infer(good, reverse);
    require(reversed.groups.size() == 1, "Reverse-strand observations must be searchable");
    for (const auto& e : reversed.groups[0].entries)
        require(e.first == db::reverse_complement(repeat), "Preserve observed repeat orientation");
    bool invalid = false;
    try { infer(row("a", "Bona-fide", "nan", repeat, all), cores); }
    catch (const std::runtime_error&) { invalid = true; }
    require(invalid, "Nonfinite database confidence must fail closed");
    for (const std::string count : {"+4", " 4", " -4", "18446744073709551616", "4x"}) {
        std::string invalid_row = good;
        const auto position = invalid_row.find("\t4\t3\t");
        invalid_row.replace(position + 1, 1, count);
        bool rejected = false;
        try { infer(invalid_row, cores); }
        catch (const std::exception&) { rejected = true; }
        require(rejected, "Malformed or overflowing database counts must fail");
    }
    std::filesystem::remove_all(root);
    std::cout << "Database confidence, co-occurrence, duplicate, ambiguity, orientation, and observation controls passed\n";
}
