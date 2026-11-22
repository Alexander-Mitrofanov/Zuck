#include "crispr_array/evaluation.h"
#include "crispr_array/post_processing_helpers.h"

#include <cmath>
#include <iostream>
#include <random>
#include <set>
#include <stdexcept>

namespace {
void require(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}

std::size_t reference_distance(const std::string& a, const std::string& b) {
    std::vector<std::vector<std::size_t>> matrix(a.size() + 1,
        std::vector<std::size_t>(b.size() + 1));
    for (std::size_t i = 0; i <= a.size(); ++i) matrix[i][0] = i;
    for (std::size_t j = 0; j <= b.size(); ++j) matrix[0][j] = j;
    for (std::size_t i = 1; i <= a.size(); ++i)
        for (std::size_t j = 1; j <= b.size(); ++j)
            matrix[i][j] = std::min({matrix[i - 1][j] + 1, matrix[i][j - 1] + 1,
                matrix[i - 1][j - 1] + (a[i - 1] != b[j - 1])});
    return matrix.back().back();
}

std::set<std::set<int>> reference_components(const std::vector<std::string>& alleles, std::size_t k) {
    std::vector<std::set<std::string>> kmers(alleles.size());
    for (std::size_t i = 0; i < alleles.size(); ++i)
        for (std::size_t pos = 0; pos + k <= alleles[i].size(); ++pos)
            kmers[i].insert(alleles[i].substr(pos, k));
    std::set<std::set<int>> result;
    std::set<int> visited;
    for (int start = 0; start < static_cast<int>(alleles.size()); ++start) {
        if (!visited.insert(start).second) continue;
        std::vector<int> queue{start};
        std::set<int> component{start};
        for (std::size_t head = 0; head < queue.size(); ++head) {
            for (int next = 0; next < static_cast<int>(alleles.size()); ++next) {
                if (visited.count(next)) continue;
                for (const auto& kmer : kmers[queue[head]]) {
                    if (!kmers[next].count(kmer)) continue;
                    visited.insert(next);
                    component.insert(next);
                    queue.push_back(next);
                    break;
                }
            }
        }
        result.insert(std::move(component));
    }
    return result;
}
}

int main() {
    using namespace Zuck::post_processing;
    require(bounded_edit_distance("A", "C", std::numeric_limits<int>::max() - 1) == 1,
            "A large edit bound must not overflow its sentinel arithmetic");
    require(get_string_similarity("", "") == 1.0f, "Empty sequences are identical");
    require(get_string_similarity("", "ACGT") == 0.0f, "One empty sequence has zero similarity");
    require(get_levenshtein_distance("", std::string(70000, 'A')) == 70000,
            "Edit distance must not truncate at 16 bits");
    require(get_levenshtein_distance("A", std::string(70000, 'A')) == 69999,
            "Dynamic programming cells must not truncate at 16 bits");
    require(get_number_of_duplicate_spacers({"", "AA"}, "AAAA") == 2,
            "Ignore empty spacers and retain overlapping occurrences");
    const std::vector<std::string> choices{"AC", "AG"};
    require(get_most_similar_sequence("AT", choices) == "AC", "Similarity ties use first input");
    require(choices.size() == 2, "Similarity lookup must not remove choices");
    require(glue_kmers("acGT CgTa GTAC") == "ACGTAC", "Normalize case before gluing overlaps");
    require(glue_kmers("ACGTAC") == "ACGTAC", "Preserve compact input");
    require(glue_kmers("   \t\r") == "", "Whitespace has no sequence");
    require(glue_kmers("AAAA AAAA") == "AAAA", "Preserve full-overlap legacy behavior");
    bool invalid = false;
    try { glue_kmers("ACGT >bad"); } catch (const std::invalid_argument&) { invalid = true; }
    require(invalid, "Reject invalid cycle characters");

    std::mt19937 random(409712);
    const auto dna = [&](std::size_t length) {
        std::string result(length, 'A');
        for (char& base : result) base = "ACGT"[random() % 4];
        return result;
    };
    for (int trial = 0; trial < 5000; ++trial) {
        const auto a = dna(random() % 30), b = dna(random() % 30);
        const auto expected = reference_distance(a, b);
        require(get_levenshtein_distance(a, b) == expected, "Linear-memory distance differs from oracle");
        for (int limit = 0; limit <= 8; ++limit)
            require(bounded_edit_distance(a, b, limit) ==
                static_cast<int>(std::min(expected, static_cast<std::size_t>(limit + 1))),
                "Banded distance differs from exhaustive oracle");
    }
    for (int trial = 0; trial < 1000; ++trial) {
        std::vector<std::string> alleles;
        for (unsigned i = 0, count = random() % 20; i < count; ++i)
            alleles.push_back(dna(random() % 15));
        const std::size_t k = 1 + random() % 6;
        std::set<std::set<int>> actual;
        for (const auto& component : shared_kmer_components(alleles, k))
            actual.emplace(component.begin(), component.end());
        require(actual == reference_components(alleles, k), "K-mer components differ from pairwise graph oracle");
    }
    std::cout << "Sequence evaluation, cycle parsing, banded distance and component regressions passed\n";
}
