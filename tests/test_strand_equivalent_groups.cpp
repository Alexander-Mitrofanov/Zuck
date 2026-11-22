#ifdef NDEBUG
#undef NDEBUG
#endif
#include "crispr_array/strand_equivalent_groups.h"

#include <cassert>
#include <string>
#include <vector>

namespace {
struct Array {
    std::string consensus;
    Zuck::ArrayRows entries;
};

std::string reverse_complement(std::string sequence) {
    for (char& base : sequence) {
        switch (base) {
            case 'A': base = 'T'; break;
            case 'T': base = 'A'; break;
            case 'C': base = 'G'; break;
            case 'G': base = 'C'; break;
            default: break;
        }
    }
    std::reverse(sequence.begin(), sequence.end());
    return sequence;
}

Array reverse_array(const Array& array) {
    Array result{reverse_complement(array.consensus), {}};
    for (const auto& row : array.entries)
        result.entries.emplace_back(reverse_complement(row.first), reverse_complement(row.second));
    return result;
}

Zuck::ArrayRows raw_rows(const std::vector<Array>& arrays) {
    Zuck::ArrayRows rows;
    for (const auto& array : arrays)
        rows.insert(rows.end(), array.entries.begin(), array.entries.end());
    std::sort(rows.begin(), rows.end());
    return rows;
}

std::vector<std::string> literal_pairs(const std::vector<Array>& arrays) {
    std::vector<std::string> pairs;
    for (const auto& array : arrays)
        for (const auto& row : array.entries) pairs.push_back(row.first + row.second);
    std::sort(pairs.begin(), pairs.end());
    return pairs;
}

std::vector<Zuck::ExactArrayKey> group_keys(const std::vector<Array>& arrays) {
    std::vector<Zuck::ExactArrayKey> keys;
    for (const auto& array : arrays)
        keys.push_back(Zuck::exact_array_key(array.consensus, array.entries, reverse_complement));
    std::sort(keys.begin(), keys.end());
    return keys;
}
}  // namespace

int main() {
    const Array original{"AACG", {{"AACG", "CCAT"}, {"AACG", "GGAT"}}};
    Array reordered = original;
    std::reverse(reordered.entries.begin(), reordered.entries.end());
    Array reversed = reverse_array(original);
    std::vector<Array> arrays{original, reordered, reversed};
    const auto initial_rows = raw_rows(arrays);
    const auto initial_pairs = literal_pairs(arrays);
    assert(Zuck::coalesce_strand_equivalent_arrays(arrays, reverse_complement) == 2);
    assert(arrays.size() == 1);
    assert(arrays[0].consensus == original.consensus);
    assert(arrays[0].entries.size() == 6);
    assert(std::equal(original.entries.begin(), original.entries.end(), arrays[0].entries.begin()));
    assert(raw_rows(arrays) == initial_rows);
    assert(literal_pairs(arrays) == initial_pairs);

    // Same consensus/spacers is insufficient if an observed repeat differs.
    Array allele = original;
    allele.entries[1].first = "AATG";
    // Preserve extra evidence rows, even when they duplicate an existing row.
    Array multiplicity = original;
    multiplicity.entries.push_back(original.entries[0]);
    // Each row cannot choose its own strand orientation independently.
    Array mixed_strands = original;
    mixed_strands.entries[1] = reversed.entries[1];
    // Paired repeats must not be permuted among spacers.
    Array paired{"AACG", {{"AACG", "CCAT"}, {"AAGG", "GGAT"}}};
    Array unpaired = paired;
    std::swap(unpaired.entries[0].first, unpaired.entries[1].first);
    arrays = {original, allele, multiplicity, mixed_strands, paired, unpaired};
    assert(Zuck::coalesce_strand_equivalent_arrays(arrays, reverse_complement) == 0);
    assert(arrays.size() == 6);

    // Removing a middle duplicate must preserve the first surviving rows.
    arrays = {original, reversed, allele, reverse_array(allele), multiplicity};
    auto permuted = arrays;
    std::reverse(permuted.begin(), permuted.end());
    const auto before_rows = raw_rows(arrays);
    const auto before_pairs = literal_pairs(arrays);
    assert(Zuck::coalesce_strand_equivalent_arrays(arrays, reverse_complement) == 2);
    assert(arrays.size() == 3);
    assert(arrays[1].entries.size() == 4);
    assert(std::equal(allele.entries.begin(), allele.entries.end(), arrays[1].entries.begin()));
    assert(arrays[2].entries == multiplicity.entries);
    assert(raw_rows(arrays) == before_rows);
    assert(literal_pairs(arrays) == before_pairs);
    assert(Zuck::coalesce_strand_equivalent_arrays(permuted, reverse_complement) == 2);
    assert(group_keys(arrays) == group_keys(permuted));
    assert(raw_rows(permuted) == before_rows);
    assert(literal_pairs(permuted) == before_pairs);
    arrays.clear();
    assert(Zuck::coalesce_strand_equivalent_arrays(arrays, reverse_complement) == 0);
}
