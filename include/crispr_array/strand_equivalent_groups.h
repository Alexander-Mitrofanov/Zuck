#ifndef ZUCK_STRAND_EQUIVALENT_GROUPS_H
#define ZUCK_STRAND_EQUIVALENT_GROUPS_H

#include <algorithm>
#include <cstddef>
#include <iterator>
#include <map>
#include <string>
#include <utility>
#include <vector>

namespace Zuck {

using ArrayRows = std::vector<std::pair<std::string, std::string>>;
using ExactArrayKey = std::pair<std::string, ArrayRows>;

// Array rows express unordered cycle observations, not genomic spacer order.
// Compare the COMPLETE row multiset in one orientation with the complete
// reverse-complemented orientation. Never canonicalize each component or row
// independently: that could merge distinct repeat/spacer relationships.
template <class ReverseComplement>
ExactArrayKey exact_array_key(const std::string& consensus, const ArrayRows& entries,
                             ReverseComplement reverse_complement) {
    ExactArrayKey forward{consensus, entries};
    ExactArrayKey reverse{reverse_complement(consensus), {}};
    reverse.second.reserve(entries.size());
    for (const auto& entry : entries) {
        reverse.second.emplace_back(reverse_complement(entry.first),
                                    reverse_complement(entry.second));
    }
    std::sort(forward.second.begin(), forward.second.end());
    std::sort(reverse.second.begin(), reverse.second.end());
    return std::min(forward, reverse);
}

// Coalesce equivalent groups while retaining EVERY original oriented row.
// RC(repeat)+RC(spacer) is not the linear reverse complement of repeat+spacer;
// replacing one with the other would discard observed adjacency evidence.
// Keep each original row unchanged and keep the first group's consensus.
// Keys describe complete input groups before any coalescing, never a growing
// merged group's content, so the operation cannot introduce transitive joins.
template <class Array, class ReverseComplement>
std::size_t coalesce_strand_equivalent_arrays(std::vector<Array>& arrays,
                                             ReverseComplement reverse_complement) {
    std::map<ExactArrayKey, std::size_t> observed;
    std::size_t retained = 0;
    const std::size_t original = arrays.size();
    for (std::size_t i = 0; i < original; ++i) {
        auto key = exact_array_key(arrays[i].consensus, arrays[i].entries,
                                   reverse_complement);
        const auto inserted = observed.emplace(std::move(key), retained);
        if (inserted.second) {
            if (retained != i) arrays[retained] = std::move(arrays[i]);
            ++retained;
        } else {
            auto& destination = arrays[inserted.first->second].entries;
            auto& source = arrays[i].entries;
            destination.insert(destination.end(),
                               std::make_move_iterator(source.begin()),
                               std::make_move_iterator(source.end()));
        }
    }
    arrays.resize(retained);
    return original - retained;
}

}  // namespace Zuck

#endif
