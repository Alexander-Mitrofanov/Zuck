#ifndef ZUCK_TERMINAL_VARIANT_RESCUE_H
#define ZUCK_TERMINAL_VARIANT_RESCUE_H

#include <cstddef>
#include <string>
#include <vector>

namespace Zuck {
inline constexpr std::size_t TERMINAL_RESCUE_MIN_OBSERVATIONS = 3;
inline constexpr std::size_t TERMINAL_RESCUE_MAX_COMPONENT_SIZE = 8;
inline constexpr std::size_t TERMINAL_RESCUE_LENGTH_DIFFERENCE = 1;

// Recover a tied, complete observed allele rather than inventing a trimmed
// sequence. Counts are traversal observations, not independent read depth.
inline std::vector<int> co_dominant_terminal_variants(
    const std::vector<std::string>& alleles,
    const std::vector<std::size_t>& support,
    const std::vector<int>& component, int representative) {
    std::vector<int> result;
    if (component.size() > TERMINAL_RESCUE_MAX_COMPONENT_SIZE ||
        support[representative] < TERMINAL_RESCUE_MIN_OBSERVATIONS)
        return result;
    const auto& anchor = alleles[representative];
    for (int candidate : component) {
        const auto& sequence = alleles[candidate];
        if (candidate == representative ||
            support[candidate] != support[representative] ||
            sequence.size() + TERMINAL_RESCUE_LENGTH_DIFFERENCE != anchor.size())
            continue;
        if (anchor.compare(0, sequence.size(), sequence) == 0 ||
            anchor.compare(anchor.size() - sequence.size(), sequence.size(), sequence) == 0)
            result.push_back(candidate);
    }
    return result;
}
}  // namespace Zuck
#endif
