#ifndef ALTERNATIVE_REPEAT_VIEW_H
#define ALTERNATIVE_REPEAT_VIEW_H

#include <algorithm>
#include <cstddef>
#include <string>

namespace alternative_repeat_view {

struct CyclicMatch {
    double identity = 0.0;
    double reciprocal_coverage = 0.0;
};

inline std::string reverse_complement(const std::string& sequence) {
    std::string result(sequence.size(), 'N');
    auto complement = [](char base) {
        switch (base) {
            case 'A': return 'T';
            case 'C': return 'G';
            case 'G': return 'C';
            case 'T': return 'A';
            default: return 'N';
        }
    };
    for (std::size_t i = 0; i < sequence.size(); ++i)
        result[sequence.size() - 1 - i] = complement(sequence[i]);
    return result;
}

inline std::string canonical_sequence(const std::string& sequence) {
    const std::string reverse = reverse_complement(sequence);
    return std::min(sequence, reverse);
}

inline std::string canonical_spacer(const std::string& spacer) {
    return canonical_sequence(spacer);
}

// Scorer-equivalent repeat comparison: align the shorter sequence without
// gaps at every circular offset of the longer sequence, in both orientations.
// Length agreement is deliberately reported separately from identity.
inline CyclicMatch oriented_cyclic_match(const std::string& candidate,
                                         const std::string& parent) {
    if (candidate.empty() || parent.empty()) return {};

    CyclicMatch result;
    result.reciprocal_coverage =
        static_cast<double>(std::min(candidate.size(), parent.size())) /
        static_cast<double>(std::max(candidate.size(), parent.size()));

    for (const std::string& oriented_parent :
         {parent, reverse_complement(parent)}) {
        const std::string& shorter =
            candidate.size() <= oriented_parent.size()
                ? candidate
                : oriented_parent;
        const std::string& longer =
            candidate.size() <= oriented_parent.size()
                ? oriented_parent
                : candidate;
        const std::string doubled = longer + longer;
        for (std::size_t offset = 0; offset < longer.size(); ++offset) {
            std::size_t matches = 0;
            for (std::size_t pos = 0; pos < shorter.size(); ++pos)
                matches += shorter[pos] == doubled[offset + pos];
            result.identity = std::max(
                result.identity,
                static_cast<double>(matches) /
                    static_cast<double>(shorter.size()));
        }
    }
    return result;
}

inline bool eligible_allele(const std::string& candidate_consensus,
                            const std::string& parent_consensus,
                            std::size_t associated_spacers,
                            std::size_t parent_spacers) {
    if (candidate_consensus.empty() ||
        candidate_consensus == parent_consensus ||
        associated_spacers < 4 ||
        parent_spacers == 0 ||
        associated_spacers >= parent_spacers)
        return false;

    const double subset_fraction =
        static_cast<double>(associated_spacers) /
        static_cast<double>(parent_spacers);
    if (subset_fraction < (1.0 / 3.0) ||
        subset_fraction > (2.0 / 3.0))
        return false;

    const CyclicMatch match =
        oriented_cyclic_match(candidate_consensus, parent_consensus);
    return match.identity < 0.80 && match.reciprocal_coverage >= 0.90;
}

}  // namespace alternative_repeat_view

#endif  // ALTERNATIVE_REPEAT_VIEW_H
