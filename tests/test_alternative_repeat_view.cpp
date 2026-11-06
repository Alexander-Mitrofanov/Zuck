#include <cmath>
#include <iostream>
#include <set>
#include <stdexcept>
#include <string>

#include "crispr_array/alternative_repeat_view.h"

int main() {
    using alternative_repeat_view::canonical_spacer;
    using alternative_repeat_view::canonical_sequence;
    using alternative_repeat_view::eligible_allele;
    using alternative_repeat_view::oriented_cyclic_match;
    using alternative_repeat_view::reverse_complement;

    auto require = [](bool condition, const char* message) {
        if (!condition) throw std::runtime_error(message);
    };

    const std::string allele = "GTTTCTAATGTACCGTGTGGAGTTAAAAC";
    const std::string parent = "CGTTTCTAATGTACCGTAGTGGAGTTGAAAC";
    const auto accepted = oriented_cyclic_match(allele, parent);
    require(std::abs(accepted.identity - 0.6896551724137931) < 1e-12,
            "accepted identity changed");
    require(std::abs(accepted.reciprocal_coverage -
                     0.9354838709677419) < 1e-12,
            "accepted coverage changed");
    require(eligible_allele(allele, parent, 4, 9),
            "accepted fixture was rejected");

    // The subset endpoints are inclusive, while three observations and a
    // strict majority above two thirds are rejected.
    require(eligible_allele(allele, parent, 4, 12),
            "one-third boundary was rejected");
    require(eligible_allele(allele, parent, 6, 9),
            "two-thirds boundary was rejected");
    require(!eligible_allele(allele, parent, 3, 9),
            "minimum-four gate failed");
    require(!eligible_allele(allele, parent, 7, 9),
            "two-thirds upper gate failed");

    // A phase shift and reverse complement remain the same repeat family and
    // therefore cannot create a low-identity alternative view.
    const std::string phase_shifted = parent.substr(7) + parent.substr(0, 7);
    require(oriented_cyclic_match(phase_shifted, parent).identity == 1.0,
            "phase-shift identity changed");
    require(oriented_cyclic_match(reverse_complement(phase_shifted), parent)
                .identity == 1.0,
            "reverse-complement identity changed");
    require(!eligible_allele(phase_shifted, parent, 4, 9),
            "phase shift created an alternative");
    require(!eligible_allele(allele.substr(4), parent, 4, 9),
            "low-coverage allele was accepted");

    // Exact and reverse-complement spacers share one deterministic key.
    const std::string spacer = "AACCGGTTACGATCGTACGTA";
    const std::string spacer_rc = reverse_complement(spacer);
    require(canonical_spacer(spacer) == canonical_spacer(spacer_rc),
            "spacer canonicalization is orientation-dependent");
    require(canonical_sequence(allele) ==
                canonical_sequence(reverse_complement(allele)),
            "repeat canonicalization is orientation-dependent");
    std::set<std::string> unique_spacers = {
        canonical_spacer(spacer), canonical_spacer(spacer_rc)};
    require(unique_spacers.size() == 1,
            "exact/RC spacers were not deduplicated");

    std::cout << "alternative repeat-view tests passed\n";
    return 0;
}
