#ifndef ZUCK_READ_RECRUITMENT_H
#define ZUCK_READ_RECRUITMENT_H

#include "crispr_array/read_linked_groups.h"
#include "xxhash/xxhash.h"
#include <cstdint>
#include <iterator>
#include <numeric>
#include <tuple>

namespace Zuck { namespace read_recruitment {

constexpr std::size_t SEED_LENGTH = 11;
constexpr std::size_t MINIMUM_REPEAT = 20, MAXIMUM_REPEAT = 55;
constexpr std::size_t MINIMUM_SPACER = 26, MAXIMUM_SPACER = 50;
constexpr std::size_t MINIMUM_PARTIAL_REPEAT = 16;
constexpr std::size_t MAXIMUM_REPEAT_MISMATCHES = 1;
constexpr int MINIMUM_QUALITY = 10;
constexpr double MAXIMUM_EXPECTED_ERRORS = 1.0;
constexpr std::size_t MINIMUM_SEQUENCE_UNITS = 2, MINIMUM_START_OFFSETS = 2;
constexpr std::size_t MINIMUM_FAMILY_SPACERS = 2, SUPPORT_CAP = 32;
constexpr double ERROR_SCREEN_ALPHA = 0.05;
// Competitive sequence explanation, not a biological tandem probability.
// Keep the strict CRASS similarity threshold as the exact rational 41/50.
constexpr double TANDEM_MINIMUM_SIMILARITY = 0.82;
constexpr std::size_t TANDEM_MINIMUM_TRAINING_COPIES = 2;

struct Config {
    std::size_t maximum_candidates_per_family = 4096;
    std::size_t maximum_total_candidates = 65536;
    std::size_t maximum_family_memberships = 262144;
    std::size_t maximum_seed_placements_per_read = 4096;
    std::size_t maximum_candidates_per_read = 512;
    std::size_t maximum_repeat_verifications_per_read = 16384;
    std::size_t maximum_read_length = 16384;
    std::size_t maximum_opportunity_matches_per_read = 262144;
};
struct Support {
    std::size_t sequence_units = 0, start_offsets = 0;
    bool units_saturated = false, offsets_saturated = false;
    bool accepted() const {
        return sequence_units >= MINIMUM_SEQUENCE_UNITS && start_offsets >= MINIMUM_START_OFFSETS;
    }
};
struct TandemEvidence {
    bool explained = false;
    std::string model = "none", training_kind = "none", motif;
    std::string training_side = "none", model_orientation = "none";
    std::size_t period = 0, spacer_template_length = 0;
    std::size_t left_observed_bases = 0, right_observed_bases = 0;
    double repeat_similarity = 0, left_similarity = 0, spacer_similarity = 0;
    double right_similarity = 0, whole_similarity = 0;
    // CRASS-like repeat diagnostics do not independently veto a candidate.
    double repeat_spacer_similarity = 0, repeat_max_base_fraction = 0;
    double repeat_max_3mer_fraction = 0;
};
struct Candidate {
    std::set<std::size_t> family_owners;
    std::string left_repeat, spacer, right_repeat;
    Support support;
    std::string evidence_kind, reason, core;
    bool accepted = false, graph_proposed = false;
    bool left_complete = true, right_complete = true;
    std::size_t left_observed_bases = 0, right_observed_bases = 0;
    std::size_t distinct_supported_family_spacers = 0;
    bool boundary_conflict = false, resource_limited = false;
    bool pair_only = false;
    std::vector<std::string> observed_right_repeats;
    std::vector<std::size_t> source_candidate_indices; // Zero-based raw ledger rows.
    // Models belong to raw observed triplets; directed aggregates carry the
    // union of sequence units/offsets from sources without either explanation.
    TandemEvidence tandem;
    Support non_tandem_support;
    bool tandem_ambiguous = false;
    std::string tandem_status = "not_admissible";
    // Observation and automatic selection are separate: an ambiguous exact
    // allele remains admissible/inclusive, with its original evidence intact.
    bool selection_eligible = false, error_ambiguous = false;
    std::string error_status = "not_admissible";
    std::uint64_t opportunity_records = 0, error_hypotheses = 0;
    int opportunity_min_spacer_quality = -1;
    bool opportunity_complete = false;
    std::size_t error_neighbor_count = 0;
    std::size_t error_neighbor_index = std::numeric_limits<std::size_t>::max();
    std::size_t error_observed_units = 0;
    double error_probability_upper = 1, error_log_bound = 0, error_log_threshold = 0;
};
struct Result {
    std::vector<Candidate> candidates;
    std::map<std::string, std::size_t> counts;
};

using Pair = std::pair<std::string, std::string>;

// Only an actually known full upstream repeat can become an emitted pair.
// Reverse orientation swaps flanking repeats; component-wise RC is incorrect
// when the two observed repeat copies differ.
inline std::vector<Pair> oriented_pairs(const Candidate& candidate) {
    std::set<Pair> pairs;
    if (candidate.left_complete) pairs.emplace(candidate.left_repeat, candidate.spacer);
    if (candidate.right_complete && !candidate.pair_only)
        pairs.emplace(read_linked::reverse_complement(candidate.right_repeat),
                      read_linked::reverse_complement(candidate.spacer));
    return {pairs.begin(), pairs.end()};
}

namespace detail {
inline std::string rc(const std::string& sequence) { return read_linked::reverse_complement(sequence); }
inline bool dna(const std::string& sequence) {
    return !sequence.empty() && std::all_of(sequence.begin(), sequence.end(),
        [](char c) { return read_evidence::base(c) >= 0; });
}
inline std::string circular_key(const std::string& sequence) {
    if (sequence.empty()) return {};
    std::string best = sequence;
    for (std::size_t i = 1; i < sequence.size(); ++i)
        best = std::min(best, sequence.substr(i) + sequence.substr(0, i));
    return best;
}
inline std::string core(const Pair& pair) {
    const auto sequence = pair.first + pair.second;
    return std::min(circular_key(sequence), circular_key(rc(sequence)));
}
inline std::string canonical(const std::string& sequence) { return std::min(sequence, rc(sequence)); }

// Read strings may contain ambiguity symbols outside the verified span. Unit
// identity still includes them; this complement is an involution on IUPAC DNA.
inline std::string unit_rc(const std::string& sequence) {
    std::string result(sequence.rbegin(), sequence.rend());
    const std::string alphabet = "ACGTRYMKBDHVSWN", complement = "TGCAYRKMVHDBSWN";
    for (auto& base : result) {
        const auto position = alphabet.find(base);
        if (position != std::string::npos) base = complement[position];
    }
    return result;
}
using Unit = std::array<std::uint64_t, 2>;
inline Unit unit_identity(const std::string& first, const std::string& second, bool paired) {
    const auto forward = std::make_pair(first, paired ? second : "");
    const auto reverse = paired ? std::make_pair(unit_rc(second), unit_rc(first)) :
                                 std::make_pair(unit_rc(first), std::string{});
    const auto& ordered = std::min(forward, reverse);
    const auto sequence = ordered.first + '\0' + ordered.second;
    // Two seeded 64-bit fingerprints bound retained evidence memory. They are
    // not molecule IDs or a cryptographic identity guarantee. A collision can
    // only merge distinct sequence units and therefore under-count support.
    return {{XXH64(sequence.data(), sequence.size(), 0),
             XXH64(sequence.data(), sequence.size(), UINT64_C(0x9e3779b97f4a7c15))}};
}
inline std::size_t edit_distance(const std::string& a, const std::string& b) {
    std::vector<std::size_t> row(b.size() + 1);
    std::iota(row.begin(), row.end(), 0);
    for (std::size_t i = 0; i < a.size(); ++i) {
        auto diagonal = row[0]; row[0] = i + 1;
        for (std::size_t j = 0; j < b.size(); ++j) {
            const auto above = row[j + 1];
            row[j + 1] = std::min({row[j] + 1, above + 1, diagonal + (a[i] != b[j])});
            diagonal = above;
        }
    }
    return row.back();
}
inline bool verifies(const std::string& sequence, std::size_t start, std::size_t length,
                     const std::string& repeat) {
    if (length > repeat.size() || start > sequence.size() || length > sequence.size() - start) return false;
    std::size_t differences = 0;
    for (std::size_t i = 0; i < length; ++i) {
        if (read_evidence::base(sequence[start + i]) < 0) return false;
        differences += sequence[start + i] != repeat[i];
        if (differences > MAXIMUM_REPEAT_MISMATCHES) return false;
    }
    return true;
}
// A strand is relative to one stable family anchor, not an independently
// canonicalized allele. Close repeat alleles can canonicalize on opposite strands.
using Owner = std::pair<std::size_t, bool>;
using Owners = std::set<Owner>;
inline Owners common(const Owners& a, const Owners& b) {
    Owners result;
    std::set_intersection(a.begin(), a.end(), b.begin(), b.end(), std::inserter(result, result.end()));
    return result;
}
inline std::set<std::size_t> families(const Owners& owners) {
    std::set<std::size_t> result;
    for (const auto& owner : owners) result.insert(owner.first);
    return result;
}
using Key = std::tuple<std::string, std::string, std::string, std::string, std::size_t, std::size_t>;
inline Key key(const Candidate& candidate) {
    return {candidate.evidence_kind, candidate.left_repeat, candidate.spacer, candidate.right_repeat,
            candidate.left_observed_bases, candidate.right_observed_bases};
}
struct Pending {
    Candidate candidate;
    std::set<std::size_t> offsets;
    std::map<Pair, std::set<std::size_t>> pair_offsets;
};
struct Accumulator {
    Candidate candidate;
    std::set<Unit> units;
    std::set<std::size_t> offsets;
    std::map<Pair, std::set<std::size_t>> pair_offsets;
    std::set<Unit> non_tandem_units;
    std::set<std::size_t> non_tandem_offsets;
};
struct Seed { std::string sequence; Owners owners; };
struct Anchor { std::size_t seed, offset, length; };
struct Hit { std::size_t start, end; std::string sequence; Owners owners; };
struct GraphPattern { Candidate candidate; std::string repeat, spacer; bool reversed = false; };
} // namespace detail

// Fit one training direction. The physical-triplet wrapper below decides which
// full repeat copies are actually available to train a competing explanation.
// The exact spacer is held out from motif learning.
inline TandemEvidence tandem_evidence_from_left(const Candidate& candidate) {
    using namespace detail;
    if (candidate.pair_only || candidate.left_repeat.size() < MINIMUM_REPEAT ||
        candidate.left_repeat.size() > MAXIMUM_REPEAT || candidate.spacer.size() < MINIMUM_SPACER ||
        candidate.spacer.size() > MAXIMUM_SPACER || candidate.right_repeat.size() > MAXIMUM_REPEAT ||
        candidate.left_observed_bases < MINIMUM_PARTIAL_REPEAT ||
        candidate.right_observed_bases < MINIMUM_PARTIAL_REPEAT ||
        candidate.left_observed_bases > candidate.left_repeat.size() ||
        candidate.right_observed_bases > candidate.right_repeat.size() ||
        !dna(candidate.left_repeat) || !dna(candidate.spacer) || !dna(candidate.right_repeat))
        throw std::invalid_argument("Invalid raw tandem evidence geometry");
    const auto& repeat = candidate.left_repeat;
    const auto& spacer = candidate.spacer;
    const auto left = repeat.substr(repeat.size() - candidate.left_observed_bases);
    const auto right = candidate.right_repeat.substr(0, candidate.right_observed_bases);
    TandemEvidence result;
    result.training_kind = candidate.left_observed_bases == repeat.size() ?
        "observed_full_repeat" : "graph_proposed_repeat";
    result.left_observed_bases = left.size();
    result.right_observed_bases = right.size();
    auto similarity = [](std::size_t differences, std::size_t length) {
        return 1.0 - static_cast<double>(differences) / static_cast<double>(length);
    };
    auto fits = [](std::size_t differences, std::size_t length) {
        return 50 * differences < 9 * length; // Strictly greater than .82.
    };
    result.repeat_spacer_similarity = similarity(edit_distance(repeat, spacer),
        std::max(repeat.size(), spacer.size()));
    std::array<std::size_t, 4> bases{};
    std::map<std::string, std::size_t> trimers;
    for (char base : repeat) ++bases[read_evidence::base(base)];
    for (std::size_t i = 0; i + 3 <= repeat.size(); ++i) ++trimers[repeat.substr(i, 3)];
    result.repeat_max_base_fraction = static_cast<double>(*std::max_element(bases.begin(), bases.end())) / repeat.size();
    std::size_t maximum_trimer = 0;
    for (const auto& item : trimers) maximum_trimer = std::max(maximum_trimer, item.second);
    // Unlike the historical CRASS loop, count every valid 3-mer window; this is
    // an explicitly normalized diagnostic, not a byte-for-byte CRASS flag.
    result.repeat_max_3mer_fraction = static_cast<double>(maximum_trimer) / (repeat.size() - 2);
    auto mismatches = [](const std::string& sequence, const std::string& motif, std::size_t phase) {
        std::size_t differences = 0;
        for (std::size_t i = 0; i < sequence.size(); ++i)
            differences += sequence[i] != motif[(phase + i) % motif.size()];
        return differences;
    };
    auto assess = [&](const std::string& motif, std::size_t template_length, std::size_t spacer_edits,
                      std::size_t repeat_edits, const std::string& model) {
        const auto left_edits = mismatches(left, motif, repeat.size() - left.size());
        const auto right_edits = mismatches(right, motif, repeat.size() + template_length);
        const auto spacer_span = std::max(spacer.size(), template_length);
        const auto whole_span = left.size() + spacer_span + right.size();
        const auto whole_edits = left_edits + spacer_edits + right_edits;
        if (!fits(left_edits, left.size()) || !fits(spacer_edits, spacer_span) ||
            !fits(right_edits, right.size()) || !fits(whole_edits, whole_span)) return false;
        result.explained = true;
        result.model = model;
        result.motif = motif;
        result.period = motif.size();
        result.spacer_template_length = template_length;
        result.repeat_similarity = similarity(repeat_edits, repeat.size());
        result.left_similarity = similarity(left_edits, left.size());
        result.spacer_similarity = similarity(spacer_edits, spacer_span);
        result.right_similarity = similarity(right_edits, right.size());
        result.whole_similarity = similarity(whole_edits, whole_span);
        return true;
    };
    // Require two periods in the full training repeat and spacer. Clipped
    // flanks validate continuation without changing the model universe: a
    // period-nine full witness must not become "non-tandem" merely because
    // the same path is also seen with a 16-base terminal repeat fragment.
    const auto maximum_period = std::min(repeat.size(), spacer.size()) /
        TANDEM_MINIMUM_TRAINING_COPIES;
    // Any passing spacer edit alignment has template length < |S|/.82.
    // One DP row supplies distances to every permitted periodic prefix.
    const auto maximum_template = (50 * spacer.size() + 40) / 41;
    for (std::size_t period = 1; period <= maximum_period; ++period) {
        std::vector<std::array<std::size_t, 4>> counts(period);
        for (std::size_t i = 0; i < repeat.size(); ++i) ++counts[i % period][read_evidence::base(repeat[i])];
        std::string motif;
        for (const auto& count : counts)
            motif += "ACGT"[std::distance(count.begin(), std::max_element(count.begin(), count.end()))];
        const auto repeat_edits = mismatches(repeat, motif, 0);
        if (!fits(repeat_edits, repeat.size())) continue;
        std::vector<std::size_t> distances(maximum_template + 1);
        std::iota(distances.begin(), distances.end(), 0);
        for (std::size_t i = 0; i < spacer.size(); ++i) {
            auto diagonal = distances[0]; distances[0] = i + 1;
            for (std::size_t j = 0; j < maximum_template; ++j) {
                const auto above = distances[j + 1];
                distances[j + 1] = std::min({distances[j] + 1, above + 1,
                    diagonal + (spacer[i] != motif[(repeat.size() + j) % period])});
                diagonal = above;
            }
        }
        for (std::size_t length = 1; length <= maximum_template; ++length)
            if (assess(motif, length, distances[length], repeat_edits, "short_period")) return result;
    }
    // A full repeat-sized spacer has another simple three-copy explanation,
    // even when the repeat has no shorter internal period. This transfers the
    // CRASS repeat/spacer similarity test with explicit flanking geometry.
    assess(repeat, repeat.size(), edit_distance(spacer, repeat), 0, "repeat_copy");
    return result;
}

// A physical witness is tandem-compatible if either actually full observed
// repeat trains a motif explaining BOTH observed repeat segments and all of S.
// Evaluate canonical raw left, then right (whole-triplet RC with swapped flanks).
// One valid competing explanation belongs to the whole physical witness and is
// intentionally shared by both emitted directions. This does not claim that
// both direction-specific motif learners pass. For graph junctions with neither
// copy fully observed, the known graph left repeat supplies an explicit proposal.
inline TandemEvidence tandem_evidence(const Candidate& candidate) {
    using namespace detail;
    const bool left_full = candidate.left_complete && candidate.left_observed_bases == candidate.left_repeat.size();
    const bool right_full = candidate.right_complete && candidate.right_observed_bases == candidate.right_repeat.size();
    auto forward = tandem_evidence_from_left(candidate);
    const bool graph_proposal = !left_full && !right_full && candidate.evidence_kind == "graph_junction";
    if ((left_full || graph_proposal) && forward.explained) {
        forward.training_side = "left"; forward.model_orientation = "forward";
        return forward;
    }
    if (right_full) {
        Candidate reversed = candidate;
        reversed.left_repeat = rc(candidate.right_repeat);
        reversed.spacer = rc(candidate.spacer);
        reversed.right_repeat = rc(candidate.left_repeat);
        std::swap(reversed.left_complete, reversed.right_complete);
        std::swap(reversed.left_observed_bases, reversed.right_observed_bases);
        auto reverse = tandem_evidence_from_left(reversed);
        if (reverse.explained) {
            reverse.training_side = "right"; reverse.model_orientation = "reverse";
            // Repeat-only diagnostics always describe the canonical raw LEFT
            // repeat, whereas model fits/spans describe the declared fit frame.
            reverse.repeat_spacer_similarity = forward.repeat_spacer_similarity;
            reverse.repeat_max_base_fraction = forward.repeat_max_base_fraction;
            reverse.repeat_max_3mer_fraction = forward.repeat_max_3mer_fraction;
            return reverse;
        }
    }
    TandemEvidence unexplained;
    unexplained.left_observed_bases = candidate.left_observed_bases;
    unexplained.right_observed_bases = candidate.right_observed_bases;
    unexplained.repeat_spacer_similarity = forward.repeat_spacer_similarity;
    unexplained.repeat_max_base_fraction = forward.repeat_max_base_fraction;
    unexplained.repeat_max_3mer_fraction = forward.repeat_max_3mer_fraction;
    return unexplained;
}

// A union bound on at least K independent error successes among N opportunities.
// K is a lower bound when unit support saturates. No subtraction from a nearly
// unit CDF, exponentiation of tiny probabilities, or iterative allele removal.
inline double independent_error_log_bound(std::uint64_t opportunities, std::size_t units, double probability) {
    if (!units || opportunities < units || probability <= 0 || probability > 1 || !std::isfinite(probability))
        throw std::invalid_argument("Invalid independent-error bound inputs");
    long double value = static_cast<long double>(units) * std::log(static_cast<long double>(probability));
    for (std::size_t i = 0; i < units; ++i)
        value += std::log(static_cast<long double>(opportunities - i)) -
                 std::log(static_cast<long double>(i + 1));
    return static_cast<double>(std::min(static_cast<long double>(0), value));
}

// Conditional independent, calibrated-substitution error screen, NOT a
// calibrated biological-allele probability or unconditional familywise bound.
// The observed-parent hypothesis universe is fixed before applying any veto.
// Broad suffix-anchored opportunity windows deliberately include other spacer
// interiors and arbitrary following bases; they can make this screen loose.
inline void annotate_error_uncertainty(Result& result, const std::vector<std::string>& files,
                                       const Config& config) {
    using namespace detail;
    using BucketKey = std::pair<std::string, std::size_t>;
    struct Opportunity {
        std::string repeat;
        std::size_t length = 0;
        std::uint64_t records = 0;
        int minimum_quality = 94;
        bool complete = true;
    };
    std::map<BucketKey, std::size_t> bucket_ids;
    std::vector<Opportunity> opportunities;
    std::vector<std::size_t> eligible;
    std::uint64_t hypotheses = 0;
    bool complete = true;
    for (std::size_t i = 0; i < result.candidates.size(); ++i) {
        auto& candidate = result.candidates[i];
        candidate.selection_eligible = candidate.accepted && !candidate.tandem_ambiguous;
        candidate.error_status = candidate.pair_only ? "not_admissible" : "provenance_only";
        if (!candidate.accepted) continue;
        if (!candidate.pair_only) throw std::logic_error("Only directed pair evidence can pass the error screen");
        eligible.push_back(i);
        const auto length = candidate.spacer.size();
        if (hypotheses > std::numeric_limits<std::uint64_t>::max() - 3 * length) complete = false;
        else hypotheses += 3 * length;
        const auto bucket = bucket_ids.emplace(BucketKey{candidate.left_repeat, length}, opportunities.size());
        if (bucket.second) opportunities.push_back({candidate.left_repeat, length, 0, 94, true});
    }
    // One wildcard position at a time bounds indexing memory by the number of
    // exact pairs, rather than pairs times spacer length. Each DNA bucket has
    // at most four distinct alleles. R orientation is part of the exact key.
    for (std::size_t position = 0; position < MAXIMUM_SPACER; ++position) {
        std::map<std::pair<std::string, std::string>, std::vector<std::size_t>> neighbors;
        for (auto i : eligible) {
            const auto& candidate = result.candidates[i];
            if (position >= candidate.spacer.size()) continue;
            auto masked = candidate.spacer;
            masked[position] = '*';
            neighbors[{candidate.left_repeat, std::move(masked)}].push_back(i);
        }
        for (const auto& entry : neighbors) for (auto i : entry.second) for (auto j : entry.second) {
            if (i == j || result.candidates[i].spacer == result.candidates[j].spacer) continue;
            auto& candidate = result.candidates[i];
            ++candidate.error_neighbor_count;
            candidate.error_neighbor_index = std::min(candidate.error_neighbor_index, j);
        }
    }
    read_evidence::ExactMatcher matcher;
    for (std::size_t i = 0; i < opportunities.size(); ++i) {
        const auto& repeat = opportunities[i].repeat;
        matcher.add(repeat.substr(repeat.size() - MINIMUM_PARTIAL_REPEAT), i);
    }
    matcher.build();
    result.counts["read_recruitment_error_hypotheses"] = hypotheses;
    result.counts["read_recruitment_opportunity_read_records"] = 0;
    result.counts["read_recruitment_opportunity_windows"] = 0;
    result.counts["read_recruitment_opportunity_resource_limited_records"] = 0;
    result.counts["read_recruitment_error_ambiguous_candidates"] = 0;
    result.counts["read_recruitment_selection_eligible_candidates"] = 0;
    auto scan = [&](const klibpp::KSeq& read) {
        ++result.counts["read_recruitment_opportunity_read_records"];
        if (read.qual.empty() || opportunities.empty() || !complete) return;
        if (read.seq.size() > config.maximum_read_length) {
            complete = false;
            ++result.counts["read_recruitment_opportunity_resource_limited_records"];
            return;
        }
        std::size_t matches = 0;
        auto oriented = [&](const std::string& sequence, const std::string& quality_string) {
            read_evidence::QualitySpans quality(quality_string, MINIMUM_QUALITY);
            matcher.scan(sequence, [&](std::size_t id, std::size_t spacer_start) {
                if (!complete) return;
                if (++matches > config.maximum_opportunity_matches_per_read) { complete = false; return; }
                auto& opportunity = opportunities[id];
                const auto spacer_end = spacer_start + opportunity.length;
                if (spacer_end > sequence.size() || MINIMUM_PARTIAL_REPEAT > sequence.size() - spacer_end) return;
                const auto start = spacer_start - MINIMUM_PARTIAL_REPEAT;
                const auto end = spacer_end + MINIMUM_PARTIAL_REPEAT;
                if (!quality.passes(start, end, MAXIMUM_EXPECTED_ERRORS)) return;
                if (opportunity.records == std::numeric_limits<std::uint64_t>::max()) {
                    opportunity.complete = false; complete = false; return;
                }
                ++opportunity.records;
                ++result.counts["read_recruitment_opportunity_windows"];
                const auto minimum = *std::min_element(quality_string.begin() + spacer_start,
                                                       quality_string.begin() + spacer_end);
                opportunity.minimum_quality = std::min(opportunity.minimum_quality, static_cast<int>(minimum) - 33);
            });
        };
        oriented(read.seq, read.qual);
        if (complete) oriented(unit_rc(read.seq), std::string(read.qual.rbegin(), read.qual.rend()));
        if (!complete) ++result.counts["read_recruitment_opportunity_resource_limited_records"];
    };
    // Even an empty/overflowed candidate set still validates the second input
    // pass. No partial parse failure can produce apparently complete evidence.
    read_linked::Reader first(files.at(0));
    std::unique_ptr<read_linked::Reader> second;
    if (files.size() == 2) second = std::make_unique<read_linked::Reader>(files[1]);
    klibpp::KSeq a, b;
    while (first.next(a)) {
        if (second && (!second->next(b) || read_linked::read_name(a.name) != read_linked::read_name(b.name)))
            throw std::runtime_error("Error-screen read pairs have different names or lengths");
        scan(a); if (second) scan(b);
    }
    if (second && second->next(b)) throw std::runtime_error("Error-screen read files contain different record counts");
    result.counts["read_recruitment_opportunity_complete"] = complete;
    for (auto i : eligible) {
        auto& candidate = result.candidates[i];
        const auto& opportunity = opportunities.at(bucket_ids.at({candidate.left_repeat, candidate.spacer.size()}));
        candidate.error_hypotheses = hypotheses;
        candidate.opportunity_records = opportunity.records;
        candidate.opportunity_min_spacer_quality = opportunity.minimum_quality == 94 ? -1 : opportunity.minimum_quality;
        candidate.opportunity_complete = complete && opportunity.complete;
        // A tandem-compatible witness cannot make a biological-allele error
        // rejection look stronger. Keep M and neighbors fixed before either
        // screen, but use only independently non-tandem units as the K bound.
        candidate.error_observed_units = candidate.non_tandem_support.sequence_units;
        candidate.error_log_threshold = hypotheses ? std::log(ERROR_SCREEN_ALPHA / static_cast<double>(hypotheses)) : 0;
        if (candidate.opportunity_min_spacer_quality >= 0)
            candidate.error_probability_upper = std::pow(10.0, -candidate.opportunity_min_spacer_quality / 10.0);
        if (candidate.tandem_ambiguous) {
            candidate.error_status = "not_evaluated_tandem_ambiguous";
        } else if (!candidate.opportunity_complete || opportunity.records < candidate.error_observed_units ||
            candidate.opportunity_min_spacer_quality < 0 || !hypotheses) {
            candidate.opportunity_complete = false;
            candidate.error_ambiguous = true;
            candidate.error_status = "opportunity_unresolved";
        } else if (!candidate.error_neighbor_count) {
            candidate.error_status = "no_supported_neighbor";
        } else {
            candidate.error_log_bound = independent_error_log_bound(opportunity.records,
                candidate.error_observed_units, candidate.error_probability_upper);
            candidate.error_ambiguous = candidate.error_log_bound >= candidate.error_log_threshold;
            candidate.error_status = candidate.error_ambiguous ? "error_ambiguous" : "rejects_independent_error_null";
        }
        candidate.selection_eligible = candidate.accepted && !candidate.error_ambiguous && !candidate.tandem_ambiguous;
        result.counts["read_recruitment_error_ambiguous_candidates"] += candidate.error_ambiguous;
        result.counts["read_recruitment_selection_eligible_candidates"] += candidate.selection_eligible;
    }
}

// This is an additive, graph-seeded local read graph. The exact spacer interior
// is its node identity. Repeat copies may differ by one observed substitution
// from a seed; no sequence consensus or endpoint-only allele collapse is used.
// Full/partial discovery requires the complete spacer and a terminal repeat
// prefix. Partial/partial evidence is restricted to an existing graph proposal.
template<class Arrays>
Result infer(const Arrays& arrays, const std::vector<std::string>& files,
             const Config& config = {}) {
    using namespace detail;
    if (files.empty() || files.size() > 2)
        throw std::invalid_argument("Read recruitment requires one or two read files");
    if (!config.maximum_candidates_per_family || !config.maximum_total_candidates ||
        !config.maximum_family_memberships || !config.maximum_seed_placements_per_read ||
        !config.maximum_candidates_per_read || !config.maximum_repeat_verifications_per_read ||
        !config.maximum_read_length || !config.maximum_opportunity_matches_per_read)
        throw std::invalid_argument("Read recruitment resource bounds must be positive");
    Result result;
    result.counts = {{"read_recruitment_applied", 1}, {"read_recruitment_records", 0},
        {"read_recruitment_records_without_quality", 0}, {"read_recruitment_input_units", 0},
        {"read_recruitment_resource_limited_units", 0}, {"read_recruitment_resource_limited_families", 0},
        {"read_recruitment_global_resource_limit", 0}, {"read_recruitment_candidates", 0},
        {"read_recruitment_accepted_candidates", 0}, {"read_recruitment_boundary_conflicts", 0},
        {"read_recruitment_ambiguous_owner_candidates", 0}, {"read_recruitment_capped_candidates", 0},
        {"read_recruitment_distinct_candidate_keys", 0}, {"read_recruitment_family_memberships", 0},
        {"read_recruitment_tandem_explained_sources", 0}, {"read_recruitment_tandem_ambiguous_candidates", 0}};
    std::map<std::string, Owners> seed_owners;
    std::map<Pair, std::set<std::size_t>> graph_pairs;
    std::map<std::string, std::set<std::string>> graph_cores;
    for (std::size_t family = 0; family < arrays.size(); ++family) {
        std::set<std::string> repeats;
        for (const auto& pair : arrays[family].entries) {
            if (!dna(pair.first) || !dna(pair.second))
                throw std::invalid_argument("Read recruitment graph pairs must contain only ACGT");
            graph_cores[core(pair)].insert(canonical(pair.second));
            graph_pairs[pair].insert(family);
            graph_pairs[{rc(pair.first), rc(pair.second)}].insert(family);
            if (pair.first.size() >= MINIMUM_REPEAT && pair.first.size() <= MAXIMUM_REPEAT)
                repeats.insert(canonical(pair.first));
        }
        if (repeats.empty()) continue;
        const auto reference = *std::max_element(repeats.begin(), repeats.end(),
            [](const auto& a, const auto& b) {
                return a.size() != b.size() ? a.size() < b.size() : a > b;
            });
        for (const auto& repeat : repeats) {
            const auto reverse = rc(repeat);
            const auto forward_distance = edit_distance(reference, repeat);
            const auto reverse_distance = edit_distance(reference, reverse);
            if (forward_distance <= reverse_distance) {
                seed_owners[repeat].insert({family, false});
                seed_owners[reverse].insert({family, true});
            }
            if (reverse_distance <= forward_distance) {
                seed_owners[repeat].insert({family, true});
                seed_owners[reverse].insert({family, false});
            }
        }
    }
    std::vector<Seed> seeds;
    std::vector<Anchor> anchors;
    std::map<Owner, std::vector<std::size_t>> owner_seeds;
    read_evidence::ExactMatcher seed_matcher;
    for (const auto& entry : seed_owners) {
        const auto seed = seeds.size();
        seeds.push_back({entry.first, entry.second});
        for (const auto& owner : entry.second) owner_seeds[owner].push_back(seed);
        // Disjoint positional seeds guarantee that one survives one SNP,
        // including short 20/21-base repeats for which k is ten rather than eleven.
        const auto length = std::min(SEED_LENGTH, entry.first.size() / 2);
        for (const auto offset : {std::size_t{0}, entry.first.size() - length}) {
            seed_matcher.add(entry.first.substr(offset, length), anchors.size());
            anchors.push_back({seed, offset, length});
        }
    }
    seed_matcher.build();
    result.counts["read_recruitment_oriented_repeat_seeds"] = seeds.size();

    std::vector<GraphPattern> graph_patterns;
    read_evidence::ExactMatcher graph_matcher;
    for (const auto& entry : graph_pairs) {
        const auto& pair = entry.first;
        const Pair reverse{rc(pair.first), rc(pair.second)};
        if (reverse < pair || pair.first.size() < MINIMUM_REPEAT || pair.first.size() > MAXIMUM_REPEAT ||
            pair.second.size() < MINIMUM_SPACER || pair.second.size() > MAXIMUM_SPACER) continue;
        Candidate candidate;
        candidate.family_owners = entry.second;
        candidate.left_repeat = candidate.right_repeat = pair.first;
        candidate.spacer = pair.second;
        candidate.evidence_kind = "graph_junction";
        candidate.graph_proposed = true;
        candidate.left_observed_bases = candidate.right_observed_bases = MINIMUM_PARTIAL_REPEAT;
        for (bool reversed : {false, true}) {
            const auto repeat = reversed ? reverse.first : pair.first;
            const auto spacer = reversed ? reverse.second : pair.second;
            const auto sequence = repeat.substr(repeat.size() - MINIMUM_PARTIAL_REPEAT) +
                spacer + repeat.substr(0, MINIMUM_PARTIAL_REPEAT);
            graph_matcher.add(sequence, graph_patterns.size());
            graph_patterns.push_back({candidate, repeat, spacer, reversed});
        }
    }
    graph_matcher.build();
    std::map<Key, Accumulator> evidence;
    std::set<Key> all_keys;
    std::map<Key, std::set<std::size_t>> all_owners;
    std::map<std::size_t, std::set<Key>> family_keys;
    std::set<std::size_t> overloaded_families;
    std::size_t family_memberships = 0;
    std::size_t observed_candidate_keys = 0;
    bool global_overflow = false;

    auto scan = [&](const klibpp::KSeq& read, std::map<Key, Pending>& pending) {
        if (read.seq.size() > config.maximum_read_length) return false;
        if (read.qual.empty()) return true;
        std::size_t verifications = 0;
        bool overloaded = false;
        auto add = [&](Candidate candidate, std::size_t offset,
                       const std::map<Pair, std::size_t>& pair_offsets) {
            if (candidate.family_owners.empty()) return;
            const auto candidate_key = key(candidate);
            const auto found = pending.emplace(candidate_key, Pending{candidate, {}, {}});
            auto& destination = found.first->second;
            destination.candidate.family_owners.insert(candidate.family_owners.begin(), candidate.family_owners.end());
            destination.offsets.insert(offset);
            for (const auto& pair_offset : pair_offsets)
                destination.pair_offsets[pair_offset.first].insert(pair_offset.second);
            if (pending.size() > config.maximum_candidates_per_read) overloaded = true;
        };
        auto discovery = [&](const std::string& sequence, const std::string& quality_string) {
            read_evidence::QualitySpans quality(quality_string, MINIMUM_QUALITY);
            std::set<std::pair<std::size_t, std::size_t>> placements;
            std::map<std::pair<std::size_t, std::size_t>, Hit> hits;
            seed_matcher.scan(sequence, [&](std::size_t anchor_id, std::size_t end) {
                if (overloaded) return;
                const auto& anchor = anchors[anchor_id];
                if (end < anchor.length + anchor.offset) return;
                const auto start = end - anchor.length - anchor.offset;
                const auto& seed = seeds[anchor.seed];
                if (start > sequence.size() || seed.sequence.size() > sequence.size() - start) return;
                if (!placements.emplace(anchor.seed, start).second) return;
                if (placements.size() > config.maximum_seed_placements_per_read ||
                    ++verifications > config.maximum_repeat_verifications_per_read) { overloaded = true; return; }
                if (!verifies(sequence, start, seed.sequence.size(), seed.sequence)) return;
                const auto stop = start + seed.sequence.size();
                const auto hit = hits.emplace(std::make_pair(start, stop),
                    Hit{start, stop, sequence.substr(start, seed.sequence.size()), {}});
                hit.first->second.owners.insert(seed.owners.begin(), seed.owners.end());
            });
            if (overloaded) return;
            for (const auto& left_item : hits) {
                const auto& left = left_item.second;
                for (auto right_item = hits.lower_bound({left.end + MINIMUM_SPACER, 0});
                     right_item != hits.end() && right_item->second.start <= left.end + MAXIMUM_SPACER; ++right_item) {
                    const auto& right = right_item->second;
                    const auto owners = common(left.owners, right.owners);
                    if (owners.empty() || !quality.passes(left.start, right.end, MAXIMUM_EXPECTED_ERRORS)) continue;
                    Candidate candidate;
                    candidate.family_owners = families(owners);
                    candidate.left_repeat = left.sequence;
                    candidate.spacer = sequence.substr(left.end, right.start - left.end);
                    candidate.right_repeat = right.sequence;
                    if (!dna(candidate.spacer)) continue;
                    candidate.evidence_kind = "full_pair";
                    auto offset = left.end - MINIMUM_PARTIAL_REPEAT;
                    const auto reverse_offset = sequence.size() - right.start - MINIMUM_PARTIAL_REPEAT;
                    const std::map<Pair, std::size_t> pair_offsets{
                        {{left.sequence, candidate.spacer}, offset},
                        {{rc(right.sequence), rc(candidate.spacer)}, reverse_offset}};
                    const auto forward = std::make_tuple(candidate.left_repeat, candidate.spacer, candidate.right_repeat);
                    const auto reverse = std::make_tuple(rc(candidate.right_repeat), rc(candidate.spacer), rc(candidate.left_repeat));
                    if (reverse < forward) {
                        std::tie(candidate.left_repeat, candidate.spacer, candidate.right_repeat) = reverse;
                        offset = reverse_offset;
                    }
                    candidate.left_observed_bases = candidate.left_repeat.size();
                    candidate.right_observed_bases = candidate.right_repeat.size();
                    add(std::move(candidate), offset, pair_offsets);
                    if (overloaded) return;
                }
                std::set<std::size_t> right_seeds;
                for (const auto& owner : left.owners) {
                    const auto& ids = owner_seeds.at(owner);
                    right_seeds.insert(ids.begin(), ids.end());
                }
                for (std::size_t length = MINIMUM_SPACER; length <= MAXIMUM_SPACER; ++length) {
                    const auto start = left.end + length;
                    if (start > sequence.size()) break;
                    const auto available = sequence.size() - start;
                    if (available < MINIMUM_PARTIAL_REPEAT) break;
                    for (auto id : right_seeds) {
                        const auto& seed = seeds[id];
                        if (available >= seed.sequence.size()) continue;
                        if (++verifications > config.maximum_repeat_verifications_per_read) { overloaded = true; return; }
                        if (!verifies(sequence, start, available, seed.sequence) ||
                            !quality.passes(left.start, sequence.size(), MAXIMUM_EXPECTED_ERRORS)) continue;
                        Candidate candidate;
                        candidate.family_owners = families(common(left.owners, seed.owners));
                        candidate.left_repeat = left.sequence;
                        candidate.spacer = sequence.substr(left.end, length);
                        if (!dna(candidate.spacer)) continue;
                        candidate.right_repeat = sequence.substr(start, available);
                        candidate.right_complete = false;
                        candidate.left_observed_bases = left.sequence.size();
                        candidate.right_observed_bases = available;
                        candidate.evidence_kind = "full_partial";
                        const auto offset = left.end - MINIMUM_PARTIAL_REPEAT;
                        const std::map<Pair, std::size_t> pair_offsets{{{candidate.left_repeat, candidate.spacer}, offset}};
                        add(std::move(candidate), offset, pair_offsets);
                        if (overloaded) return;
                    }
                }
            }
        };
        discovery(read.seq, read.qual);
        if (!overloaded) discovery(unit_rc(read.seq), std::string(read.qual.rbegin(), read.qual.rend()));
        if (overloaded) return false;

        read_evidence::QualitySpans quality(read.qual, MINIMUM_QUALITY);
        graph_matcher.scan(read.seq, [&](std::size_t id, std::size_t end) {
            if (overloaded) return;
            const auto& pattern = graph_patterns[id];
            const auto span = pattern.spacer.size() + 2 * MINIMUM_PARTIAL_REPEAT;
            const auto start = end - span;
            const auto spacer_start = start + MINIMUM_PARTIAL_REPEAT;
            const auto spacer_end = spacer_start + pattern.spacer.size();
            const auto left_length = std::min(spacer_start, pattern.repeat.size());
            const auto right_length = std::min(read.seq.size() - spacer_end, pattern.repeat.size());
            // Check every available repeat base, not just convenient anchors.
            if (read.seq.compare(spacer_start - left_length, left_length, pattern.repeat,
                    pattern.repeat.size() - left_length, left_length) != 0 ||
                read.seq.compare(spacer_end, right_length, pattern.repeat, 0, right_length) != 0 ||
                !quality.passes(spacer_start - left_length, spacer_end + right_length, MAXIMUM_EXPECTED_ERRORS)) return;
            const auto reverse_offset = read.seq.size() - end;
            const std::map<Pair, std::size_t> pair_offsets{
                {{pattern.repeat, pattern.spacer}, start},
                {{rc(pattern.repeat), rc(pattern.spacer)}, reverse_offset}};
            auto candidate = pattern.candidate;
            candidate.left_observed_bases = pattern.reversed ? right_length : left_length;
            candidate.right_observed_bases = pattern.reversed ? left_length : right_length;
            add(std::move(candidate), pattern.reversed ? reverse_offset : start, pair_offsets);
        });
        return !overloaded;
    };

    read_linked::Reader first(files[0]);
    std::unique_ptr<read_linked::Reader> second;
    if (files.size() == 2) second = std::make_unique<read_linked::Reader>(files[1]);
    klibpp::KSeq a, b;
    while (first.next(a)) {
        if (second && (!second->next(b) || read_linked::read_name(a.name) != read_linked::read_name(b.name)))
            throw std::runtime_error("Recruitment read pairs have different names or lengths");
        ++result.counts["read_recruitment_input_units"];
        result.counts["read_recruitment_records"] += second ? 2 : 1;
        result.counts["read_recruitment_records_without_quality"] += a.qual.empty() + (second && b.qual.empty());
        if (global_overflow) continue; // Still validate every input record/pair.
        std::map<Key, Pending> pending;
        bool within_bounds = scan(a, pending);
        if (second && !scan(b, pending)) within_bounds = false;
        if (!within_bounds) { ++result.counts["read_recruitment_resource_limited_units"]; continue; }
        if (pending.empty()) continue;
        const auto unit = unit_identity(a.seq, second ? b.seq : "", bool(second));
        for (auto& item : pending) {
            if (all_keys.insert(item.first).second && all_keys.size() > config.maximum_total_candidates) {
                global_overflow = true; break;
            }
            auto& candidate = item.second.candidate;
            bool family_limited = false;
            for (auto family : candidate.family_owners) {
                if (all_owners[item.first].insert(family).second) ++family_memberships;
                if (family_memberships > config.maximum_family_memberships) { global_overflow = true; break; }
                if (overloaded_families.count(family)) { family_limited = true; continue; }
                auto& keys = family_keys[family];
                keys.insert(item.first);
                if (keys.size() > config.maximum_candidates_per_family) {
                    overloaded_families.insert(family);
                    family_limited = true;
                    for (auto old = evidence.begin(); old != evidence.end();) {
                        if (old->second.candidate.family_owners.count(family)) old = evidence.erase(old);
                        else ++old;
                    }
                }
            }
            if (global_overflow) break;
            for (auto family : all_owners.at(item.first))
                family_limited |= overloaded_families.count(family) != 0;
            if (family_limited) { evidence.erase(item.first); continue; }
            const auto found = evidence.emplace(item.first, Accumulator{candidate, {}, {}, {}, {}, {}});
            auto& accumulated = found.first->second;
            accumulated.candidate.family_owners.insert(candidate.family_owners.begin(), candidate.family_owners.end());
            if (accumulated.units.size() < SUPPORT_CAP) accumulated.units.insert(unit);
            for (auto offset : item.second.offsets) {
                if (accumulated.offsets.size() >= SUPPORT_CAP) break;
                accumulated.offsets.insert(offset);
            }
            for (const auto& pair_offsets : item.second.pair_offsets) {
                auto& destination = accumulated.pair_offsets[pair_offsets.first];
                for (auto offset : pair_offsets.second) {
                    if (destination.size() >= SUPPORT_CAP) break;
                    destination.insert(offset);
                }
            }
        }
        observed_candidate_keys = all_keys.size();
        if (global_overflow) { evidence.clear(); all_keys.clear(); all_owners.clear(); family_keys.clear(); }
    }
    if (second && second->next(b))
        throw std::runtime_error("Recruitment read files contain different record counts");
    result.counts["read_recruitment_global_resource_limit"] = global_overflow;
    result.counts["read_recruitment_resource_limited_families"] = overloaded_families.size();
    result.counts["read_recruitment_family_memberships"] = family_memberships;
    result.counts["read_recruitment_distinct_candidate_keys"] = global_overflow ? observed_candidate_keys : all_keys.size();
    std::vector<Accumulator*> ordered;
    for (auto& item : evidence) {
        auto& accumulated = item.second;
        auto& candidate = accumulated.candidate;
        candidate.family_owners = all_owners.at(item.first);
        for (auto owner : candidate.family_owners)
            candidate.resource_limited |= overloaded_families.count(owner) != 0;
        candidate.core = core({candidate.left_repeat, candidate.spacer});
        const auto spacer = canonical(candidate.spacer);
        for (const auto& pair : oriented_pairs(candidate)) {
            const auto old = graph_cores.find(core(pair));
            if (old != graph_cores.end() && !old->second.count(spacer)) candidate.boundary_conflict = true;
        }
        candidate.tandem = tandem_evidence(candidate);
        candidate.tandem_status = "provenance_only";
        if (!candidate.resource_limited) ordered.push_back(&accumulated);
    }
    const auto raw_count = ordered.size();
    std::map<Pair, Accumulator> pooled;
    for (std::size_t source = 0; source < raw_count && !global_overflow; ++source) {
        const auto& raw = *ordered[source];
        for (const auto& pair_offsets : raw.pair_offsets) {
            const auto& pair = pair_offsets.first;
            auto& aggregate = pooled[pair];
            auto& candidate = aggregate.candidate;
            if (all_keys.size() + pooled.size() > config.maximum_total_candidates) {
                global_overflow = true; break;
            }
            const bool forward = pair == Pair{raw.candidate.left_repeat, raw.candidate.spacer};
            const auto next_repeat = forward ? raw.candidate.right_repeat : rc(raw.candidate.left_repeat);
            const auto left_observed = forward ? raw.candidate.left_observed_bases : raw.candidate.right_observed_bases;
            const auto right_observed = forward ? raw.candidate.right_observed_bases : raw.candidate.left_observed_bases;
            const bool right_complete = forward ? raw.candidate.right_complete : raw.candidate.left_complete;
            if (candidate.left_repeat.empty()) {
                candidate.left_repeat = pair.first;
                candidate.spacer = pair.second;
                candidate.evidence_kind = "full_pair";
                candidate.pair_only = true;
                candidate.left_observed_bases = left_observed;
                candidate.right_observed_bases = right_observed;
            }
            candidate.source_candidate_indices.push_back(source);
            // An unequal-flank triple can expose a competing cut only on its
            // opposite side. Pooling its shorter circular interpretation must
            // not erase that observed boundary ambiguity.
            candidate.boundary_conflict |= raw.candidate.boundary_conflict;
            candidate.resource_limited |= raw.candidate.resource_limited;
            candidate.observed_right_repeats.push_back(next_repeat);
            candidate.left_observed_bases = std::min(candidate.left_observed_bases, left_observed);
            candidate.right_observed_bases = std::min(candidate.right_observed_bases, right_observed);
            candidate.right_complete &= right_complete;
            if (raw.candidate.evidence_kind == "graph_junction") candidate.evidence_kind = "graph_junction";
            else if (raw.candidate.evidence_kind == "full_partial" && candidate.evidence_kind != "graph_junction")
                candidate.evidence_kind = "full_partial";
            for (auto owner : raw.candidate.family_owners) {
                if (candidate.family_owners.insert(owner).second) ++family_memberships;
                if (family_memberships > config.maximum_family_memberships) { global_overflow = true; break; }
            }
            if (global_overflow) break;
            for (const auto& unit : raw.units) {
                if (aggregate.units.size() >= SUPPORT_CAP) break;
                aggregate.units.insert(unit);
            }
            for (auto offset : pair_offsets.second) {
                if (aggregate.offsets.size() >= SUPPORT_CAP) break;
                aggregate.offsets.insert(offset);
            }
            if (!raw.candidate.tandem.explained) {
                for (const auto& unit : raw.units) {
                    if (aggregate.non_tandem_units.size() >= SUPPORT_CAP) break;
                    aggregate.non_tandem_units.insert(unit);
                }
                for (auto offset : pair_offsets.second) {
                    if (aggregate.non_tandem_offsets.size() >= SUPPORT_CAP) break;
                    aggregate.non_tandem_offsets.insert(offset);
                }
            }
        }
    }
    result.counts["read_recruitment_global_resource_limit"] = global_overflow;
    result.counts["read_recruitment_family_memberships"] = family_memberships;
    result.counts["read_recruitment_pair_aggregates"] = global_overflow ? 0 : pooled.size();
    if (global_overflow) { annotate_error_uncertainty(result, files, config); return result; }
    for (auto& item : pooled) {
        auto& candidate = item.second.candidate;
        auto& repeats = candidate.observed_right_repeats;
        std::sort(repeats.begin(), repeats.end());
        repeats.erase(std::unique(repeats.begin(), repeats.end()), repeats.end());
        candidate.right_repeat = repeats.front(); // Explicit representative; never licenses reverse output.
        ordered.push_back(&item.second);
    }
    std::map<std::size_t, std::set<std::string>> supported_spacers;
    for (auto* pointer : ordered) {
        auto& accumulated = *pointer;
        auto& candidate = accumulated.candidate;
        candidate.support = {accumulated.units.size(), accumulated.offsets.size(),
            accumulated.units.size() >= SUPPORT_CAP, accumulated.offsets.size() >= SUPPORT_CAP};
        candidate.non_tandem_support = {accumulated.non_tandem_units.size(), accumulated.non_tandem_offsets.size(),
            accumulated.non_tandem_units.size() >= SUPPORT_CAP, accumulated.non_tandem_offsets.size() >= SUPPORT_CAP};
        candidate.core = core({candidate.left_repeat, candidate.spacer});
        const auto spacer = canonical(candidate.spacer);
        for (const auto& pair : oriented_pairs(candidate)) {
            const auto graph = graph_pairs.find(pair);
            if (graph != graph_pairs.end())
                for (auto owner : candidate.family_owners)
                    if (graph->second.count(owner)) candidate.graph_proposed = true;
            const auto old = graph_cores.find(core(pair));
            if (old != graph_cores.end() && !old->second.count(spacer)) candidate.boundary_conflict = true;
        }
        if (candidate.support.accepted() && !candidate.boundary_conflict && !candidate.resource_limited)
            for (auto owner : candidate.family_owners) supported_spacers[owner].insert(spacer);
    }
    for (auto* pointer : ordered) {
        auto& candidate = pointer->candidate;
        // Every possible owner must have independently supported diversity.
        // Ambiguous ownership is retained; callers must not invent a locus.
        candidate.distinct_supported_family_spacers = std::numeric_limits<std::size_t>::max();
        for (auto owner : candidate.family_owners)
            candidate.distinct_supported_family_spacers = std::min(
                candidate.distinct_supported_family_spacers, supported_spacers[owner].size());
        if (candidate.family_owners.empty()) candidate.distinct_supported_family_spacers = 0;
        // Raw triplets remain provenance only. All callable claims use the
        // explicit directed pair aggregate, whose offsets belong to that exact
        // orientation; canonical-triplet offsets cannot license its reverse.
        candidate.accepted = candidate.pair_only && candidate.support.accepted() && !candidate.boundary_conflict &&
            !candidate.resource_limited && candidate.distinct_supported_family_spacers >= MINIMUM_FAMILY_SPACERS;
        candidate.tandem_ambiguous = candidate.accepted && !candidate.non_tandem_support.accepted();
        if (candidate.pair_only) candidate.tandem_status = !candidate.accepted ? "not_admissible" :
            candidate.tandem_ambiguous ? "tandem_ambiguous" : "independent_non_tandem_support";
        candidate.reason = candidate.resource_limited ? "family_resource_limit" :
            candidate.boundary_conflict ? "existing_core_alternative_boundary" :
            candidate.support.sequence_units < MINIMUM_SEQUENCE_UNITS ? "insufficient_sequence_units" :
            candidate.support.start_offsets < MINIMUM_START_OFFSETS ? "insufficient_start_offsets" :
            candidate.distinct_supported_family_spacers < MINIMUM_FAMILY_SPACERS ? "insufficient_family_diversity" :
            candidate.family_owners.size() > 1 ? "supported_ambiguous_family" : "supported";
        if (!candidate.pair_only && candidate.support.accepted() && !candidate.boundary_conflict &&
            !candidate.resource_limited && candidate.distinct_supported_family_spacers >= MINIMUM_FAMILY_SPACERS)
            candidate.reason = "triplet_evidence_only";
        if (candidate.accepted && candidate.pair_only) candidate.reason += "_pair_aggregate";
        result.counts["read_recruitment_accepted_candidates"] += candidate.accepted;
        result.counts["read_recruitment_boundary_conflicts"] += candidate.boundary_conflict;
        result.counts["read_recruitment_ambiguous_owner_candidates"] += candidate.family_owners.size() > 1;
        result.counts["read_recruitment_capped_candidates"] += candidate.support.units_saturated || candidate.support.offsets_saturated;
        result.counts["read_recruitment_tandem_explained_sources"] += !candidate.pair_only && candidate.tandem.explained;
        result.counts["read_recruitment_tandem_ambiguous_candidates"] += candidate.tandem_ambiguous;
        result.candidates.push_back(std::move(candidate));
    }
    result.counts["read_recruitment_candidates"] = result.candidates.size();
    annotate_error_uncertainty(result, files, config);
    return result;
}

}} // namespace Zuck::read_recruitment
#endif
