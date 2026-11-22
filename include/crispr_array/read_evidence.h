#ifndef ZUCK_READ_EVIDENCE_H
#define ZUCK_READ_EVIDENCE_H

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <queue>
#include <stdexcept>
#include <string>
#include <vector>

namespace Zuck { namespace read_evidence {

inline int base(char nucleotide) {
    switch (nucleotide) {
        case 'A': return 0;
        case 'C': return 1;
        case 'G': return 2;
        case 'T': return 3;
        default: return -1;
    }
}

// Shared Aho–Corasick matcher for exact DNA evidence. Add both orientations
// when needed, then build once before scanning any reads. Match ends are
// exclusive; ambiguous read bases break a match instead of bridging it.
class ExactMatcher {
    // Compact node indices avoid doubling the automaton footprint on 64-bit
    // hosts; pattern IDs and read coordinates still use the full size_t range.
    using State = std::uint32_t;
    static constexpr State missing = std::numeric_limits<State>::max();
    struct Node {
        std::array<State, 4> next{{missing, missing, missing, missing}};
        State failure = 0;
        State output = missing;
        std::vector<std::size_t> ids;
    };
    std::vector<Node> nodes_{1};
    bool built_ = false;

public:
    void add(const std::string& sequence, std::size_t id) {
        if (built_) throw std::logic_error("Cannot add patterns to a built evidence matcher");
        if (sequence.empty()) throw std::invalid_argument("Empty read evidence pattern");
        State state = 0;
        for (char nucleotide : sequence) {
            const int symbol = base(nucleotide);
            if (symbol < 0) throw std::invalid_argument("Read evidence patterns must contain only ACGT");
            if (nodes_[state].next[symbol] == missing) {
                if (nodes_.size() >= missing)
                    throw std::length_error("Read evidence matcher exceeds node index capacity");
                nodes_[state].next[symbol] = static_cast<State>(nodes_.size());
                nodes_.emplace_back();
            }
            state = nodes_[state].next[symbol];
        }
        nodes_[state].ids.push_back(id);
    }

    void build() {
        if (built_) throw std::logic_error("Evidence matcher was already built");
        std::queue<State> pending;
        for (auto& child : nodes_[0].next) {
            if (child == missing) child = 0;
            else pending.push(child);
        }
        while (!pending.empty()) {
            const auto state = pending.front();
            pending.pop();
            const auto failure = nodes_[state].failure;
            nodes_[state].output = nodes_[failure].ids.empty() ? nodes_[failure].output : failure;
            for (std::size_t symbol = 0; symbol < 4; ++symbol) {
                auto& child = nodes_[state].next[symbol];
                if (child == missing) child = nodes_[failure].next[symbol];
                else {
                    nodes_[child].failure = nodes_[failure].next[symbol];
                    pending.push(child);
                }
            }
        }
        built_ = true;
    }

    template<class OnMatch>
    void scan(const std::string& sequence, OnMatch on_match) const {
        if (!built_) throw std::logic_error("Build the evidence matcher before scanning");
        State state = 0;
        for (std::size_t position = 0; position < sequence.size(); ++position) {
            const int symbol = base(sequence[position]);
            if (symbol < 0) { state = 0; continue; }
            state = nodes_[state].next[symbol];
            for (auto node = state; node != missing; node = nodes_[node].output)
                for (auto id : nodes_[node].ids) on_match(id, position + 1);
        }
    }
};

// Lazily index one read's quality string: reads without a relevant sequence
// match allocate nothing, while repeated span checks take constant time.
// Reader validates Phred+33 input before constructing this helper.
class QualitySpans {
    const std::string& quality_;
    int minimum_quality_;
    std::vector<double> errors_;
    std::vector<std::size_t> low_quality_;

    static const std::array<double, 94>& error_probabilities() {
        static const auto values = [] {
            std::array<double, 94> result{};
            for (std::size_t q = 0; q < result.size(); ++q)
                result[q] = std::pow(10.0, -static_cast<double>(q) / 10.0);
            return result;
        }();
        return values;
    }

    void index() {
        errors_.resize(quality_.size() + 1);
        low_quality_.resize(quality_.size() + 1);
        const auto& probability = error_probabilities();
        for (std::size_t position = 0; position < quality_.size(); ++position) {
            const int q = static_cast<unsigned char>(quality_[position]) - 33;
            if (q < 0 || q >= static_cast<int>(probability.size()))
                throw std::invalid_argument("Invalid Phred+33 read quality");
            errors_[position + 1] = errors_[position] + probability[q];
            low_quality_[position + 1] = low_quality_[position] + (q < minimum_quality_);
        }
    }

public:
    QualitySpans(const std::string& quality, int minimum_quality)
        : quality_(quality), minimum_quality_(minimum_quality) {}

    bool passes(std::size_t start, std::size_t end, double maximum_errors) {
        if (start >= end || end > quality_.size())
            throw std::out_of_range("Invalid read quality span");
        if (errors_.empty()) index();
        if (low_quality_[end] != low_quality_[start]) return false;
        double expected_errors = errors_[end] - errors_[start];
        // Each addition within this span can round at the magnitude of the
        // entire prefix. A fixed tolerance is insufficient for long reads with
        // a large error sum before an otherwise good span. Bound those rounding
        // steps and subtract again directly whenever they could cross the cut.
        const double rounding_bound = std::numeric_limits<double>::epsilon() *
            (static_cast<double>(end - start) + 2.0) *
            (std::abs(errors_[end]) + std::abs(errors_[start]) + std::abs(maximum_errors));
        if (std::abs(expected_errors - maximum_errors) <= std::max(1e-10, rounding_bound)) {
            expected_errors = 0;
            const auto& probability = error_probabilities();
            for (auto position = start; position < end; ++position)
                expected_errors += probability[static_cast<unsigned char>(quality_[position]) - 33];
        }
        return expected_errors <= maximum_errors;
    }
};

}} // namespace Zuck::read_evidence
#endif
