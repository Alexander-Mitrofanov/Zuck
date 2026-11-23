#ifndef ZUCK_READ_FIRST_CLIPPED_JUNCTION_H
#define ZUCK_READ_FIRST_CLIPPED_JUNCTION_H

#include "crispr_array/read_first_seeds.h"
#include "crispr_array/read_recruitment.h"

namespace Zuck { namespace read_first { namespace clipped {

using Key = std::pair<std::size_t, std::string>;
using Unit = read_recruitment::detail::Unit;

struct Candidate {
    std::size_t strand, begin, end, tip;
};

// A complete repeat establishes one junction, and a 4–9 base prefix or suffix
// of its partner extends all the way to the physical read end. Require evidence
// for this length from other spacers, on other fragments, in the same family.
inline bool supports_length(const Family& family, const std::string& spacer,
                            std::size_t tip, std::uint64_t fragment) {
    if (family.templates < 5 || family.contexts < 4) return false;
    const auto focal = read_linked::canonical(spacer);
    std::map<std::string, std::set<std::uint64_t>> near;
    std::set<std::string> others;
    for (const auto& b : family.boundaries) {
        if (b.fragment == fragment || !seeds::different_interiors(b.spacer, focal)) continue;
        others.insert(b.spacer);
        const auto l = b.spacer.size(), n = spacer.size();
        const bool compatible = tip <= 6 ? l == n : l + 1 >= n && n + 1 >= l;
        if (compatible) near[b.spacer].insert(b.fragment);
    }
    // Short hooks occur by chance; use only the modal spacer length when it
    // accounts for at least three quarters of independent observed contexts.
    if (near.size() < 2 || (tip <= 6 && 4 * near.size() < 3 * others.size())) return false;
    for (auto a = near.begin(); a != near.end(); ++a) {
        for (auto b = near.begin(); b != a; ++b) {
            if (seeds::different_interiors(a->first, b->first) &&
                (a->second.size() > 1 || b->second.size() > 1 || *a->second.begin() != *b->second.begin()))
                return true;
        }
    }
    return false;
}

inline bool near_known(std::string spacer, const std::set<std::string>& pool) {
    for (std::size_t i = 0; i <= spacer.size(); ++i) {
        if (i < spacer.size() && pool.count(spacer.substr(0,i) + spacer.substr(i+1))) return true;
        for (const char b : {'A','C','G','T'}) {
            if (pool.count(spacer.substr(0,i) + b + spacer.substr(i))) return true;
            if (i < spacer.size() && b != spacer[i] &&
                pool.count(spacer.substr(0,i) + b + spacer.substr(i+1))) return true;
        }
    }
    return false;
}

class Collector {
    const std::vector<Family>& families_;
    std::map<Key,std::set<Unit>> units_;
    std::size_t memberships_ = 0;
    bool limited_ = false;

public:
    explicit Collector(const std::vector<Family>& families) : families_(families) {}

    // Reuse the strict scan's strand-aware complete-repeat hits. Every candidate
    // spacer is present on one read; opposite mates can never count twice.
    template<class Sites, class Strands, class Placements>
    void scan(const klibpp::KSeq& read, const Sites& left, const Sites& right,
              const Strands& strands, const Placements& strict, std::vector<Candidate>& out) const {
        if (limited_) return;
        const auto n = read.seq.size();
        std::vector<Candidate> candidates;
        const auto clipped_end = [&](const Sites& hits, bool downstream) {
            for (const auto& hit : hits) {
                const auto strand = hit.first.first;
                const auto& repeat = strands[strand].repeat;
                const auto pos = hit.first.second;
                if (hit.second != repeat.size()) continue;
                for (std::size_t tip = 4; tip < 10; ++tip) {
                    if (tip > n || (downstream ? n-tip <= pos : pos <= tip)) continue;
                    const auto begin = downstream ? pos : tip;
                    const auto end = downstream ? n - tip : pos;
                    if (end - begin < read_recruitment::MINIMUM_SPACER ||
                        end - begin > read_recruitment::MAXIMUM_SPACER) continue;
                    if (read.seq.compare(downstream ? end : 0, tip, repeat,
                                         downstream ? 0 : repeat.size()-tip, tip) != 0) continue;
                    candidates.push_back({strand,begin,end,tip});
                    if (candidates.size() > 64) return false;
                }
            }
            return true;
        };
        if (!clipped_end(left,true) || !clipped_end(right,false)) return;
        read_evidence::QualitySpans quality(read.qual, 15);
        for (const auto& c : candidates) {
            const auto overlaps = [&](const auto& p) {
                return c.begin < p.end && p.begin < c.end;
            };
            if (std::any_of(strict.begin(),strict.end(),overlaps) ||
                std::any_of(candidates.begin(),candidates.end(),[&](const auto& p) {
                    return &p != &c && overlaps(p);
                })) continue;
            const auto& strand = strands[c.strand];
            const auto flank = strand.repeat.size();
            const auto left_span = c.begin == c.tip ? c.tip : flank;
            const auto right_span = n-c.end == c.tip ? c.tip : flank;
            if (!quality.passes(c.begin-left_span,c.end+right_span,0.6)) continue;
            // Both the six terminal spacer bases and the short repeat must be Q20.
            const bool start_clip = c.begin == c.tip;
            const auto lo = start_clip ? 0 : c.end-6;
            const auto hi = start_clip ? c.begin+6 : n;
            if (std::any_of(read.qual.begin()+lo,read.qual.begin()+hi,[](char q) { return q < 53; }))
                continue;
            out.push_back(c);
        }
    }

    template<class Strands>
    void store(const klibpp::KSeq& read, const Strands& strands,
               const std::vector<Candidate>& candidates, const Unit& unit, std::uint64_t fragment) {
        if (limited_) return;
        for (const auto& c : candidates) {
            const auto& strand = strands[c.strand];
            auto spacer = read.seq.substr(c.begin,c.end-c.begin);
            if (!read_recruitment::detail::dna(spacer) ||
                !seeds::admissible_interior(strand.repeat,spacer)) continue;
            if (strand.reverse) spacer = read_linked::reverse_complement(spacer);
            const auto& family = families_[strand.family];
            if (family.repeat == read_linked::reverse_complement(family.repeat) ||
                !supports_length(family,spacer,c.tip,fragment)) continue;
            memberships_ += units_[{strand.family,spacer}].insert(unit).second;
            if (memberships_ > 250000 || units_.size() > 65536) {
                limited_ = true;
                units_.clear();
                return;
            }
        }
    }

    const auto& candidates() const { return units_; }
    bool limited() const { return limited_; }
};

}}} // namespace Zuck::read_first::clipped
#endif
