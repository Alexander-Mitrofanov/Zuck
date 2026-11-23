#ifndef ZUCK_READ_FIRST_SPARSE_H
#define ZUCK_READ_FIRST_SPARSE_H

#include "crispr_array/read_first.h"

namespace Zuck { namespace read_first { namespace sparse_junctions {

using Counts = std::map<std::string, std::size_t>;
using Key = extraction::Key;
using Unit = extraction::Unit;
constexpr int seed_size = 8;

struct Side {
    int length = 0;
    int mismatches = 0;
};

// Verify all available repeat bases up to the physical read boundary. The
// nearest four bases fix the junction. A short hook has to be entirely exact.
inline Side side(const std::string& seq, const std::string& repeat, int cut, bool upstream) {
    const int repeat_size = static_cast<int>(repeat.size());
    const int length = std::min(repeat_size, upstream ? cut : static_cast<int>(seq.size()) - cut);
    if (length < 4) return {};
    Side result{length, 0};
    for (int d = 0; d < length; ++d) {
        const char base = seq[upstream ? cut-1-d : cut+d];
        if (base == repeat[upstream ? repeat_size-1-d : d]) continue;
        if (length < 17 || d < 4 || read_evidence::base(base) < 0 || ++result.mismatches > 1) return {};
    }
    return result;
}

inline bool geometric(const Side& left, const Side& right, int size) {
    if (!left.length || !right.length) return false;
    const auto shorter = std::min(left.length, right.length);
    const auto longer = std::max(left.length, right.length);
    if (shorter < 10) return longer == size && left.mismatches + right.mismatches == 0;
    return longer >= 16 && left.length + right.length >= 32;
}

// The seed stage records raw interior observations on physical fragments,
// even when the strict caller ultimately emits no spacer for this family.
inline bool supported(const Family& family, const std::string& spacer, int tip, std::uint64_t focal) {
    const auto canonical = read_linked::canonical(spacer);
    std::set<std::string> other_types;
    FragmentContexts nearby;
    const int radius = tip <= 6 ? 0 : tip < 10 ? 1 : 4;
    for (const auto& b : family.boundaries) {
        if (b.fragment == focal || !seeds::different_interiors(canonical, b.spacer)) continue;
        other_types.insert(b.spacer);
        if (std::abs(static_cast<int>(b.spacer.size()) - static_cast<int>(spacer.size())) <= radius)
            nearby[b.spacer].insert(b.fragment);
    }
    if (nearby.size() < 2 || (tip <= 6 && 4 * nearby.size() < 3 * other_types.size())) return false;
    // A repeated interior and its sequencing errors cannot provide two votes.
    return independent_contexts(nearby) >= 2;
}

struct Placement {
    std::size_t strand;
    int begin, end;
    Side left, right;
};

class Scanner {
    const std::vector<Family>& families_;
    std::vector<extraction::Strand> strands_;
    std::vector<std::pair<std::size_t,int>> seeds_;
    read_evidence::ExactMatcher matcher_;

public:
    explicit Scanner(const std::vector<Family>& families) : families_(families) {
        for (std::size_t i = 0; i < families.size(); ++i) {
            const auto& family = families[i];
            if (family.templates < 2 || family.contexts < 2 || family.repeat.size() < 20 ||
                family.repeat.size() > 55 || family.boundaries.size() < 2) continue;
            const auto& repeat = family.repeat;
            const auto rc = read_linked::reverse_complement(repeat);
            if (repeat == rc) continue;
            const int n = static_cast<int>(repeat.size());
            for (bool reverse : {false, true}) {
                const auto& seq = reverse ? rc : repeat;
                const auto id = strands_.size();
                strands_.push_back({i, reverse, seq});
                // Disjoint pairs of junction-facing seeds retain a proposal
                // when a single internal repeat substitution hits one seed.
                for (const int p : std::set<int>{0, seed_size, n-2*seed_size, n-seed_size}) {
                    matcher_.add(seq.substr(p,seed_size), seeds_.size());
                    seeds_.emplace_back(id,p);
                }
            }
        }
        matcher_.build();
    }

    bool empty() const { return strands_.empty(); }

    bool scan(const klibpp::KSeq& read, const klibpp::KSeq& first, const klibpp::KSeq* mate,
              std::map<Key,bool>& pending, Counts& counts) const {
        if (read.qual.empty()) return true;
        const int n = static_cast<int>(read.seq.size());
        if (read.seq.size() > read_recruitment::Config{}.maximum_read_length) return false;
        std::set<std::pair<std::size_t,int>> starts;
        std::size_t hits = 0;
        matcher_.scan(read.seq, [&](std::size_t id, std::size_t end) {
            if (++hits <= 4096)
                starts.emplace(seeds_[id].first, static_cast<int>(end)-seed_size-seeds_[id].second);
        });
        if (hits > 4096 || starts.size() > 256) return false;
        if (starts.empty()) return true;

        read_evidence::QualitySpans quality(read.qual, 0);
        std::set<std::tuple<std::size_t,int,int>> seen;
        std::vector<Placement> plausible;
        for (const auto& start : starts) {
            const auto& strand = strands_[start.first];
            const int size = static_cast<int>(strand.repeat.size());
            for (bool upstream : {true,false}) {
                const int cut = upstream ? start.second + size : start.second;
                if (cut < 4 || n - cut < 4) continue;
                const auto anchor = side(read.seq, strand.repeat, cut, upstream);
                if (!anchor.length) continue;
                for (int gap = seeds::min_spacer; gap <= seeds::max_spacer; ++gap) {
                    const int begin = upstream ? cut : cut-gap;
                    const int end = begin + gap;
                    if (begin < 4 || end + 4 > n || !seen.emplace(start.first,begin,end).second) continue;
                    if (seen.size() > 8192) return false;
                    const auto opposite = side(read.seq, strand.repeat, upstream ? end : begin, !upstream);
                    const auto left = upstream ? anchor : opposite;
                    const auto right = upstream ? opposite : anchor;
                    if (!geometric(left,right,size)) continue;
                    const int tip = std::min(left.length,right.length);
                    if (read.seq.find_first_not_of("ACGT", begin) < static_cast<std::size_t>(end)) continue;
                    // Trust aggregate Phred evidence in the actual spacer and
                    // adjoining two 8-base junction windows. A lone low Q
                    // deep inside a matching repeat need not discard the read.
                    if (!quality.passes(begin,end,tip < 10 ? 0.9 : 1.5) ||
                        !quality.passes(begin-4,begin+4,0.6) || !quality.passes(end-4,end+4,0.6)) {
                        ++counts["sparse_quality_rejected"];
                        continue;
                    }
                    plausible.push_back({start.first,begin,end,left,right});
                    if (plausible.size() > 256) return false;
                }
            }
        }
        if (plausible.empty()) return true;
        const auto fragment = seeds::seed_fragment(first,mate);
        std::vector<Placement> certified;
        for (const auto& p : plausible) {
            const auto& model = strands_[p.strand];
            const auto& family = families_[model.family];
            const auto s = read.seq.substr(p.begin,p.end-p.begin);
            if (!seeds::admissible_interior(model.repeat,s) ||
                !supported(family,s,std::min(p.left.length,p.right.length),fragment)) {
                ++counts["sparse_family_rejected"];
                continue;
            }
            certified.push_back(p);
        }
        for (const auto& p : certified) {
            // Compare all independently supported models. A different cut or
            // repeat family at this interval makes the assignment ambiguous.
            const auto ambiguous = std::any_of(certified.begin(),certified.end(),[&](const auto& q) {
                return &p != &q && p.begin < q.end && q.begin < p.end;
            });
            if (ambiguous) { ++counts["sparse_ambiguous_placements"]; continue; }
            const auto& model = strands_[p.strand];
            auto seq = read.seq.substr(p.begin,p.end-p.begin);
            if (model.reverse) seq = read_linked::reverse_complement(seq);
            if (extraction::tandem_explains(families_[model.family].repeat,seq)) continue;
            const int size = static_cast<int>(model.repeat.size());
            auto& full = pending[{model.family,seq}];
            full |= p.left.length == size && p.right.length == size;
            ++counts[std::min(p.left.length,p.right.length) < 10 ? "sparse_hook_placements" : "sparse_bilateral_placements"];
        }
        return true;
    }
};
} // namespace sparse_junctions

inline std::vector<Call> recruit_sparse(const std::vector<Family>& families,
        const std::vector<std::string>& paths, std::map<std::string,std::size_t>& counts) {
    using namespace sparse_junctions;
    Scanner scanner(families);
    if (paths.empty() || scanner.empty()) return {};
    std::map<Key,std::map<Unit,bool>> inventory;
    std::size_t memberships = 0;
    bool limited = false;
    seeds::each_template(paths,[&](const klibpp::KSeq& first, const klibpp::KSeq* second) {
        if (limited) return;
        std::map<Key,bool> pending;
        if (!scanner.scan(first,first,second,pending,counts) ||
            (second && !scanner.scan(*second,first,second,pending,counts))) {
            ++counts["sparse_limited_templates"];
            return;
        }
        if (pending.empty()) return;
        const auto unit = read_recruitment::detail::unit_identity(first.seq,second ? second->seq : "",bool(second));
        for (const auto& [key,full] : pending) {
            auto added = inventory[key].emplace(unit,full);
            memberships += added.second;
            added.first->second |= full;
        }
        if (inventory.size() > 65536 || memberships > 250000) {
            limited = true;
            inventory.clear();
        }
    });
    counts["sparse_global_resource_limit"] = limited;
    counts["sparse_candidate_types"] = inventory.size();
    std::vector<Call> result;
    for (const auto& [key,units] : inventory) {
        const auto full = std::count_if(units.begin(),units.end(),[](const auto& entry) { return entry.second; });
        result.push_back({families[key.first].repeat,key.second,units.size(),static_cast<std::size_t>(full)});
    }
    return result;
}

}} // namespace Zuck::read_first
#endif
