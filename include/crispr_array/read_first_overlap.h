#ifndef ZUCK_READ_FIRST_OVERLAP_H
#define ZUCK_READ_FIRST_OVERLAP_H

#include "crispr_array/read_first.h"

namespace Zuck { namespace read_first {
namespace overlap {

using Unit = read_recruitment::detail::Unit;
constexpr std::size_t minimum_overlap = 18;
constexpr std::size_t maximum_tip = 9;
constexpr std::size_t max_flank = read_recruitment::MAXIMUM_SPACER + maximum_tip;
constexpr std::size_t max_per_family = 4096, max_comparisons = 250000;

// Flank is in family orientation, immediately downstream (left=true) or
// upstream (left=false) of a complete repeat, extending to the actual read end.
using ContextKey = std::tuple<bool,std::string,Unit>;
using Bank = std::map<ContextKey,std::uint64_t>;
struct Context {
    bool left;
    std::string bases;
    Unit unit;
    std::uint64_t high_quality;
};

inline bool matches(const Context& c, const std::string& s, const std::string& repeat) {
    const auto length = std::min(s.size(), c.bases.size());
    const auto tip = c.bases.size() - length;
    if (tip > maximum_tip) return false;
    if (c.left)
        return c.bases.compare(0,length,s,0,length) == 0 &&
               c.bases.compare(length,tip,repeat,0,tip) == 0;
    return c.bases.compare(tip,length,s,s.size()-length,length) == 0 &&
           c.bases.compare(0,tip,repeat,repeat.size()-tip,tip) == 0;
}

inline bool corroborates(const Context& left, const Context& right, std::size_t n) {
    if (left.unit == right.unit) return false;
    const auto a = left.bases.size(), b = right.bases.size();
    // With <4 distal repeat bases neither read by itself locates that weak
    // boundary. Require the whole identical spacer in BOTH independent reads.
    if (!(a >= n && b >= n) && a < n+4 && b < n+4) return false;
    const auto end = std::min(a,n), start = n - std::min(b,n);
    if (start + minimum_overlap > end) return false;
    const auto right_offset = b - (n - start);
    const auto jointly_high = (left.high_quality >> start) & (right.high_quality >> right_offset);
    std::size_t run = 0;
    for (std::size_t i = 0; i < end-start; ++i) {
        run = ((jointly_high >> i) & 1U) ? run+1 : 0;
        if (run >= minimum_overlap) return true;
    }
    return false;
}

struct Hit { std::size_t owner, start, end; bool reverse; };
struct Pending { std::size_t family; bool left; std::string bases; std::uint64_t quality; };

inline bool collect(const klibpp::KSeq& read, const std::vector<Family>& families,
                    const std::vector<bool>& active, const read_evidence::ExactMatcher& matcher,
                    std::vector<Pending>& pending) {
    if (read.qual.empty()) return true;
    if (read.seq.size() > read_recruitment::Config{}.maximum_read_length) return false;
    std::vector<Hit> hits;
    matcher.scan(read.seq, [&](std::size_t id, std::size_t end) {
        if (hits.size() <= 128)
            hits.push_back({id/2,end-families[id/2].repeat.size(),end,bool(id%2)});
    });
    if (hits.size() > 128) return false;
    read_evidence::QualitySpans quality(read.qual,10);
    for (const auto& h : hits) {
        if (!active[h.owner]) continue;
        if (std::any_of(hits.begin(),hits.end(),[&](const auto& other) {
            return &h != &other && h.start < other.end && other.start < h.end;
        })) continue;
        for (bool downstream : {false,true}) {
            const auto begin = downstream ? h.end : 0;
            const auto end = downstream ? read.seq.size() : h.start;
            const auto length = end-begin;
            if (length < minimum_overlap || length > max_flank) continue;
            if (!quality.passes(downstream ? h.start : 0, downstream ? end : h.end,1.0)) continue;
            auto bases = read.seq.substr(begin,length);
            if (!read_recruitment::detail::dna(bases)) continue;
            if (h.reverse) bases = read_linked::reverse_complement(bases);
            std::uint64_t mask = 0;
            for (std::size_t i = 0; i < length; ++i)
                if (read.qual[h.reverse ? end-1-i : begin+i] >= 53) mask |= (std::uint64_t{1} << i);
            pending.push_back({h.owner,downstream != h.reverse,std::move(bases),mask});
        }
    }
    return true;
}

using Edge = std::pair<std::size_t,std::size_t>;
struct Candidate {
    std::vector<Edge> support;
};

inline bool has_close_allele(const std::string& spacer, const std::set<std::string>& pool) {
    if (pool.empty()) return false;
    for (std::size_t i = 0; i <= spacer.size(); ++i) {
        if (i < spacer.size() && pool.count(spacer.substr(0,i) + spacer.substr(i+1))) return true;
        for (char base : {'A','C','G','T'}) {
            if (pool.count(spacer.substr(0,i) + base + spacer.substr(i))) return true;
            if (i < spacer.size() && base != spacer[i] &&
                pool.count(spacer.substr(0,i) + base + spacer.substr(i+1))) return true;
        }
    }
    return false;
}

inline std::vector<Call> family_calls(const std::string& repeat, const Bank& bank,
                                      const std::set<std::string>& known, std::map<std::string,std::size_t>& counts) {
    std::vector<Context> contexts;
    std::map<std::string,std::vector<std::size_t>> left, right;
    std::map<std::string,Candidate> candidates;
    std::vector<bool> assigned_known;
    for (const auto& entry : bank) {
        const auto& key = entry.first;
        contexts.push_back({std::get<0>(key),std::get<1>(key),std::get<2>(key),entry.second});
        const auto& c = contexts.back();
        auto word = c.left ? c.bases.substr(0,minimum_overlap) : c.bases.substr(c.bases.size()-minimum_overlap);
        (c.left ? left : right)[word].push_back(contexts.size()-1);
        assigned_known.push_back(false);
        for (std::size_t tip = 0; tip <= maximum_tip; ++tip) {
            const auto n = c.bases.size()-tip;
            if (n < read_recruitment::MINIMUM_SPACER || n > read_recruitment::MAXIMUM_SPACER) continue;
            auto spacer = c.left ? c.bases.substr(0,n) : c.bases.substr(tip);
            if (matches(c,spacer,repeat) && !known.count(spacer)) candidates.try_emplace(std::move(spacer));
        }
    }
    // A shared end also belonging to a known spacer cannot establish a new allele.
    // Each whole S still originates in a single observed flank.
    for (const auto& s : known) {
        if (s.size() < minimum_overlap) continue;
        for (bool is_left : {false,true}) {
            const auto& index = is_left ? left : right;
            const auto word = is_left ? s.substr(0,minimum_overlap) : s.substr(s.size()-minimum_overlap);
            const auto found = index.find(word);
            if (found != index.end())
                for (auto i : found->second) assigned_known[i] = assigned_known[i] || matches(contexts[i],s,repeat);
        }
    }
    std::vector<std::size_t> owners(contexts.size());
    std::size_t comparisons = 0;
    for (auto& [s, candidate] : candidates) {
        const auto a = left.find(s.substr(0,minimum_overlap));
        const auto b = right.find(s.substr(s.size()-minimum_overlap));
        if (a == left.end() || b == right.end()) continue;
        std::set<std::size_t> retained;
        for (auto i : a->second) {
            if (assigned_known[i] || !matches(contexts[i],s,repeat)) continue;
            for (auto j : b->second) {
                if (++comparisons > max_comparisons) {
                    ++counts["read_first_overlap_pairing_limited_families"];
                    return {};
                }
                if (assigned_known[j] || !matches(contexts[j],s,repeat) ||
                    !corroborates(contexts[i],contexts[j],s.size())) continue;
                candidate.support.emplace_back(i,j);
                retained.insert(i);
                retained.insert(j);
            }
        }
        for (auto i : retained) ++owners[i];
    }
    std::set<std::string> corroborated;
    for (const auto& [s,candidate] : candidates)
        if (!candidate.support.empty()) corroborated.insert(s);
    std::vector<Call> calls;
    for (const auto& [s,candidate] : candidates) {
        if (candidate.support.empty()) continue;
        std::set<Unit> units, complete;
        for (const auto& edge : candidate.support) {
            if (owners[edge.first] != 1 || owners[edge.second] != 1) continue;
            for (auto i : {edge.first,edge.second}) {
                units.insert(contexts[i].unit);
                if (contexts[i].bases.size() >= s.size()) complete.insert(contexts[i].unit);
            }
        }
        if (units.size() < 2) { ++counts["read_first_overlap_ambiguous_paths"]; continue; }
        // If the only change could be an error in the single full sequence,
        // require an independent complete observation of that exact allele.
        if (complete.size() < 2 && (has_close_allele(s,known) || has_close_allele(s,corroborated))) {
            ++counts["read_first_overlap_uncorroborated_alleles"];
            continue;
        }
        if (!seeds::admissible_interior(repeat,s) || extraction::tandem_explains(repeat,s)) {
            ++counts["read_first_overlap_tandem_or_low_complexity"];
            continue;
        }
        calls.push_back({repeat,s,units.size(),0});
        ++counts[complete.size() >= 2 ? "read_first_overlap_paths_with_two_intact_templates" :
                                       "read_first_overlap_paths_with_one_intact_template"];
    }
    return calls;
}

} // namespace overlap

inline std::vector<Call> recruit_overlap(const std::vector<Family>& families,
                                         const std::vector<std::string>& paths,
                                         const std::vector<Call>& existing,
                                         std::map<std::string,std::size_t>& counts) {
    using namespace overlap;
    counts["read_first_overlap_novel_spacers"] = 0;
    if (families.empty() || paths.empty() || existing.empty()) return {};
    std::map<std::string,std::set<std::string>> known;
    std::set<std::string> old_spacers, old_units;
    for (const auto& c : existing) {
        auto repeat = read_linked::canonical(c.repeat);
        if (repeat.empty() || !read_recruitment::detail::dna(c.spacer)) continue;
        known[repeat].insert(c.repeat == repeat ? c.spacer : read_linked::reverse_complement(c.spacer));
        old_spacers.insert(read_linked::canonical(c.spacer));
        old_units.insert(read_recruitment::detail::core({c.repeat,c.spacer}));
    }
    std::vector<bool> active(families.size());
    read_evidence::ExactMatcher matcher;
    for (std::size_t i = 0; i < families.size(); ++i) {
        const auto& family = families[i];
        const auto& r = family.repeat;
        const auto reverse = read_linked::reverse_complement(r);
        if (r.size() < 20 || reverse.empty() || reverse == r) continue;
        matcher.add(r,2*i);
        matcher.add(reverse,2*i+1);
        if (family.contexts < 4 || family.templates < 5) continue;
        const auto& spacers = known[r];  // Discovered families are canonically oriented.
        for (auto b = spacers.begin(); b != spacers.end() && !active[i]; ++b)
            for (auto a = spacers.begin(); a != b && !active[i]; ++a)
                active[i] = seeds::different_interiors(*a,*b);
    }
    counts["read_first_overlap_eligible_families"] = std::count(active.begin(),active.end(),true);
    if (!counts["read_first_overlap_eligible_families"]) return {};
    matcher.build();
    std::vector<Bank> banks(families.size());
    std::set<std::size_t> limited;
    constexpr std::size_t max_total = 250000;
    std::size_t stored = 0;
    bool global_limit = false;
    seeds::each_template(paths,[&](const klibpp::KSeq& a,const klibpp::KSeq* b) {
        if (global_limit) return;
        std::vector<Pending> pending;
        if (!collect(a,families,active,matcher,pending) || (b && !collect(*b,families,active,matcher,pending))) {
            ++counts["read_first_overlap_limited_templates"];
            return;
        }
        if (pending.empty()) return;
        const auto id = read_recruitment::detail::unit_identity(a.seq,b ? b->seq : "",bool(b));
        for (auto& item : pending) {
            if (limited.count(item.family)) continue;
            auto& bank = banks[item.family];
            auto inserted = bank.emplace(ContextKey{item.left,std::move(item.bases),id},item.quality);
            if (!inserted.second) {
                // Retain one observed quality mask; combining bases from
                // separate copies can invent a jointly high-quality interval.
                auto& old = inserted.first->second;
                const auto rank = [](std::uint64_t mask) {
                    unsigned n = 0;
                    for (auto bits = mask; bits; bits &= bits-1) ++n;
                    return std::make_pair(n,mask);
                };
                if (rank(item.quality) > rank(old)) old = item.quality;
            } else if (++stored > max_total) {
                global_limit = true;
                return;
            }
            if (bank.size() > max_per_family) {
                limited.insert(item.family);
                bank.clear();
            }
        }
    });
    counts["read_first_overlap_context_limited_families"] = limited.size();
    counts["read_first_overlap_global_resource_limit"] = global_limit;
    counts["read_first_overlap_distinct_flanks_seen_before_limits"] = stored;
    if (global_limit) return {};
    std::vector<Call> candidates, result;
    std::map<std::string,std::set<std::string>> cuts;
    std::map<std::string,std::set<std::string>> owners;
    for (std::size_t i = 0; i < families.size(); ++i) {
        if (!active[i] || limited.count(i)) continue;
        for (auto& call : family_calls(families[i].repeat,banks[i],known[families[i].repeat],counts)) {
            const auto spacer = read_linked::canonical(call.spacer);
            const auto unit = read_recruitment::detail::core({call.repeat,call.spacer});
            cuts[unit].insert(spacer);
            owners[spacer].insert(call.repeat);
            candidates.push_back(std::move(call));
        }
    }
    for (const auto& call : candidates) {
        const auto spacer = read_linked::canonical(call.spacer);
        const auto unit = read_recruitment::detail::core({call.repeat,call.spacer});
        if (old_spacers.count(spacer) || old_units.count(unit) || cuts[unit].size() != 1 || owners[spacer].size() != 1)
            continue;
        result.push_back(call);
    }
    counts["read_first_overlap_novel_spacers"] = result.size();
    return result;
}

}} // namespace Zuck::read_first
#endif
