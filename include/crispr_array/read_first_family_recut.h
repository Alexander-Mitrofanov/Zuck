#ifndef ZUCK_READ_FIRST_FAMILY_RECUT_H
#define ZUCK_READ_FIRST_FAMILY_RECUT_H

#include "crispr_array/read_first.h"
#include <deque>

namespace Zuck { namespace read_first { namespace family_recut {

using Geometry = std::tuple<std::string, std::size_t, std::size_t>;  // Canonical graph core, left/right extension.

// With long repeat P+G+S and read spacer X, the exact short-model spacer must
// be S+X+P. A circular-unit match alone does not establish this placement.
inline bool extension(std::string graph, std::string graph_spacer, const Call& call,
                      Geometry& geometry, std::string& aligned) {
    if (graph.empty() || call.repeat.size() <= graph.size() ||
        graph_spacer.size() != call.spacer.size() + call.repeat.size() - graph.size()) return false;
    const auto reverse = read_linked::reverse_complement(graph);
    if (reverse == graph) return false;  // The core cannot orient the family comparisons.
    if (reverse < graph) {
        graph = reverse;
        graph_spacer = read_linked::reverse_complement(graph_spacer);
    }
    bool matched = false;
    std::size_t placements = 0;
    for (bool rc : {false, true}) {
        const auto repeat = rc ? read_linked::reverse_complement(call.repeat) : call.repeat;
        const auto spacer = rc ? read_linked::reverse_complement(call.spacer) : call.spacer;
        for (auto pos = repeat.find(graph); pos != std::string::npos; pos = repeat.find(graph, pos+1)) {
            ++placements;
            const auto right = repeat.size() - pos - graph.size();
            if (graph_spacer != repeat.substr(pos + graph.size()) + spacer + repeat.substr(0, pos)) continue;
            geometry = {graph, pos, right};
            aligned = repeat;
            matched = true;
        }
    }
    return matched && placements == 1;
}

inline bool dissimilar(const std::string& a, const std::string& b) {
    const auto reverse = read_linked::reverse_complement(b);
    return seeds::different_interiors(a, b) && seeds::different_interiors(a, reverse) &&
           4 * std::min(read_recruitment::detail::edit_distance(a, b),
                        read_recruitment::detail::edit_distance(a, reverse)) >= std::max(a.size(), b.size());
}

struct Branch {
    const Family* model;
    FragmentContexts interiors;
    std::set<std::uint64_t> fragments;

    explicit Branch(const Family& family) : model(&family) {
        // A seed witness may clip the opposite repeat at the read end. Only
        // two wholly observed exact copies certify an independent interior.
        for (const auto& b : family.boundaries) {
            if (b.full_repeat_pair) interiors[b.spacer].insert(b.fragment);
            fragments.insert(b.fragment);
        }
    }
};

// Require two divergent nonfocal interiors assignable to distinct fragments.
// A template shared with the competing branch or a focal-like interior cannot
// corroborate either branch. Alleles or shifted copies of the same interior
// also cannot establish both family identities.
inline bool other_interiors(const Branch& own, const Branch& rival, const std::string& focal,
                            const std::set<std::uint64_t>& focal_fragments) {
    std::set<std::uint64_t> excluded = rival.fragments;
    excluded.insert(focal_fragments.begin(), focal_fragments.end());
    for (const auto& entry : own.interiors)
        if (!seeds::different_interiors(entry.first, focal))
            excluded.insert(entry.second.begin(), entry.second.end());
    FragmentContexts qualified;
    for (const auto& entry : own.interiors) {
        if (!dissimilar(entry.first, focal)) continue;
        if (std::any_of(rival.interiors.begin(), rival.interiors.end(), [&](const auto& other) {
                return !seeds::different_interiors(entry.first, other.first);
            })) continue;
        for (auto fragment : entry.second)
            if (!excluded.count(fragment)) qualified[entry.first].insert(fragment);
    }
    for (auto a = qualified.begin(); a != qualified.end(); ++a)
        for (auto b = qualified.begin(); b != a; ++b)
            if (dissimilar(a->first, b->first) &&
                (a->second.size() > 1 || b->second.size() > 1 || *a->second.begin() != *b->second.begin())) return true;
    return false;
}

// The core and extension geometry are identical. Sequence divergence must
// therefore lie in the added terminal sequence, away from the shared core.
inline bool distinct_signatures(const std::string& a, const std::string& b) {
    std::size_t mismatches = 0;
    for (std::size_t i = 0; i < a.size(); ++i) mismatches += a[i] != b[i];
    return mismatches >= 2;
}

struct Cohort {
    std::map<std::string, Branch> branches;  // Key is the long repeat in graph orientation.
    bool crossed = false;
};

struct Decision {
    const char* reason = "no_read_alternative";
    std::size_t other_contexts = 0;  // Verified lower bound, excluding the focal spacer.
    std::string competing_repeat;
};

struct FocalEvidence {
    Geometry geometry;
    std::string aligned;
    const Family* family;
    std::set<read_recruitment::detail::Unit> units;
    std::set<std::uint64_t> seed_fragments;
    bool truncated = false;
};

class Index {
    std::map<Geometry, Cohort> cohorts_;
    std::map<const Call*, FocalEvidence> focal_;
    std::set<Geometry> incomplete_geometries_;
    bool limited_ = false;
    bool raw_limited_ = false;
    static constexpr std::size_t branch_limit = 32;
    static constexpr std::size_t query_limit = 4096, per_query_limit = 4096, identifier_limit = 262144;

public:
    void register_call(const std::string& graph, const std::string& spacer, const Call& call, const Family& family) {
        Geometry geometry;
        std::string aligned;
        if (!extension(graph, spacer, call, geometry, aligned)) return;
        if (focal_.size() >= query_limit && !focal_.count(&call)) { raw_limited_ = true; return; }
        focal_.emplace(&call, FocalEvidence{geometry,aligned,&family,{},{},false});
    }

    // The strict anchor extractor deliberately rejects overlapping alternative
    // cuts on the same read. Directly query the full physical sequence for
    // already eligible, uniquely nested, within-read candidates. Whole-unit
    // presence supplies fragment counts and cannot discriminate the two cuts.
    void check_focal_templates(const std::vector<std::string>& paths) {
        if (focal_.empty() || raw_limited_) return;
        struct Pattern { FocalEvidence* focal; std::size_t length; };
        std::vector<Pattern> patterns;
        read_evidence::ExactMatcher matcher;
        for (auto& entry : focal_) {
            const auto& call = *entry.first;
            const auto whole = call.repeat + call.spacer + call.repeat;
            const auto reverse = read_linked::reverse_complement(whole);
            const auto id = patterns.size();
            matcher.add(whole,id);
            if (reverse != whole) matcher.add(reverse,id);
            patterns.push_back({&entry.second,whole.size()});
        }
        matcher.build();
        std::size_t identifiers = 0;
        seeds::each_template(paths, [&](const klibpp::KSeq& a, const klibpp::KSeq* b) {
            if (raw_limited_) return;
            std::set<std::size_t> pending;
            const auto inspect = [&](const klibpp::KSeq& read) {
                if (read.qual.empty() || raw_limited_) return;
                if (read.seq.size() != read.qual.size()) throw std::runtime_error("Mismatched focal read sequence and quality lengths");
                if (read.seq.size() > read_recruitment::Config{}.maximum_read_length) { raw_limited_ = true; return; }
                std::size_t hits = 0;
                matcher.scan(read.seq, [&](std::size_t id, std::size_t end) {
                    if (raw_limited_) return;
                    if (++hits > 1024) { raw_limited_ = true; return; }
                    const auto& pattern = patterns[id];
                    if (!pattern.focal->truncated &&
                        std::all_of(read.qual.begin() + (end-pattern.length), read.qual.begin() + end,
                                    [](unsigned char q) { return q >= 33+20; })) pending.insert(id);
                });
            };
            inspect(a);
            if (b) inspect(*b);
            if (pending.empty() || raw_limited_) return;
            const auto unit = read_recruitment::detail::unit_identity(a.seq,b ? b->seq : "",b != nullptr);
            const auto fragment = seeds::seed_fragment(a,b);
            for (auto id : pending) {
                auto& e = *patterns[id].focal;
                if (e.units.size() >= per_query_limit && !e.units.count(unit)) { e.truncated = true; continue; }
                if (e.seed_fragments.size() >= per_query_limit && !e.seed_fragments.count(fragment)) {
                    e.truncated = true;
                    continue;
                }
                identifiers += e.units.insert(unit).second;
                identifiers += e.seed_fragments.insert(fragment).second;
                if (identifiers > identifier_limit) { raw_limited_ = true; return; }
            }
        });
        if (raw_limited_) return;
        for (const auto& entry : focal_) {
            const auto& candidate = *entry.first;
            const auto& e = entry.second;
            if (e.truncated) {
                // A truncated branch could conceal an additional competitor or mixed transition.
                incomplete_geometries_.insert(e.geometry);
                continue;
            }
            if (!supported(candidate,e)) continue;
            auto& branches = cohorts_[e.geometry].branches;
            // One extra entry marks the saturated geometry ineligible.
            if (branches.size() <= branch_limit) branches.try_emplace(e.aligned, *e.family);
        }
    }

    static bool supported(const Call& candidate, const FocalEvidence& e) {
        return (candidate.templates >= 2 && candidate.full_pairs >= 1) || e.units.size() >= 2;
    }

    void log_focal(std::ostream& stream, const Call* candidate) const {
        const auto item = focal_.find(candidate);
        if (item == focal_.end()) {
            stream << "\tNA\t" << (raw_limited_ ? "global_limit" : "not_queried") << "\tNA";
            return;
        }
        const auto& e = item->second;
        const char* status = raw_limited_ ? "global_limit" : e.truncated ? "per_query_limit" : "complete";
        const char* route = raw_limited_ || e.truncated ? "NA" :
                            candidate->templates >= 2 && candidate->full_pairs ? "strict_recruitment" :
                            e.units.size() >= 2 ? "direct_full_templates" : "insufficient";
        stream << '\t' << e.units.size() << '\t' << status << '\t' << route;
    }

    // Detect changes in the terminal family signature on *one physical read*.
    // Two consecutive different long alleles are compatible with a common
    // short repeat flanked by blockwise conserved spacer ends. Abstain for the
    // entire geometry when even one high-quality mixed transition is observed.
    void check_transitions(const std::vector<std::string>& paths) {
        struct Pattern { Cohort* cohort; const std::string* repeat; bool reverse; };
        std::vector<Pattern> patterns;
        read_evidence::ExactMatcher matcher;
        for (auto& entry : cohorts_) {
            auto& cohort = entry.second;
            if (cohort.branches.size() < 2 || cohort.branches.size() > branch_limit) continue;
            for (const auto& branch : cohort.branches) {
                for (bool reverse : {false, true}) {
                    matcher.add(reverse ? read_linked::reverse_complement(branch.first) : branch.first, patterns.size());
                    patterns.push_back({&cohort, &branch.first, reverse});
                }
            }
        }
        if (patterns.empty()) return;
        matcher.build();
        const auto scan = [&](const klibpp::KSeq& read) {
            if (read.qual.empty() || limited_) return;
            if (read.seq.size() > read_recruitment::Config{}.maximum_read_length) { limited_ = true; return; }
            std::deque<std::pair<std::size_t,std::size_t>> recent;  // Pattern ID and exclusive read end.
            read_evidence::QualitySpans quality(read.qual, 20);
            matcher.scan(read.seq, [&](std::size_t id, std::size_t end) {
                const auto& p = patterns[id];
                if (p.cohort->crossed || limited_) return;
                while (!recent.empty() && recent.front().second + seeds::max_repeat + seeds::max_spacer < end)
                    recent.pop_front();
                if (recent.size() >= 1024) { limited_ = true; return; }
                const auto begin = end - p.repeat->size();
                for (const auto& previous : recent) {
                    const auto& q = patterns[previous.first];
                    const auto last = previous.second;
                    if (p.cohort == q.cohort && p.reverse == q.reverse && p.repeat != q.repeat &&
                        last + seeds::min_spacer <= begin && begin <= last + seeds::max_spacer &&
                        quality.passes(last - q.repeat->size(), end, 1.0)) {
                        p.cohort->crossed = true;
                        break;
                    }
                }
                recent.emplace_back(id,end);
            });
        };
        seeds::each_template(paths, [&](const klibpp::KSeq& a, const klibpp::KSeq* b) {
            scan(a);
            if (b) scan(*b);
        });
    }

    Decision decide(const std::string& graph, const std::string& spacer, const Call& observed) const {
        Geometry geometry;
        std::string aligned;
        if (!extension(graph, spacer, observed, geometry, aligned)) return {"not_unique_nested_recut",0,{}};
        if (raw_limited_) return {"raw_template_global_limit",0,{}};
        const auto focal = focal_.find(&observed);
        if (focal == focal_.end()) return {"unregistered_read_family",0,{}};
        const auto& e = focal->second;
        if (e.truncated) return {"raw_template_query_limit",0,{}};
        if (incomplete_geometries_.count(geometry)) return {"raw_template_cohort_limit",0,{}};
        if (!supported(observed,e)) return {"insufficient_focal_fragment_evidence",0,{}};
        const auto found = cohorts_.find(geometry);
        if (found == cohorts_.end() || found->second.branches.size() < 2) return {"no_competing_family",0,{}};
        const auto& cohort = found->second;
        if (limited_ || cohort.branches.size() > branch_limit) return {"family_evidence_limit",0,{}};
        if (cohort.crossed) return {"mixed_family_transition",0,{}};
        const auto own = cohort.branches.find(aligned);
        if (own == cohort.branches.end()) return {"unregistered_read_family",0,{}};
        const auto target = read_linked::canonical(observed.spacer);
        for (const auto& rival : cohort.branches) {
            if (!distinct_signatures(aligned,rival.first)) continue;
            if (other_interiors(own->second,rival.second,target,e.seed_fragments) &&
                other_interiors(rival.second,own->second,target,e.seed_fragments))
                return {"family_recut",2,rival.second.model->repeat};
        }
        return {"no_independent_competing_families",0,{}};
    }
};

}}} // namespace Zuck::read_first::family_recut
#endif
