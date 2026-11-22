#ifndef ZUCK_BOUNDARY_CANDIDATES_H
#define ZUCK_BOUNDARY_CANDIDATES_H

#include "crispr_array/element_confidence.h"
#include "io/result_bundle.h"

namespace Zuck { namespace boundary_candidates {
using Pair = std::pair<std::string, std::string>;
namespace db = array_database;
namespace ec = element_confidence;

// Coordinates refer to the original oriented R+S circle: [a, |R|+b)
// is the new repeat and [|R|+b, |R|+|S|+a) is the new spacer.
// Transfer bases across BOTH junctions; never invent flanking bases.
inline Pair shifted(const Pair& pair, int a, int b) {
    const std::string core = pair.first + pair.second;
    const int n = static_cast<int>(core.size());
    const int nr = static_cast<int>(pair.first.size()) + b - a;
    if (n == 0 || nr <= 0 || nr >= n) return {};
    const int start = (a % n + n) % n;
    const auto rotated = core.substr(start) + core.substr(0, start);
    return {rotated.substr(0, nr), rotated.substr(nr)};
}

struct Assessment {
    double same_array_score = 0, anchored_score = 0;
    std::size_t other_spacers = 0, other_cores = 0;
    std::string source, anchor_source;
};

class Search {
    int radius_;
    std::vector<Pair> bases_;
    std::vector<std::string> cores_;
    std::vector<std::vector<std::size_t>> contexts_;
    struct Context {
        std::string sequence;
        std::set<std::string> units;
        std::set<std::size_t> offsets;
    };
    std::vector<Context> patterns_;
    std::size_t maximum_read_length_ = 0;
    bool has_reads_ = false;
    std::map<std::string, double> repeat_scores_, spacer_scores_;
    std::map<Pair, Assessment> assessments_;
    std::map<Pair, std::size_t> supported_;
    std::map<Pair, int> radii_;
    std::map<std::string, std::set<std::string>> repeat_cores_;
public:
    std::map<std::string, std::size_t> counts;
    explicit Search(int radius) : radius_(radius) {
        if (radius < 0 || radius > 3) throw std::runtime_error("Boundary radius must be in [0,3]");
    }
    ec::Support support(std::size_t context) const {
        const auto& p = patterns_.at(context);
        return {p.units.size(), p.offsets.size(), p.units.size() >= ec::SUPPORT_CAP,
            p.offsets.size() >= ec::SUPPORT_CAP,
            has_reads_ && maximum_read_length_ >= p.sequence.size()};
    }
    template<class Visitor> void each(Visitor visit) const {
        for (std::size_t i = 0; i < bases_.size(); ++i)
            for (int a = -radius_; a <= radius_; ++a)
                for (int b = -radius_; b <= radius_; ++b) {
                    auto pair = shifted(bases_[i], a, b);
                    visit(i, a, b, pair, contexts_[i][a + radius_]);
                }
    }
    static bool valid(const Pair& p) {
        // The complete lattice is written, including out-of-range cuts.
        // Discovery limits apply only to promotion, not enumeration.
        return p.first.size() >= db::MIN_REPEAT && p.first.size() <= db::MAX_REPEAT &&
               p.second.size() >= db::MIN_SPACER && p.second.size() <= db::MAX_SPACER;
    }
    void prepare(const std::set<Pair>& bases, const std::vector<std::string>& files,
                 const std::vector<db::Record>& records) {
        if (files.size() > 2) throw std::runtime_error("Boundary evidence accepts at most two read files");
        bases_.assign(bases.begin(), bases.end());
        std::map<std::string, std::size_t> ids;
        for (const auto& pair : bases_) {
            if (!db::dna(pair.first) || !db::dna(pair.second) || pair.first.size() < 20 || pair.second.size() < 20)
                throw std::runtime_error("Invalid boundary seed");
            cores_.push_back(ec::circular_identity(pair.first + pair.second));
            contexts_.emplace_back();
            for (int a = -radius_; a <= radius_; ++a) {
                const auto p = shifted(pair, a, a);
                const auto span = db::canonical(p.first + p.second + p.first.substr(0, 16));
                const auto found = ids.emplace(span, patterns_.size());
                if (found.second) patterns_.push_back({span, {}, {}});
                contexts_.back().push_back(found.first->second);
            }
        }
        counts["boundary_seed_pairs"] = bases_.size();
        counts["boundary_combinations"] = bases_.size() * (2 * radius_ + 1) * (2 * radius_ + 1);
        counts["boundary_context_patterns"] = patterns_.size();
        ids.clear();
        has_reads_ = !files.empty();
        if (has_reads_) {
            read_evidence::ExactMatcher matcher;
            for (std::size_t i = 0; i < patterns_.size(); ++i) {
                matcher.add(patterns_[i].sequence, i);
                const auto rc = db::reverse_complement(patterns_[i].sequence);
                if (rc != patterns_[i].sequence) matcher.add(rc, i);
            }
            matcher.build();
            auto scan = [&](const klibpp::KSeq& r, std::map<std::size_t, std::set<std::size_t>>& hits) {
                ++counts["boundary_read_records"];
                if (r.qual.empty()) return;
                maximum_read_length_ = std::max(maximum_read_length_, r.seq.size());
                read_evidence::QualitySpans quality(r.qual, database_reads::MINIMUM_QUALITY);
                matcher.scan(r.seq, [&](std::size_t id, std::size_t end) {
                    const auto& p = patterns_[id];
                    if (p.units.size() >= ec::SUPPORT_CAP && p.offsets.size() >= ec::SUPPORT_CAP) return;
                    const auto start = end - p.sequence.size();
                    if (quality.passes(start, end, database_reads::MAXIMUM_EXPECTED_ERRORS)) hits[id].insert(start);
                });
            };
            read_linked::Reader first(files[0]);
            std::unique_ptr<read_linked::Reader> second;
            if (files.size() == 2) second = std::make_unique<read_linked::Reader>(files[1]);
            klibpp::KSeq a, b;
            while (first.next(a)) {
                if (second && (!second->next(b) || read_linked::read_name(a.name) != read_linked::read_name(b.name)))
                    throw std::runtime_error("Boundary read pairs have different names or lengths");
                std::map<std::size_t, std::set<std::size_t>> hits;
                scan(a, hits); if (second) scan(b, hits);
                if (hits.empty()) continue;
                const auto unit = a.seq + '\0' + (second ? b.seq : "");
                for (const auto& hit : hits) {
                    auto& p = patterns_[hit.first];
                    if (p.units.size() < ec::SUPPORT_CAP) p.units.insert(unit);
                    for (auto offset : hit.second) {
                        if (p.offsets.size() >= ec::SUPPORT_CAP) break;
                        p.offsets.insert(offset);
                    }
                }
            }
            if (second && second->next(b)) throw std::runtime_error("Boundary read files have different record counts");
        }
        each([&](std::size_t i, int a, int b, const Pair& p, std::size_t context) {
            if (!valid(p) || !support(context).accepted()) return;
            supported_.emplace(p, context);
            const auto radius = std::max(std::abs(a),std::abs(b));
            const auto inserted = radii_.emplace(p,radius);
            if (!inserted.second) inserted.first->second = std::min(inserted.first->second,radius);
            repeat_cores_[db::canonical(p.first)].insert(cores_[i]);
        });
        counts["boundary_supported_pairs"] = supported_.size();
        std::map<std::string, std::vector<Pair>> by_repeat;
        for (const auto& item : supported_) by_repeat[db::canonical(item.first.first)].push_back(item.first);
        std::map<Pair, std::string> core_for;
        for (const auto& item : supported_) core_for[item.first] = ec::circular_identity(item.first.first + item.first.second);
        for (const auto& record : records) {
            const auto minimum = record.bona_fide ? db::MIN_BONA_FIDE_SPACERS : db::MIN_POSSIBLE_SPACERS;
            if (record.confidence < (record.bona_fide ? db::MIN_BONA_FIDE_CONFIDENCE : db::MIN_POSSIBLE_CONFIDENCE)) continue;
            for (const auto& r : record.repeats) repeat_scores_[r] = std::max(repeat_scores_[r], record.weight());
            for (const auto& s : record.spacers) spacer_scores_[s] = std::max(spacer_scores_[s], record.weight());
            std::vector<Pair> all;
            std::set<std::pair<std::string, std::string>> anchors;
            for (const auto& repeat : record.repeats) {
                const auto found = by_repeat.find(repeat);
                if (found == by_repeat.end()) continue;
                for (const auto& p : found->second) {
                    all.push_back(p);
                    if (record.spacers.count(db::canonical(p.second))) {
                        anchors.emplace(db::canonical(p.second), core_for.at(p));
                    }
                }
            }
            for (const auto& p : all) {
                const bool is_known = record.spacers.count(db::canonical(p.second));
                std::set<std::string> other_spacers, other_cores;
                for (const auto& anchor : anchors)
                    if (anchor.first != db::canonical(p.second) && anchor.second != core_for.at(p)) {
                        other_spacers.insert(anchor.first); other_cores.insert(anchor.second);
                    }
                const auto others = std::min(other_spacers.size(), other_cores.size());
                if (others < minimum - (is_known ? 1 : 0)) continue;
                const double score = record.weight() * std::min<std::size_t>(others, 3) / 3.0;
                auto& result = assessments_[p];
                if (is_known && (score > result.same_array_score ||
                    (score == result.same_array_score && record.id < result.source))) {
                    result.same_array_score = score; result.source = record.id;
                    result.other_spacers = other_spacers.size(); result.other_cores = other_cores.size();
                }
                if (score > result.anchored_score || (score == result.anchored_score && record.id < result.anchor_source)) {
                    result.anchored_score = score; result.anchor_source = record.id;
                }
            }
        }
    }
    Assessment assessment(const Pair& p) const {
        const auto found = assessments_.find(p);
        return found == assessments_.end() ? Assessment{} : found->second;
    }
    double repeat_score(const Pair& p) const {
        const auto it = repeat_scores_.find(db::canonical(p.first));
        return it == repeat_scores_.end() ? 0 : it->second;
    }
    double spacer_score(const Pair& p) const {
        const auto it = spacer_scores_.find(db::canonical(p.second));
        return it == spacer_scores_.end() ? 0 : it->second;
    }
    std::set<Pair> select(const std::string& policy) const {
        std::set<Pair> result;
        if (policy == "off") return result;
        if (policy == "catalogue") {
            // Retain unique winners at every nested radius. Enlarging the
            // candidate search cannot erase a previously supported local cut.
            // 0.5 requires >=2 other anchors in a Bona-fide source (or maximal
            // eligible Possible evidence); the focal spacer/core never votes.
            for (int radius=std::min(1,radius_); radius<=radius_; ++radius) {
                std::map<std::string,std::vector<Pair>> windows;
                for (const auto& item : supported_)
                    if (radii_.at(item.first)<=radius && assessment(item.first).same_array_score>=0.5)
                        windows[ec::circular_identity(item.first.first+item.first.second)].push_back(item.first);
                for (const auto& window : windows) {
                    double best=0;
                    std::set<Pair> top, identities;
                    for (const auto& p : window.second) {
                        const double score=assessment(p).same_array_score;
                        if (score>best) {best=score;top.clear();identities.clear();}
                        if (score==best) {top.insert(p);identities.emplace(db::canonical(p.first),db::canonical(p.second));}
                    }
                    if (identities.size()==1) result.insert(top.begin(),top.end());
                }
            }
            return result;
        }
        const bool anchored = policy == "anchored" || policy == "anchored-best";
        const bool best_only = policy == "catalogue-best" || policy == "anchored-best";
        std::map<std::string, std::vector<Pair>> groups;
        for (const auto& item : supported_) {
            const auto e = assessment(item.first);
            if ((anchored ? e.anchored_score : e.same_array_score) <= 0) continue;
            if (!best_only) result.insert(item.first);
            else groups[ec::circular_identity(item.first.first + item.first.second)].push_back(item.first);
        }
        for (const auto& group : groups) {
            std::pair<bool,double> best{false,0};
            std::set<Pair> top, identities;
            for (const auto& p : group.second) {
                const auto e = assessment(p);
                const auto rank = std::make_pair(e.same_array_score > 0, anchored ? e.anchored_score : e.same_array_score);
                if (rank > best) { best=rank; top.clear(); identities.clear(); }
                if (rank == best) {
                    top.insert(p); identities.emplace(db::canonical(p.first), db::canonical(p.second));
                }
            }
            // Best policies abstain on equal-scoring distinct decompositions.
            if (identities.size()==1) result.insert(top.begin(), top.end());
        }
        return result;
    }
    template<class Arrays, class Lengths>
    void append(Arrays& arrays, const Lengths& lengths, const std::string& policy,
                const std::filesystem::path& path, std::map<std::string,std::size_t>& output_counts) const {
        output_counts.insert(counts.begin(), counts.end());
        output_counts["boundary_previous_groups"] = arrays.size();
        std::set<Pair> existing;
        for (const auto& array : arrays) existing.insert(array.entries.begin(),array.entries.end());
        std::map<std::string,std::vector<Pair>> groups;
        for (const auto& p : select(policy))
            if (lengths.pair_allowed(p.first,p.second)) groups[db::canonical(p.first)].push_back(p);
        std::ofstream out(path);
        if (!out) throw std::runtime_error("Cannot write boundary evidence: " + path.string());
        out << "array_id\tentry_index\trepeat_variant\tspacer\tdatabase_source\tanchor_source"
               "\tsame_array_score\tanchored_score\tcontext_sequence_units\tcontext_start_offsets\tnew_oriented_pair\n";
        out << std::setprecision(17);
        for (const auto& group : groups) {
            std::set<std::string> spacers;
            bool adds = false;
            for (const auto& p : group.second) { spacers.insert(db::canonical(p.second)); adds |= !existing.count(p); }
            if (!adds || spacers.size()<2 || !lengths.repeat_allowed(group.first.size())) continue;
            arrays.push_back({group.first,group.second});
            ++output_counts["boundary_groups_added"];
            output_counts["boundary_rows_added"] += group.second.size();
            for (std::size_t i=0; i<group.second.size(); ++i) {
                const auto& p=group.second[i]; const auto e=assessment(p); const auto s=support(supported_.at(p));
                out << "Array_" << arrays.size() << '\t' << i+1 << '\t' << p.first << '\t' << p.second
                    << '\t' << e.source << '\t' << e.anchor_source << '\t' << e.same_array_score << '\t' << e.anchored_score
                    << '\t' << s.sequence_units << '\t' << s.start_offsets << '\t' << !existing.count(p) << '\n';
                if (!existing.count(p)) ++output_counts["boundary_new_oriented_pairs"];
            }
        }
        result_bundle::close_checked(out,path);
    }
    void write(const std::filesystem::path& path) const {
        std::ofstream out(path);
        if (!out) throw std::runtime_error("Cannot write boundary candidates: " + path.string());
        out << "seed_id\trepeat_start_shift\trepeat_end_shift\trepeat_variant\tspacer\twithin_discovery_lengths"
               "\tcontext_sequence_units\tcontext_start_offsets\tcontext_observable\tcontext_score"
               "\trepeat_database_score\tspacer_database_score\tsame_array_score\tanchored_score"
               "\tother_spacers\tother_cores\tdatabase_source\tanchor_source\trepeat_context_cores\tcircular_core"
               "\tcontext_units_saturated\tcontext_offsets_saturated\n";
        out << std::setprecision(17);
        each([&](std::size_t i, int a, int b, const Pair& p, std::size_t context) {
            const auto s = support(context); const auto e = assessment(p);
            const auto rc = repeat_cores_.find(db::canonical(p.first));
            out << i + 1 << '\t' << a << '\t' << b << '\t' << p.first << '\t' << p.second
                << '\t' << valid(p) << '\t' << s.sequence_units << '\t' << s.start_offsets
                << '\t' << s.observable << '\t';
            ec::write_score(out, s.index());
            out << '\t' << repeat_score(p) << '\t' << spacer_score(p)
                << '\t' << e.same_array_score << '\t' << e.anchored_score
                << '\t' << e.other_spacers << '\t' << e.other_cores
                << '\t' << e.source << '\t' << e.anchor_source
                << '\t' << (rc == repeat_cores_.end() ? 0 : rc->second.size())
                << '\t' << cores_[i] << '\t' << s.units_saturated << '\t' << s.offsets_saturated << '\n';
        });
        result_bundle::close_checked(out, path);
    }
};
}} // namespace Zuck::boundary_candidates
#endif
