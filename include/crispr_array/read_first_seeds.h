#ifndef ZUCK_READ_FIRST_SEEDS_H
#define ZUCK_READ_FIRST_SEEDS_H

#include "crispr_array/read_linked_groups.h"
#include <algorithm>
#include <array>
#include <cstdint>
#include <map>
#include <set>
#include <string>
#include <tuple>
#include <vector>

namespace Zuck { namespace read_first {

// Two consecutive copies in one physical read, in canonical repeat orientation.
// Each array stores bases at start-2, start-1, end, end+1 (end exclusive).
// Missing/out-of-quality positions are N and never count as disagreements.
struct BoundaryObservation {
    std::string spacer;
    std::uint64_t fragment;
    std::array<char,4> first, second;
    bool full_repeat_pair = false;
};

struct Family {
    std::string repeat;
    std::size_t templates = 0, contexts = 0;
    std::vector<BoundaryObservation> boundaries;
};

namespace seeds {

constexpr int core_length = 15;
constexpr int min_repeat = 20, max_repeat = 55, min_spacer = 26, max_spacer = 50;
constexpr int radius = 120;
constexpr std::size_t max_cores = 16384, windows_per_core = 64;

inline std::uint64_t dna_hash(const std::string& dna) {
    std::uint64_t hash = 14695981039346656037ULL;
    for (unsigned char c : dna) { hash ^= c; hash *= 1099511628211ULL; }
    return hash;
}

// Use the same fragment identifier when excluding a focal read from seed
// corroboration as was attached to that read's original seed witnesses.
inline std::uint64_t seed_fragment(const klibpp::KSeq& a, const klibpp::KSeq* b) {
    auto first = dna_hash(a.seq), second = b ? dna_hash(b->seq) : 0;
    if (second < first) std::swap(first, second);
    return first * 1099511628211ULL ^ second;
}

// Reject low complexity and near-periodic sequence before registering a seed.
inline bool diverse(const std::string& dna) {
    std::array<bool, 64> words{};
    unsigned code = 0;
    std::size_t unique = 0;
    for (std::size_t i = 0; i < dna.size(); ++i) {
        const int symbol = read_evidence::base(dna[i]);
        if (symbol < 0) return false;
        code = ((code << 2) | static_cast<unsigned>(symbol)) & 63;
        if (i >= 2 && !words[code]) { words[code] = true; ++unique; }
    }
    if (unique < 7) return false;
    for (std::size_t period = 1; period <= std::min<std::size_t>(12, dna.size()/2); ++period) {
        std::size_t matches = 0;
        for (std::size_t i = period; i < dna.size(); ++i) matches += dna[i] == dna[i-period];
        if (10 * matches >= 9 * (dna.size() - period)) return false;
    }
    return true;
}

template<class Callback>
void each_template(const std::vector<std::string>& paths, Callback callback) {
    if (paths.empty() || paths.size() > 2) throw std::invalid_argument("Read-first discovery requires one or two read files");
    read_linked::Reader first(paths[0]);
    std::unique_ptr<read_linked::Reader> second;
    if (paths.size() == 2) second = std::make_unique<read_linked::Reader>(paths[1]);
    klibpp::KSeq a, b;
    while (first.next(a)) {
        if (second && (!second->next(b) || read_linked::read_name(a.name) != read_linked::read_name(b.name)))
            throw std::runtime_error("Read-first input pairs have different names or lengths");
        callback(a, second ? &b : nullptr);
    }
    if (second && second->next(b)) throw std::runtime_error("Read-first input pair counts differ");
}

class PairedSeeds {
    std::array<std::uint64_t, 65536> last_{};
    std::uint64_t origin_ = 1;
    std::vector<int> words_;
public:
    template<class Emit> void scan(const klibpp::KSeq& read, Emit emit) {
        const auto& seq = read.seq;
        if (seq.size() < 65 || read.qual.empty()) return;
        words_.assign(seq.size(), -1);
        unsigned word = 0, valid = 0;
        read_evidence::QualitySpans qualities(read.qual, 10);
        for (std::size_t pos = 0; pos < seq.size(); ++pos) {
            const int base = read_evidence::base(seq[pos]);
            if (base < 0) { valid = 0; word = 0; }
            else { ++valid; word = ((word << 2) | static_cast<unsigned>(base)) & 65535; }
            if (pos < 7) continue;
            const std::size_t j = pos - 7;
            if (valid >= 8) words_[j] = static_cast<int>(word);
            // Delay insertion until a legal second copy can start. An 18-base
            // exact repeat tract always contains a sampled eight-base word.
            if (j >= 46 && (j - 46) % 8 == 0 && words_[j-46] >= 0)
                last_[words_[j-46]] = origin_ + j - 46;
            if (words_[j] < 0 || last_[word] < origin_) continue;
            const std::size_t i = last_[word] - origin_, delta = j - i;
            if (delta > max_repeat + max_spacer) continue;
            std::size_t left = 0, right = 8;
            while (left < i && left < 8 && seq[i-left-1] == seq[j-left-1] &&
                   read_evidence::base(seq[i-left-1]) >= 0) ++left;
            if (left == 8) continue;  // Earlier sampled seed describes this tract.
            while (j + right < seq.size() && right + left <= 60 &&
                   seq[i+right] == seq[j+right] && read_evidence::base(seq[j+right]) >= 0) ++right;
            const auto len = left + right;
            if (len < 18 || len > 60 || delta < len + 20 || delta > len + 58) continue;
            if (!qualities.passes(i-left, i+right, 0.7) || !qualities.passes(j-left, j+right, 0.7)) continue;
            const std::string tract = seq.substr(i-left, len);
            if (!diverse(tract) || !diverse(seq.substr(i+right, delta-len))) continue;
            // Contained and overextended repeat observations can nominate the
            // same strand-independent internal word. A word only proposes a
            // pileup; the physical reads determine its final repeat boundaries.
            std::string core;
            for (std::size_t start = 1; start + core_length < len; ++start) {
                auto kmer = read_linked::canonical(tract.substr(start, core_length));
                if (core.empty() || kmer < core) core = std::move(kmer);
            }
            emit(core);
        }
        origin_ += seq.size() + 1;
    }
};

struct Window {
    std::string sequence;
    std::size_t anchor;
    std::uint64_t fragment;
    bool at_start, at_end;
};
using Windows = std::map<std::uint64_t, Window>;

inline void retain(Windows& windows, const klibpp::KSeq& read, std::size_t start, bool reverse,
                   std::uint64_t fragment) {
    const std::size_t left = start > radius ? start - radius : 0;
    const std::size_t right = std::min(read.seq.size(), start + core_length + radius);
    auto seq = read.seq.substr(left, right-left);
    for (std::size_t i = 0; i < seq.size(); ++i) {
        const char c = seq[i];
        const int b = read_evidence::base(c);
        seq[i] = b < 0 || read.qual[left+i] < 43 ? 'N' : (reverse ? "TGCA"[b] : c);
    }
    if (reverse) std::reverse(seq.begin(), seq.end());
    const std::size_t anchor = reverse ? right-start-core_length : start-left;
    if (seq.find('N', anchor) < anchor + core_length) return;
    // Deterministic bottom-k by sequence, anchor and template, independent of
    // read order. Identical template sequences cannot create independent proof.
    std::uint64_t key = (dna_hash(seq) ^ anchor) * 1099511628211ULL ^ fragment;
    if (windows.size() == windows_per_core && key >= windows.rbegin()->first) return;
    windows.emplace(key, Window{std::move(seq), anchor, fragment,
                    reverse ? right == read.seq.size() : left == 0,
                    reverse ? left == 0 : right == read.seq.size()});
    if (windows.size() > windows_per_core) windows.erase(std::prev(windows.end()));
}

inline char flank(const Window& w, int offset) {
    const auto p = static_cast<std::ptrdiff_t>(w.anchor) + offset;
    return p < 0 || static_cast<std::size_t>(p) >= w.sequence.size() ? 'N' : w.sequence[p];
}

inline bool contrast(const std::vector<std::pair<char, std::uint64_t>>& flanks) {
    for (std::size_t i = 0; i < flanks.size(); ++i)
        for (std::size_t j = i+1; j < flanks.size(); ++j)
            if (flanks[i].first != flanks[j].first && flanks[i].second != flanks[j].second)
                return true;
    return false;
}

inline bool different_interiors(const std::string& a, const std::string& b) {
    // Exclude SNP neighbours and small shifts/near-contained tracts when
    // deciding whether two reads witness distinct spacer contexts.
    const auto& small = a.size() <= b.size() ? a : b;
    const auto& large = a.size() <= b.size() ? b : a;
    for (std::size_t offset = 0; offset <= large.size()-small.size(); ++offset) {
        std::size_t differences = 0;
        for (std::size_t i = 0; i < small.size(); ++i) differences += small[i] != large[offset+i];
        if (differences <= 2) return false;
    }
    return true;
}

inline bool admissible_interior(const std::string& repeat, const std::string& spacer) {
    if (!diverse(spacer) || !different_interiors(repeat, spacer)) return false;
    // Reject spacers explained by another phase of a near-exact tandem repeat.
    for (std::size_t shift = 0; shift < repeat.size(); ++shift) {
        std::size_t differences = 0;
        for (std::size_t i = 0; i < spacer.size() && differences <= 2; ++i)
            differences += spacer[i] != repeat[(shift+i) % repeat.size()];
        if (differences <= 2) return false;
    }
    return true;
}

inline Family infer(const std::string& core, const Windows& windows) {
    if (windows.size() < 2) return {};
    auto conserved = [&](int offset) {
        std::array<unsigned, 4> bases{};
        unsigned total = 0;
        for (const auto& entry : windows) {
            const auto b = read_evidence::base(flank(entry.second, offset));
            if (b >= 0) { ++bases[b]; ++total; }
        }
        const auto best = std::max_element(bases.begin(), bases.end());
        return total >= 2 && 5 * *best >= 4 * total ? "ACGT"[best-bases.begin()] : 'N';
    };
    std::string left, right;
    for (int i = -1; i >= -max_repeat; --i) {
        const auto c = conserved(i);
        if (c == 'N') break;
        left.push_back(c);
    }
    for (int i = core_length; i <= core_length+max_repeat; ++i) {
        const auto c = conserved(i);
        if (c == 'N') break;
        right.push_back(c);
    }
    std::reverse(left.begin(), left.end());
    const std::string repeat = left + core + right;
    if (repeat.size() < min_repeat || repeat.size() > max_repeat || !diverse(repeat)) return {};
    const auto canonical = read_linked::canonical(repeat);
    Family family{canonical, 0, 0, {}};
    std::vector<std::pair<char, std::uint64_t>> front, back;
    std::map<std::string, std::set<std::uint64_t>> interiors;
    for (const auto& entry : windows) {
        const auto& w = entry.second;
        if (w.anchor < left.size()) continue;
        const auto start = w.anchor - left.size(), end = start + repeat.size();
        if (w.sequence.compare(start, repeat.size(), repeat) != 0) continue;
        if (start && w.sequence[start-1] != 'N') front.emplace_back(w.sequence[start-1], w.fragment);
        if (end < w.sequence.size() && w.sequence[end] != 'N') back.emplace_back(w.sequence[end], w.fragment);
        const auto witness = [&](std::size_t begin, std::size_t finish) {
            auto spacer = w.sequence.substr(begin, finish - begin);
            if (!admissible_interior(repeat, spacer)) return;
            interiors[spacer].insert(w.fragment);
            const auto adjacent = [&](int repeat_start) {
                const int s = repeat_start - static_cast<int>(w.anchor);
                const int e = s + static_cast<int>(repeat.size());
                return std::array<char,4>{{flank(w,s-2), flank(w,s-1), flank(w,e), flank(w,e+1)}};
            };
            auto first = adjacent(static_cast<int>(begin) - static_cast<int>(repeat.size()));
            auto second = adjacent(static_cast<int>(finish));
            if (canonical != repeat) {
                std::swap(first, second);
                for (auto* bases : {&first, &second}) {
                    std::reverse(bases->begin(), bases->end());
                    for (char& c : *bases) {
                        const auto b = read_evidence::base(c);
                        c = b < 0 ? 'N' : "TGCA"[b];
                    }
                }
            }
            const bool full_pair = begin >= repeat.size() && finish + repeat.size() <= w.sequence.size() &&
                w.sequence.compare(begin-repeat.size(), repeat.size(), repeat) == 0 &&
                w.sequence.compare(finish, repeat.size(), repeat) == 0;
            family.boundaries.push_back({read_linked::canonical(spacer), w.fragment, first, second, full_pair});
        };
        for (std::size_t gap = min_spacer; gap <= max_spacer; ++gap) {
            if (end + gap + 12 <= w.sequence.size()) {
                const auto anchor = end + gap, n = std::min(repeat.size(), w.sequence.size()-anchor);
                if ((n == repeat.size() || w.at_end) && w.sequence.compare(anchor, n, repeat, 0, n) == 0) {
                    witness(end, anchor);
                }
            }
            if (start >= gap + 12) {
                const auto anchor = start - gap, n = std::min(repeat.size(), anchor);
                if ((n == repeat.size() || w.at_start) &&
                    w.sequence.compare(anchor-n, n, repeat, repeat.size()-n, n) == 0) {
                    witness(anchor, start);
                }
            }
        }
    }
    if (!contrast(front) || !contrast(back)) return {};
    std::set<std::uint64_t> templates;
    bool independent = false;
    for (auto a = interiors.begin(); a != interiors.end(); ++a) {
        templates.insert(a->second.begin(), a->second.end());
        for (auto b = interiors.begin(); !independent && b != a; ++b)
            if (different_interiors(a->first, b->first) &&
                (a->second.size() > 1 || b->second.size() > 1 || *a->second.begin() != *b->second.begin()))
                independent = true;
    }
    if (!independent) return {};
    family.templates = templates.size();
    family.contexts = interiors.size();
    return family;
}

} // namespace seeds

// Map fully comparable flanks to a proposed extension. Keep unknown separate
// from measured contradictory DNA, including when a mate ends in a repeat.
inline int extension_state(const std::array<char,4>& bases, std::size_t start,
                           const std::string& extension) {
    bool complete = true;
    for (std::size_t i = 0; i < extension.size(); ++i) {
        const char observed = bases[start+i];
        if (observed == 'N') complete = false;
        else if (observed != extension[i]) return -1;
    }
    return complete ? 1 : 0;
}

using FragmentContexts = std::map<std::string, std::set<std::uint64_t>>;

// A deterministic lower bound on dissimilar spacer contexts assignable to
// distinct templates. Counting a repeated context on eight reads still adds one.
inline std::size_t independent_contexts(const FragmentContexts& observations) {
    std::set<std::uint64_t> used;
    std::vector<const std::string*> accepted;
    for (const auto& entry : observations) {
        if (std::any_of(accepted.begin(), accepted.end(), [&](const auto* previous) {
                return !seeds::different_interiors(entry.first, *previous);
            })) continue;
        const auto fragment = std::find_if(entry.second.begin(), entry.second.end(),
                                          [&](auto id) { return !used.count(id); });
        if (fragment == entry.second.end()) continue;
        used.insert(*fragment);
        accepted.push_back(&entry.first);
    }
    return accepted.size();
}

// Positive only for a single unambiguous, exactly nested one/two-base cut.
// The observations come from the shorter repeat's wholly internal seed core;
// the hypothesized added bases were never required to enter this pileup.
inline std::size_t shorter_cut_contexts(const Family& shorter, const std::string& longer,
                                        const std::string& excluded_spacer) {
    if (longer.size() <= shorter.repeat.size() || longer.size() > shorter.repeat.size()+2 ||
        shorter.boundaries.empty()) return 0;
    std::set<std::pair<std::string,std::string>> overhangs;
    for (const auto& strand : {longer, read_linked::reverse_complement(longer)}) {
        for (auto p = strand.find(shorter.repeat); p != std::string::npos;
             p = strand.find(shorter.repeat, p+1)) {
            overhangs.emplace(strand.substr(0,p), strand.substr(p+shorter.repeat.size()));
        }
    }
    if (overhangs.size() != 1) return 0;
    const auto& prefix = overhangs.begin()->first;
    const auto& suffix = overhangs.begin()->second;
    const auto target = read_linked::canonical(excluded_spacer);
    std::set<std::uint64_t> excluded;
    if (!target.empty()) {
        for (const auto& b : shorter.boundaries)
            if (!seeds::different_interiors(b.spacer, target)) excluded.insert(b.fragment);
    }
    FragmentContexts contradict, fully_extended, mixed;
    for (const auto& b : shorter.boundaries) {
        if (excluded.count(b.fragment)) continue;
        const int fp = extension_state(b.first,  2-prefix.size(), prefix);
        const int sp = extension_state(b.second, 2-prefix.size(), prefix);
        const int fs = extension_state(b.first,  2, suffix);
        const int ss = extension_state(b.second, 2, suffix);
        if ((prefix.empty() || fp < 0 || sp < 0) && (suffix.empty() || fs < 0 || ss < 0))
            contradict[b.spacer].insert(b.fragment);
        if (fp == 1 && sp == 1 && fs == 1 && ss == 1)
            fully_extended[b.spacer].insert(b.fragment);
        if (fp*sp == -1 || fs*ss == -1)
            mixed[b.spacer].insert(b.fragment);
    }
    const auto contrary = independent_contexts(contradict);
    if (contrary < 2) return 0;
    const auto positive = independent_contexts(fully_extended);
    const auto crossed = independent_contexts(mixed);
    // Preserve a phased extension cohort: consecutive copies with the long
    // ending on several dissimilar spacers and little evidence of mixing.
    if (positive >= 2 && (crossed < 2 || positive >= 2*crossed)) return 0;
    FragmentContexts recurrent;
    for (const auto& p : fully_extended)
        if (p.second.size() >= 2 && !contradict.count(p.first)) recurrent.insert(p);
    if (independent_contexts(recurrent) >= 2) return 0;
    return contrary;
}

inline std::vector<Family> discover(const std::vector<std::string>& paths,
                                    std::map<std::string,std::size_t>& counts) {
    seeds::PairedSeeds finder;
    std::set<std::string> cores;
    bool saturated = false;
    std::size_t pairs = 0;
    seeds::each_template(paths, [&](const klibpp::KSeq& first, const klibpp::KSeq* second) {
        auto add = [&](const std::string& core) {
            ++pairs;
            if (cores.size() >= seeds::max_cores && !cores.count(core)) { saturated = true; return; }
            cores.insert(core);
        };
        finder.scan(first, add);
        if (second) finder.scan(*second, add);
    });
    counts["read_first_seed_pairs"] = pairs;
    counts["read_first_seed_cores"] = cores.size();
    counts["read_first_seed_capacity_exceeded"] = saturated;
    if (cores.empty() || saturated) return {};
    std::vector<std::string> names(cores.begin(), cores.end());
    std::vector<seeds::Windows> windows(names.size());
    read_evidence::ExactMatcher matcher;
    for (std::size_t i = 0; i < names.size(); ++i) {
        matcher.add(names[i], 2*i);
        matcher.add(read_linked::reverse_complement(names[i]), 2*i+1);
    }
    matcher.build();
    seeds::each_template(paths, [&](const klibpp::KSeq& a, const klibpp::KSeq* b) {
        std::uint64_t fragment = 0;
        bool identified = false;
        auto inspect = [&](const klibpp::KSeq& r) {
            if (r.qual.empty()) return;
            matcher.scan(r.seq, [&](std::size_t id, std::size_t end) {
                if (!identified) {
                    fragment = seeds::seed_fragment(a,b);
                    identified = true;
                }
                seeds::retain(windows[id/2], r, end-seeds::core_length, id%2, fragment);
            });
        };
        inspect(a);
        if (b) inspect(*b);
    });
    std::map<std::string,Family> families;
    std::size_t capped = 0;
    for (std::size_t i = 0; i < names.size(); ++i) {
        capped += windows[i].size() == seeds::windows_per_core;
        const Family family = seeds::infer(names[i], windows[i]);
        if (family.repeat.empty()) continue;
        auto& best = families[family.repeat];
        if (std::tie(family.contexts,family.templates) > std::tie(best.contexts,best.templates))
            best = family;
    }
    counts["read_first_at_window_cap"] = capped;
    std::vector<const Family*> ordered;
    for (const auto& f : families) ordered.push_back(&f.second);
    std::sort(ordered.begin(), ordered.end(), [](const Family* a, const Family* b) {
        return a->repeat.size() != b->repeat.size() ? a->repeat.size() < b->repeat.size() : a->repeat < b->repeat;
    });
    std::set<std::string> suppressed;
    for (const auto* child : ordered) {
        std::set<std::string> owners;
        for (std::size_t trim = 1; trim <= 2; ++trim) {
            for (std::size_t front = 0; front <= trim; ++front) {
                const auto cut = read_linked::canonical(child->repeat.substr(front, child->repeat.size()-trim));
                auto it = families.find(cut);
                if (it == families.end() || suppressed.count(cut)) continue;
                const auto& parent = it->second;
                if (2*parent.contexts < 3*child->contexts || parent.contexts < child->contexts+2 ||
                    2*parent.templates < 3*child->templates || parent.templates < child->templates+2) continue;
                if (shorter_cut_contexts(parent, child->repeat, "") >= 2) owners.insert(cut);
            }
        }
        if (owners.size() == 1) suppressed.insert(child->repeat);
    }
    counts["read_first_suppressed_terminal_extensions"] = suppressed.size();
    std::vector<Family> result;
    for (auto& f : families)
        if (!suppressed.count(f.first)) result.push_back(std::move(f.second));
    counts["read_first_discovered_families"] = result.size();
    return result;
}

}} // namespace Zuck::read_first
#endif
