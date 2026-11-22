#ifndef ZUCK_READ_FIRST_VARIANTS_H
#define ZUCK_READ_FIRST_VARIANTS_H

#include "crispr_array/read_first.h"

namespace Zuck { namespace read_first {

struct VariantDetail {
    std::string repeat, spacer, proof;
    std::size_t templates = 0, full_pairs = 0, exact_templates = 0, short_templates = 0;
    // Placement-specific, canonical family orientation, one-based positions.
    // A dash means both observed repeat segments agree with the model.
    std::set<std::string> substitution_signatures;
};

namespace variant_junctions {
using Unit = read_recruitment::detail::Unit;
using Key = extraction::Key;
constexpr int seed_length = 12;

struct Model {
    std::size_t family;
    bool reverse;
    std::string repeat;
};

struct Flank {
    int length = 0, change = -1, exact_junction = 0;
    char observed = 'N';
};

struct Placement {
    std::size_t strand;
    int begin, end;
    Flank left, right;
};

struct Witness {
    bool full_pair = false, well_bounded = false, single = false, exact = false, clipped = false;
    int minimum_offset = std::numeric_limits<int>::max(), maximum_offset = -1;
    std::set<std::string> signatures;
};

inline void merge(Witness& a, const Witness& b) {
    a.full_pair |= b.full_pair; a.well_bounded |= b.well_bounded;
    a.single |= b.single; a.exact |= b.exact; a.clipped |= b.clipped;
    a.minimum_offset = std::min(a.minimum_offset, b.minimum_offset);
    a.maximum_offset = std::max(a.maximum_offset, b.maximum_offset);
    a.signatures.insert(b.signatures.begin(), b.signatures.end());
}

struct Lengths {
    int minimum = 51, maximum = 25;
    bool contains(int size) const { return minimum <= size && size <= maximum; }
};

inline std::vector<Lengths> supported_lengths(const std::vector<Family>& families,
                                              const std::vector<Call>& accepted) {
    std::map<std::string, std::vector<const Call*>> by_family;
    for (const auto& call : accepted) by_family[call.repeat].push_back(&call);
    std::vector<Lengths> lengths(families.size());
    for (std::size_t i = 0; i < families.size(); ++i) {
        if (families[i].templates < 5 || families[i].contexts < 4) continue;
        const auto& calls = by_family[families[i].repeat];
        bool dissimilar = false;
        for (const auto* a : calls) {
            if (read_recruitment::detail::edit_distance(a->spacer,calls.front()->spacer) * 5 >=
                std::max(a->spacer.size(),calls.front()->spacer.size())) dissimilar = true;
        }
        if (!dissimilar) continue;
        for (const auto* a : calls) {
            const int size = static_cast<int>(a->spacer.size());
            lengths[i].minimum = std::min(lengths[i].minimum, std::max(26,size-4));
            lengths[i].maximum = std::max(lengths[i].maximum, std::min(50,size+4));
        }
    }
    return lengths;
}

// Compare the entire repeat tract that the physical read covers. A read-end
// clip can shorten a tract; an internal mismatch cannot silently shorten it.
inline Flank confirm(const klibpp::KSeq& read, const std::string& repeat,
                     int boundary, bool left) {
    const int size = static_cast<int>(repeat.size());
    const int n = std::min(size, left ? boundary : static_cast<int>(read.seq.size())-boundary);
    if (n < 6) return {};
    Flank flank{n,-1,n,'N'};
    for (int d = 0; d < n; ++d) {
        const int r = left ? size-1-d : d;
        const int p = left ? boundary-1-d : boundary+d;
        const char base = read.seq[p];
        if (read_evidence::base(base) < 0) return {};
        if (base == repeat[r]) continue;
        // Require long flanks for variants; the spacer-facing eight bases
        // and every base of a short terminal repeat must agree literally.
        if (flank.change >= 0 || d < 8 || n < 16 || read.qual[p] < 48) return {};
        flank.change = r; flank.observed = base; flank.exact_junction = d;
    }
    return flank;
}

inline std::string substitution(const Flank& f, const Model& m, char side) {
    if (f.change < 0) return {};
    const int pos = m.reverse ? static_cast<int>(m.repeat.size())-1-f.change : f.change;
    char reference = m.repeat[f.change], observed = f.observed;
    if (m.reverse) {
        reference = "TGCA"[read_evidence::base(reference)];
        observed = "TGCA"[read_evidence::base(observed)];
    }
    return std::string(1,side) + std::to_string(pos+1) + reference + ">" + observed;
}

class Scanner {
    std::vector<Model> strands_;
    // A perfect internal seed proposes an oriented repeat coordinate. Different
    // seeds giving the same coordinate trigger only one physical verification.
    std::vector<std::pair<std::size_t,int>> offsets_;
    read_evidence::ExactMatcher matcher_;
    const std::vector<Lengths>& lengths_;
public:
    Scanner(const std::vector<Family>& families, const std::vector<Lengths>& lengths) : lengths_(lengths) {
        for (std::size_t i = 0; i < families.size(); ++i) {
            const auto& repeat = families[i].repeat;
            if (repeat.size() < 20 || repeat.size() > 55) continue;
            const auto add = [&](const std::string& seq, bool reverse) {
                const auto strand = strands_.size();
                strands_.push_back({i,reverse,seq});
                const int last = static_cast<int>(seq.size()) - seed_length - 1;
                for (int p = 1;; p = std::min(p+6,last)) {
                    matcher_.add(seq.substr(p,seed_length), offsets_.size());
                    offsets_.emplace_back(strand,p);
                    if (p == last) break;
                }
            };
            add(repeat,false);
            const auto rc = read_linked::reverse_complement(repeat);
            if (rc != repeat) add(rc,true);
        }
        matcher_.build();
    }

    bool scan(const klibpp::KSeq& read, std::map<Key,Witness>& pending,
              std::map<std::string,std::size_t>& counts) const {
        if (read.qual.empty() || read.seq.size() < 26+12) return true;
        if (read.seq.size() > read_recruitment::Config{}.maximum_read_length) return false;
        const int length = static_cast<int>(read.seq.size());
        std::set<std::pair<std::size_t,int>> starts;
        std::size_t hits = 0;
        matcher_.scan(read.seq, [&](std::size_t id, std::size_t end) {
            if (++hits > 4096) return;
            const auto& seed = offsets_[id];
            starts.emplace(seed.first, static_cast<int>(end)-seed_length-seed.second);
        });
        if (hits > 4096 || starts.size() > 256) return false;
        read_evidence::QualitySpans quality(read.qual,10), high(read.qual,15);
        std::set<std::tuple<std::size_t,int,int>> seen;
        std::vector<Placement> placements;
        for (const auto& start : starts) {
            const auto& model = strands_[start.first];
            const int size = static_cast<int>(model.repeat.size());
            for (bool upstream : {true,false}) for (int gap = 26; gap <= 50; ++gap) {
                const int begin = upstream ? start.second+size : start.second-gap;
                const int end = begin+gap;
                if (begin < 6 || end+6 > length || !seen.emplace(start.first,begin,end).second) continue;
                if (seen.size() > 8192) return false;
                const auto a = confirm(read,model.repeat,begin,true);
                if (!a.length) continue;
                const auto b = confirm(read,model.repeat,end,false);
                if (!b.length) continue;
                // A 6-9 nt terminal hook can only corroborate an exact, complete
                // repeat on the other side; it never establishes a cut alone.
                if ((a.length < 10 && !(b.length == size && b.change < 0)) ||
                    (b.length < 10 && !(a.length == size && a.change < 0))) continue;
                if (a.length < size && b.length < size && std::min(a.length,b.length) < 12) continue;
                if (read.seq.find_first_not_of("ACGT",begin) < static_cast<std::size_t>(end)) continue;
                if (!quality.passes(begin-std::min(12,a.length),end+std::min(12,b.length),0.75)) continue;
                placements.push_back({start.first,begin,end,a,b});
                if (placements.size() > 256) return false;
            }
        }
        for (const auto& p : placements) {
            const auto& model = strands_[p.strand];
            if (!lengths_[model.family].contains(p.end-p.begin)) continue;
            if (std::any_of(placements.begin(),placements.end(),[&](const auto& q) {
                return &p != &q && p.begin < q.end && q.begin < p.end;
            })) {
                ++counts["read_first_variant_ambiguous_placements"];
                continue;
            }
            const int size = static_cast<int>(model.repeat.size());
            const bool full_a = p.left.length == size, full_b = p.right.length == size;
            const auto sp = read.seq.substr(p.begin,p.end-p.begin);
            const bool palindrome = model.repeat == read_linked::reverse_complement(model.repeat);
            // A self-complementary repeat cannot assign the strand of an
            // altered copy consistently enough to pool directional alleles.
            if (palindrome) continue;
            auto& w = pending[{model.family,model.reverse ? read_linked::reverse_complement(sp) : sp}];
            w.full_pair |= full_a && full_b;
            const int shorter = std::min(p.left.length,p.right.length);
            w.well_bounded |= (full_a || full_b) && shorter >= 10;
            const bool pristine_full = (full_a && p.left.change < 0) || (full_b && p.right.change < 0);
            w.single |= pristine_full && shorter >= 16 && p.left.exact_junction >= 12 &&
                        p.right.exact_junction >= 12 && high.passes(p.begin-12,p.end+12,0.5);
            w.exact |= p.left.change < 0 && p.right.change < 0;
            w.clipped |= shorter < 10;
            const int offset = model.reverse ? length-p.end : p.begin;
            w.minimum_offset = std::min(w.minimum_offset,offset);
            w.maximum_offset = std::max(w.maximum_offset,offset);
            std::string left = substitution(model.reverse ? p.right : p.left,model,'L');
            std::string right = substitution(model.reverse ? p.left : p.right,model,'R');
            if (!left.empty() && !right.empty()) left += ',';
            w.signatures.insert(left.empty() && right.empty() ? "-" : left+right);
        }
        return true;
    }
};

} // namespace variant_junctions

inline std::vector<Call> recruit_variants(const std::vector<Family>& families,
        const std::vector<std::string>& paths, const std::vector<Call>& strict,
        std::map<std::string,std::size_t>& counts, std::vector<VariantDetail>* details = nullptr) {
    using namespace variant_junctions;
    if (details) details->clear();
    const auto lengths = supported_lengths(families,strict);
    const auto eligible = std::count_if(lengths.begin(),lengths.end(),[](const auto& l) { return l.minimum <= l.maximum; });
    counts["read_first_variant_eligible_families"] = eligible;
    if (!eligible || paths.empty()) return {};
    Scanner scanner(families,lengths);
    std::map<Key,std::map<Unit,Witness>> inventory;
    std::map<std::pair<std::string,std::string>,std::size_t> known;
    for (const auto& call : strict) known[{call.repeat,call.spacer}] = call.templates;
    bool overflow = false;
    std::size_t memberships = 0, signature_slots = 0;
    seeds::each_template(paths,[&](const klibpp::KSeq& a,const klibpp::KSeq* b) {
        if (overflow) return;
        std::map<Key,Witness> pending;
        if (!scanner.scan(a,pending,counts) || (b && !scanner.scan(*b,pending,counts))) {
            ++counts["read_first_variant_limited_fragments"];
            return;
        }
        if (pending.empty()) return;
        const auto id = read_recruitment::detail::unit_identity(a.seq,b ? b->seq : "",b != nullptr);
        for (const auto& item : pending) {
            auto& family_calls = inventory[item.first];
            auto observed = family_calls.emplace(id,Witness{});
            memberships += observed.second;
            auto& w = observed.first->second;
            const auto ns = w.signatures.size();
            merge(w,item.second);
            signature_slots += w.signatures.size()-ns;
        }
        if (inventory.size() > 32768 || memberships > 250000 || signature_slots > 400000) {
            overflow = true;
            inventory.clear();
        }
    });
    counts["read_first_variant_global_resource_limit"] = overflow;
    counts["read_first_variant_candidate_types"] = inventory.size();
    std::vector<Call> result;
    for (const auto& item : inventory) {
        const auto& repeat = families[item.first.first].repeat;
        const auto& spacer = item.first.second;
        if (known.count({repeat,spacer})) continue;
        VariantDetail detail;
        detail.repeat = repeat; detail.spacer = spacer; detail.templates = item.second.size();
        Witness pooled;
        for (const auto& unit : item.second) {
            const auto& w = unit.second;
            detail.full_pairs += w.full_pair;
            detail.exact_templates += w.exact;
            detail.short_templates += w.clipped;
            merge(pooled,w);
        }
        if ((detail.templates == 1 && !pooled.single) || !pooled.well_bounded ||
            (detail.templates > 1 && pooled.minimum_offset == pooled.maximum_offset)) {
            ++counts["read_first_variant_insufficient_geometry"];
            continue;
        }
        // Variants of a well-supported spacer differing by one nucleotide can
        // be substitution errors. Preserve recurrent alternatives only when
        // neither strict nor rescued observations dominate their support.
        bool neighbor = false;
        auto key = item.first;
        for (std::size_t pos = 0; !neighbor && pos < spacer.size(); ++pos) {
            for (char alt : std::string("ACGT")) {
                if (alt == spacer[pos]) continue;
                key.second[pos] = alt;
                const auto anchor = known.find({repeat,key.second});
                const auto other = inventory.find(key);
                std::size_t support = anchor == known.end() ? 0 : anchor->second;
                if (other != inventory.end()) support = std::max(support,other->second.size());
                if (support && (detail.templates == 1 || support >= 2*detail.templates)) {
                    neighbor = true;
                    break;
                }
            }
            key.second[pos] = spacer[pos];
        }
        if (neighbor || !seeds::admissible_interior(repeat,spacer) || extraction::tandem_explains(repeat,spacer)) {
            ++counts["read_first_variant_neighbor_or_tandem"];
            continue;
        }
        detail.proof = detail.templates == 1 ? "one_exact_full_repeat_high_quality" :
                                              "independent_bilateral_spacer_templates";
        detail.substitution_signatures = std::move(pooled.signatures);
        result.push_back({repeat,spacer,detail.templates,detail.full_pairs});
        if (details) details->push_back(std::move(detail));
    }
    counts["read_first_variant_additions"] = result.size();
    return result;
}

}} // namespace Zuck::read_first
#endif
