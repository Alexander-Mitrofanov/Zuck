#ifndef ZUCK_READ_FIRST_H
#define ZUCK_READ_FIRST_H

#include "crispr_array/read_first_clipped_junction.h"
#include "crispr_array/read_recruitment.h"

namespace Zuck { namespace read_first {

struct Call {
    std::string repeat, spacer;
    std::size_t templates, full_pairs;
    // Spacer positions observed at Q30 on an intact repeat-spacer-repeat read.
    std::uint64_t whole_q30_bases = 0;
};

namespace extraction {

constexpr std::size_t anchor = 10;
constexpr std::size_t double_partial_anchor = 12;
using Unit = read_recruitment::detail::Unit;
using Key = std::pair<std::size_t, std::string>;  // Family and family-oriented spacer.

struct Strand {
    std::size_t family;
    bool reverse;
    std::string repeat;
};

struct Placement {
    std::size_t strand, begin, end, left, right;
};

struct Witness {
    bool full_pair = false;
    std::uint64_t whole_q30_bases = 0;
    bool one_full = false;
    bool high_quality = false;
    std::size_t minimum_offset = std::numeric_limits<std::size_t>::max(), maximum_offset = 0;
    std::string quality;  // In family orientation; retained only for singleton error screening.
};

struct Evidence {
    std::map<Unit, bool> units;  // True if either mate contains both complete repeat copies.
    std::uint64_t whole_q30_bases = 0;
    bool one_full = false;
    bool high_quality = false;
    std::size_t minimum_offset = std::numeric_limits<std::size_t>::max(), maximum_offset = 0;
    std::string singleton_quality;
};

inline void merge_quality(std::string& target, const std::string& value) {
    if (target.empty()) { target = value; return; }
    for (std::size_t i = 0; i < value.size(); ++i)
        target[i] = std::max(target[i], value[i]);
}

class Anchors {
    read_evidence::ExactMatcher matcher_;
    std::vector<Strand> strands_;

public:
    explicit Anchors(const std::vector<Family>& families) {
        for (std::size_t i = 0; i < families.size(); ++i) {
            const auto& family = families[i];
            if (family.templates < 2 || family.contexts < 2 || family.repeat.size() < 2 * anchor)
                continue;
            const auto reverse = read_linked::reverse_complement(family.repeat);
            if (reverse.empty()) continue;
            const auto add = [&](const std::string& repeat, bool reversed) {
                const auto id = strands_.size();
                strands_.push_back({i, reversed, repeat});
                matcher_.add(repeat.substr(repeat.size() - anchor), 2 * id);
                matcher_.add(repeat.substr(0, anchor), 2 * id + 1);
            };
            add(family.repeat, false);
            if (reverse != family.repeat) add(reverse, true);
        }
        matcher_.build();
    }

    const auto& strands() const { return strands_; }

    bool scan(const klibpp::KSeq& read, std::map<Key, Witness>& pending,
              std::map<std::string, std::size_t>& counts, const clipped::Collector* clipped_collector,
              std::vector<clipped::Candidate>& clipped_candidates) const {
        const auto& sequence = read.seq;
        if (read.qual.empty() || sequence.size() < 2 * anchor + read_recruitment::MINIMUM_SPACER)
            return true;
        if (sequence.size() > read_recruitment::Config{}.maximum_read_length) return false;
        using Site = std::pair<std::size_t, std::size_t>;  // Strand ID and spacer-facing cut.
        std::map<Site, std::size_t> left, right;
        std::size_t hits = 0;
        matcher_.scan(sequence, [&](std::size_t id, std::size_t end) {
            ++hits;
            if (hits > 4096) return;
            const auto strand = id / 2;
            const auto& repeat = strands_[strand].repeat;
            if (id % 2 == 0) {
                const auto length = std::min(end, repeat.size());
                if (sequence.compare(end - length, length, repeat, repeat.size() - length, length) == 0)
                    left.emplace(Site{strand, end}, length);
            } else {
                const auto begin = end - anchor;
                const auto length = std::min(sequence.size() - begin, repeat.size());
                if (sequence.compare(begin, length, repeat, 0, length) == 0)
                    right.emplace(Site{strand, begin}, length);
            }
        });
        if (hits > 4096) return false;
        if (left.empty() || right.empty()) return true;
        std::vector<Placement> placements;
        read_evidence::QualitySpans quality(read.qual, read_recruitment::MINIMUM_QUALITY);
        for (const auto& item : left) {
            const auto strand = item.first.first;
            const auto begin = item.first.second;
            const auto size = strands_[strand].repeat.size();
            auto it = right.lower_bound({strand, begin + read_recruitment::MINIMUM_SPACER});
            for (; it != right.end() && it->first.first == strand; ++it) {
                const auto end = it->first.second;
                if (end - begin > read_recruitment::MAXIMUM_SPACER) break;
                const auto l = item.second, r = it->second;
                if (l < size && r < size && std::min(l, r) < double_partial_anchor) continue;
                if (sequence.find_first_not_of("ACGT", begin) < end ||
                    !quality.passes(begin - l, end + r, read_recruitment::MAXIMUM_EXPECTED_ERRORS)) continue;
                placements.push_back({strand, begin, end, l, r});
                if (placements.size() > 512) return false;
            }
        }

        if (clipped_collector)
            clipped_collector->scan(read,left,right,strands_,placements,clipped_candidates);

        read_evidence::QualitySpans singleton_quality(read.qual, 15);
        for (const auto& p : placements) {
            // Competing overlapping spacer intervals and ambiguous family ownership
            // cannot be resolved by assigning the first matching repeat.
            bool competing = false;
            for (const auto& q : placements) {
                if (&p == &q || p.end <= q.begin || q.end <= p.begin) continue;
                competing = true;
                break;
            }
            if (competing) {
                ++counts["read_first_ambiguous_placements"];
                continue;
            }
            const auto& strand = strands_[p.strand];
            auto spacer = sequence.substr(p.begin, p.end - p.begin);
            auto qualities = read.qual.substr(p.begin, p.end - p.begin);
            const bool reverse = strand.reverse || (strand.repeat == read_linked::reverse_complement(strand.repeat) &&
                                                     read_linked::reverse_complement(spacer) < spacer);
            if (reverse) {
                spacer = read_linked::reverse_complement(spacer);
                std::reverse(qualities.begin(), qualities.end());
            }
            auto& witness = pending[{strand.family, spacer}];
            const auto size = strand.repeat.size();
            if (p.left == size && p.right == size) {
                witness.full_pair = true;
                static_assert(read_recruitment::MAXIMUM_SPACER <= 64);
                for (std::size_t i = 0; i < qualities.size(); ++i)
                    if (static_cast<unsigned char>(qualities[i]) >= 33 + 30)
                        witness.whole_q30_bases |= std::uint64_t{1} << i;
            }
            const bool one_full = p.left == size || p.right == size;
            witness.one_full |= one_full;
            witness.high_quality |= one_full && singleton_quality.passes(p.begin, p.end, 0.5);
            const auto offset = reverse ? sequence.size() - p.end : p.begin;
            witness.minimum_offset = std::min(witness.minimum_offset, offset);
            witness.maximum_offset = std::max(witness.maximum_offset, offset);
            merge_quality(witness.quality, qualities);
        }
        return true;
    }
};

inline bool likely_substitution(const Key& key, const Evidence& evidence,
                                const std::map<Key, Evidence>& inventory) {
    if (evidence.units.size() != 1) return false;
    auto neighbor = key;
    for (std::size_t i = 0; i < key.second.size(); ++i) {
        const int q = static_cast<unsigned char>(evidence.singleton_quality[i]) - 33;
        const double per_base = std::pow(10.0, -0.1 * q) / 3.0;
        for (char base : std::string("ACGT")) {
            if (base == key.second[i]) continue;
            neighbor.second[i] = base;
            const auto found = inventory.find(neighbor);
            if (found != inventory.end() && found->second.units.size() >= 2 &&
                found->second.units.size() * per_base >= 0.05) return true;
        }
        neighbor.second[i] = key.second[i];
    }
    return false;
}

inline bool tandem_explains(const std::string& repeat, const std::string& spacer) {
    // Test the proposed family repeat and spacer against native tandem models.
    // This hypothesis check does not contribute observed repeat bases or template counts.
    read_recruitment::Candidate proposed;
    proposed.left_repeat = repeat;
    proposed.right_repeat = repeat;
    proposed.spacer = spacer;
    proposed.left_observed_bases = proposed.right_observed_bases = repeat.size();
    return read_recruitment::tandem_evidence(proposed).explained;
}

// A well-sampled family made entirely of closely related spacer interiors is
// also compatible with copies of a longer tandem unit. Any dissimilar observed
// spacer protects the family, including its closely related minority alleles.
inline std::set<std::size_t> repeated_units(const std::vector<Family>& families,
                                           const std::map<Key, Evidence>& inventory) {
    using Member = std::map<Key, Evidence>::value_type;
    std::map<std::size_t, std::vector<const Member*>> groups;
    for (const auto& entry : inventory)
        if (families[entry.first.first].templates >= 8) groups[entry.first.first].push_back(&entry);
    std::set<std::size_t> ambiguous;
    for (const auto& group : groups) {
        const auto& members = group.second;
        // Avoid quadratic comparisons in large families; this gate is only
        // for small homogeneous interiors supported across many templates.
        if (members.size() < 2 || members.size() > 16) continue;
        std::set<Unit> full;
        for (const auto* member : members) {
            for (const auto& u : member->second.units) {
                if (u.second) full.insert(u.first);
                if (full.size() >= 4) break;
            }
            if (full.size() >= 4) break;
        }
        if (full.size() < 4) continue;
        bool dissimilar = false;
        for (std::size_t i = 0; i < members.size() && !dissimilar; ++i) {
            const auto& a = members[i]->first.second;
            for (std::size_t j = i + 1; j < members.size(); ++j) {
                const auto& b = members[j]->first.second;
                if (50 * read_recruitment::detail::edit_distance(a, b) >= 9 * std::max(a.size(), b.size())) {
                    dissimilar = true;
                    break;
                }
            }
        }
        if (!dissimilar) ambiguous.insert(group.first);
    }
    return ambiguous;
}

} // namespace extraction

inline std::vector<Call> recruit(const std::vector<Family>& families, const std::vector<std::string>& paths,
                                 std::map<std::string, std::size_t>& counts, bool shallow_rescue = false,
                                 std::vector<Call>* clipped_calls = nullptr) {
    using namespace extraction;
    if (families.empty() || paths.empty()) return {};
    if (paths.size() > 2) throw std::invalid_argument("Read-first recruitment expects one or two inputs");
    Anchors anchors(families);
    std::unique_ptr<clipped::Collector> collector;
    if (shallow_rescue && clipped_calls) collector = std::make_unique<clipped::Collector>(families);
    std::map<Key, Evidence> inventory;
    std::size_t memberships = 0;
    bool overflow = false;
    read_linked::Reader first(paths[0]);
    std::unique_ptr<read_linked::Reader> second;
    if (paths.size() == 2) second = std::make_unique<read_linked::Reader>(paths[1]);
    klibpp::KSeq a, b;
    while (first.next(a)) {
        if (second && (!second->next(b) || read_linked::read_name(a.name) != read_linked::read_name(b.name)))
            throw std::runtime_error("Read-first recruitment received unpaired reads");
        ++counts["read_first_recruit_input_templates"];
        counts["read_first_recruit_no_quality"] += a.qual.empty() + (second && b.qual.empty());
        if (overflow) continue;
        std::map<Key, Witness> pending;
        std::vector<clipped::Candidate> first_cuts, second_cuts;
        if (!anchors.scan(a,pending,counts,collector.get(),first_cuts) ||
            (second && !anchors.scan(b,pending,counts,collector.get(),second_cuts))) {
            ++counts["read_first_recruit_limited_templates"];
            continue;
        }
        if (pending.empty() && first_cuts.empty() && second_cuts.empty()) continue;
        const auto unit = read_recruitment::detail::unit_identity(a.seq, second ? b.seq : "", bool(second));
        if (!first_cuts.empty() || !second_cuts.empty()) {
            const auto fragment = seeds::seed_fragment(a,second ? &b : nullptr);
            collector->store(a,anchors.strands(),first_cuts,unit,fragment);
            if (second) collector->store(b,anchors.strands(),second_cuts,unit,fragment);
        }
        for (const auto& item : pending) {
            auto& evidence = inventory[item.first];
            const auto& witness = item.second;
            auto added = evidence.units.emplace(unit, witness.full_pair);
            added.first->second |= witness.full_pair;
            memberships += added.second;
            evidence.one_full |= witness.one_full;
            evidence.whole_q30_bases |= witness.whole_q30_bases;
            evidence.high_quality |= witness.high_quality;
            evidence.minimum_offset = std::min(evidence.minimum_offset, witness.minimum_offset);
            evidence.maximum_offset = std::max(evidence.maximum_offset, witness.maximum_offset);
            if (evidence.units.size() == 1) merge_quality(evidence.singleton_quality, witness.quality);
            else evidence.singleton_quality.clear();
        }
        if (inventory.size() > 65536 || memberships > 1000000) {
            overflow = true;
            inventory.clear();  // Avoid input-order-dependent calls from incomplete evidence.
        }
    }
    if (second && second->next(b)) throw std::runtime_error("Read-first recruitment received extra mate reads");

    std::vector<Call> calls;
    counts["read_first_recruit_global_resource_limit"] = overflow;
    counts["read_first_recruit_candidates"] = inventory.size();
    const auto repetitive_families = repeated_units(families, inventory);
    counts["read_first_repeated_unit_families"] = repetitive_families.size();
    for (const auto& item : inventory) {
        const auto& evidence = item.second;
        const auto& family = families[item.first.first];
        const auto& repeat = family.repeat;
        const auto& spacer = item.first.second;
        if (repetitive_families.count(item.first.first)) {
            ++counts["read_first_repeated_unit_candidates"];
            continue;
        }
        const bool strong_family = family.contexts >= 4 && family.templates >= 5;
        // Every placement has already passed Q10 and <=1 expected error over
        // the entire physical repeat-spacer-repeat span. Families with limited
        // evidence retain the additional Q15 and <=0.5 spacer error criterion.
        if (evidence.units.size() == 1 && (!evidence.one_full || (!evidence.high_quality && !strong_family))) {
            ++counts["read_first_weak_singletons"];
            continue;
        }
        if (!evidence.one_full && evidence.minimum_offset == evidence.maximum_offset) {
            ++counts["read_first_double_partial_single_offset"];
            continue;
        }
        if (likely_substitution(item.first, evidence, inventory)) {
            ++counts["read_first_substitution_neighbors"];
            continue;
        }
        if (tandem_explains(repeat, spacer)) {
            ++counts["read_first_tandem_explained"];
            continue;
        }
        const auto full = std::count_if(evidence.units.begin(), evidence.units.end(),
                                        [](const auto& u) { return u.second; });
        calls.push_back({repeat, spacer, evidence.units.size(), static_cast<std::size_t>(full),
                         evidence.whole_q30_bases});
        counts["read_first_one_template_calls"] += evidence.units.size() == 1;
        counts["read_first_family_calibrated_singletons"] += evidence.units.size() == 1 && !evidence.high_quality;
        counts["read_first_double_partial_only_calls"] += !evidence.one_full;
    }
    if (collector && !overflow) {
        std::map<std::size_t,std::set<std::string>> pool;
        for (const auto& known : inventory) pool[known.first.first].insert(known.first.second);
        for (const auto& item : collector->candidates()) pool[item.first.first].insert(item.first.second);
        for (const auto& item : collector->candidates()) {
            const auto family = item.first.first;
            const auto& spacer = item.first.second;
            if (inventory.count(item.first) || repetitive_families.count(family)) continue;
            if (item.second.size() == 1 && clipped::near_known(spacer,pool.at(family))) {
                ++counts["read_first_clipped_unconfirmed_alleles"];
                continue;
            }
            if (!tandem_explains(families[family].repeat,spacer))
                clipped_calls->push_back({families[family].repeat,spacer,item.second.size(),0});
        }
        counts["read_first_clipped_candidates"] = collector->candidates().size();
        counts["read_first_clipped_proposals"] = clipped_calls->size();
        counts["read_first_clipped_resource_limit"] = collector->limited();
    }
    counts["read_first_recruit_calls"] = calls.size();
    return calls;
}

}} // namespace Zuck::read_first
#endif
