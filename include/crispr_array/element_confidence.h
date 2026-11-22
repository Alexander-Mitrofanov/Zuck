#ifndef ZUCK_ELEMENT_CONFIDENCE_H
#define ZUCK_ELEMENT_CONFIDENCE_H

#include <limits>
#include "crispr_array/array_database.h"

namespace Zuck { namespace element_confidence {

// These are bounded evidence indices, not calibrated correctness probabilities.
// A cap is inclusive: a displayed 32 means at least 32 distinct observations.
constexpr std::size_t SUPPORT_CAP = 32;
constexpr std::size_t SCORE_UNIT_SATURATION = 4;
constexpr std::size_t SCORE_OFFSET_SATURATION = 2;

inline double unknown() { return std::numeric_limits<double>::quiet_NaN(); }
inline std::string circular_identity(const std::string& s) {
    return std::min(array_database::circular_key(s),
        array_database::circular_key(array_database::reverse_complement(s)));
}
inline std::pair<std::string, std::string> pair_identity(
        const std::pair<std::string, std::string>& p) {
    return std::min(p, std::make_pair(array_database::reverse_complement(p.first),
                                     array_database::reverse_complement(p.second)));
}

struct Support {
    std::size_t sequence_units = 0, start_offsets = 0;
    bool units_saturated = false, offsets_saturated = false, observable = false;
    bool accepted() const { return observable && sequence_units >= 2 && start_offsets >= 2; }
    double index() const {
        if (!observable) return unknown();
        return static_cast<double>(std::min(sequence_units, SCORE_UNIT_SATURATION)) /
            SCORE_UNIT_SATURATION *
            static_cast<double>(std::min(start_offsets, SCORE_OFFSET_SATURATION)) /
            SCORE_OFFSET_SATURATION;
    }
};

struct Row {
    Support repeat, spacer, context;
    double repeat_score = unknown(), spacer_score = unknown(), context_score = unknown();
    double repeat_database_score = 0, spacer_database_score = 0;
    std::size_t repeat_context_core_count = 0;
    double repeat_family_score = unknown();
    std::string database_source;
    double database_score = 0;
    std::size_t database_other_spacers = 0, database_other_cores = 0;
    std::size_t boundary_alternatives = 0, boundary_database_alternatives = 0;
    double boundary_score = unknown();
    std::string boundary_status, context_observability;
};
struct Result {
    std::vector<Row> rows;
    std::map<std::string, std::size_t> counts;
};

// Optional family-aware context validation. The complete upstream repeat and
// spacer remain exact. Only the following 16-base anchor may come from another
// observed allele in the same output family. Same-length alleles must differ
// at no more than four positions; reverse-strand rows are oriented first.
// This does not borrow repeat alleles across families or invent a consensus.
template<class Arrays>
std::vector<std::vector<std::string>> family_repeat_prefixes(const Arrays& arrays,
        const std::set<std::pair<std::size_t, std::size_t>>& excluded_allele_rows = {}) {
    for (const auto& row : excluded_allele_rows)
        if (row.first >= arrays.size() || row.second == 0 || row.second > arrays[row.first].entries.size())
            throw std::invalid_argument("Family context exclusion references an absent row");
    std::vector<std::vector<std::string>> result;
    for (std::size_t owner = 0; owner < arrays.size(); ++owner) {
        const auto& array = arrays[owner];
        std::set<std::string> alleles;
        for (std::size_t entry = 0; entry < array.entries.size(); ++entry)
            if (!excluded_allele_rows.count({owner, entry + 1})) alleles.insert(array.entries[entry].first);
        // Keep all annotation rows aligned. Evidence-only novel rows retain
        // their own exact prefix, but cannot lend a new allele to an old row.
        for (const auto& pair : array.entries) {
            const auto& repeat = pair.first;
            if (!array_database::dna(repeat) || repeat.size() < database_reads::NEXT_REPEAT_ANCHOR)
                throw std::invalid_argument("Invalid family context repeat");
            std::set<std::string> prefixes{repeat.substr(0, database_reads::NEXT_REPEAT_ANCHOR)};
            for (const auto& allele : alleles) {
                if (!array_database::dna(allele))
                    throw std::invalid_argument("Invalid family context repeat allele");
                if (allele.size() != repeat.size()) continue;
                for (const auto& oriented : {allele, array_database::reverse_complement(allele)}) {
                    std::size_t mismatches = 0;
                    for (std::size_t p = 0; p < repeat.size() && mismatches <= 4; ++p)
                        mismatches += repeat[p] != oriented[p];
                    if (mismatches <= 4)
                        prefixes.insert(oriented.substr(0, database_reads::NEXT_REPEAT_ANCHOR));
                }
            }
            result.emplace_back(prefixes.begin(), prefixes.end());
        }
    }
    return result;
}

inline Result annotate(const std::vector<std::pair<std::string, std::string>>& pairs,
                       const std::vector<std::string>& files,
                       const std::vector<array_database::Record>& records = {},
                       const std::vector<std::vector<std::string>>& next_repeat_prefixes = {}) {
    if (files.size() > 2) throw std::runtime_error("Confidence scoring accepts at most two read files");
    if (!next_repeat_prefixes.empty() && next_repeat_prefixes.size() != pairs.size())
        throw std::invalid_argument("Family context prefix row count mismatch");
    Result result;
    result.rows.resize(pairs.size());
    result.counts = {{"confidence_rows", pairs.size()}, {"confidence_read_records", 0},
        {"confidence_quality_read_records", 0}, {"confidence_input_units", 0},
        {"confidence_context_supported_rows", 0}, {"confidence_patterns", 0},
        {"confidence_capped_patterns", 0}, {"confidence_ambiguous_rows", 0}};
    struct Pattern {
        std::string sequence;
        std::set<std::string> units;
        std::set<std::size_t> offsets;
        std::vector<std::string> alternatives;
    };
    std::vector<Pattern> patterns;
    std::map<std::string, std::size_t> pattern_ids;
    auto add_pattern = [&](const std::string& sequence) {
        if (!array_database::dna(sequence))
            throw std::runtime_error("Confidence patterns must contain only ACGT");
        const auto key = array_database::canonical(sequence);
        const auto entry = pattern_ids.emplace(key, patterns.size());
        if (entry.second) patterns.push_back({key, {}, {}, {key}});
        return entry.first->second;
    };
    auto add_context = [&](const std::string& repeat, const std::string& spacer, std::size_t row) {
        const auto original = repeat + spacer + repeat.substr(0, database_reads::NEXT_REPEAT_ANCHOR);
        if (next_repeat_prefixes.empty()) return add_pattern(original);
        std::set<std::string> sequences{array_database::canonical(original)};
        for (const auto& prefix : next_repeat_prefixes[row]) {
            if (prefix.size() != database_reads::NEXT_REPEAT_ANCHOR || !array_database::dna(prefix))
                throw std::invalid_argument("Invalid family context next-repeat anchor");
            sequences.insert(array_database::canonical(repeat + spacer + prefix));
        }
        if (sequences.size() == 1) return add_pattern(*sequences.begin());
        ++result.counts["confidence_family_context_alternative_rows"];
        std::string key;
        for (const auto& sequence : sequences) { key += sequence; key += '\0'; }
        const auto entry = pattern_ids.emplace(key, patterns.size());
        if (entry.second)
            patterns.push_back({*sequences.begin(), {}, {}, {sequences.begin(), sequences.end()}});
        return entry.first->second;
    };
    std::vector<std::array<std::size_t, 3>> row_patterns;
    std::vector<std::string> cores, canonical_repeats, canonical_spacers;
    std::map<std::string, std::vector<std::size_t>> by_repeat, by_spacer;
    using Identity = std::pair<std::string, std::string>;
    std::vector<Identity> identities;
    std::map<std::string, std::map<Identity, double>> alternatives;
    for (std::size_t i = 0; i < pairs.size(); ++i) {
        const auto& [repeat, spacer] = pairs[i];
        if (repeat.size() < database_reads::NEXT_REPEAT_ANCHOR || spacer.empty())
            throw std::runtime_error("Invalid confidence repeat/spacer span");
        row_patterns.push_back({add_pattern(repeat), add_pattern(spacer), add_context(repeat, spacer, i)});
        cores.push_back(circular_identity(repeat + spacer));
        identities.push_back(pair_identity(pairs[i]));
        alternatives[cores.back()].emplace(identities.back(), 0.0);
        canonical_repeats.push_back(array_database::canonical(repeat));
        canonical_spacers.push_back(array_database::canonical(spacer));
        by_repeat[canonical_repeats.back()].push_back(i);
        by_spacer[canonical_spacers.back()].push_back(i);
    }
    result.counts["confidence_patterns"] = patterns.size();
    std::size_t max_quality_read_length = 0;
    if (!files.empty()) {
        read_evidence::ExactMatcher matcher;
        for (std::size_t i = 0; i < patterns.size(); ++i) {
            for (const auto& sequence : patterns[i].alternatives) {
                matcher.add(sequence, i);
                const auto reverse = array_database::reverse_complement(sequence);
                if (reverse != sequence) matcher.add(reverse, i);
            }
        }
        matcher.build();
        auto scan = [&](const klibpp::KSeq& read, std::map<std::size_t, std::set<std::size_t>>& hits) {
            ++result.counts["confidence_read_records"];
            if (read.qual.empty()) return;
            ++result.counts["confidence_quality_read_records"];
            max_quality_read_length = std::max(max_quality_read_length, read.seq.size());
            read_evidence::QualitySpans quality(read.qual, database_reads::MINIMUM_QUALITY);
            matcher.scan(read.seq, [&](std::size_t id, std::size_t end) {
                const auto& pattern = patterns[id];
                if (pattern.units.size() >= SUPPORT_CAP && pattern.offsets.size() >= SUPPORT_CAP)
                    return;
                const auto start = end - pattern.sequence.size();
                if (quality.passes(start, end, database_reads::MAXIMUM_EXPECTED_ERRORS))
                    hits[id].insert(start);
            });
        };
        read_linked::Reader first(files[0]);
        std::unique_ptr<read_linked::Reader> second;
        if (files.size() == 2) second = std::make_unique<read_linked::Reader>(files[1]);
        klibpp::KSeq a, b;
        while (first.next(a)) {
            if (second && (!second->next(b) || read_linked::read_name(a.name) != read_linked::read_name(b.name)))
                throw std::runtime_error("Confidence read pairs have different names or lengths");
            std::map<std::size_t, std::set<std::size_t>> hits;
            scan(a, hits); if (second) scan(b, hits);
            ++result.counts["confidence_input_units"];
            if (hits.empty()) continue;
            const std::string unit = a.seq + '\0' + (second ? b.seq : "");
            for (const auto& [id, offsets] : hits) {
                auto& pattern = patterns[id];
                if (pattern.units.size() < SUPPORT_CAP) pattern.units.insert(unit);
                for (auto offset : offsets) {
                    if (pattern.offsets.size() >= SUPPORT_CAP) break;
                    pattern.offsets.insert(offset);
                }
            }
        }
        if (second && second->next(b))
            throw std::runtime_error("Confidence read files contain different record counts");
    }
    auto support = [&](std::size_t id) {
        const auto& p = patterns[id];
        return Support{p.units.size(), p.offsets.size(), p.units.size() >= SUPPORT_CAP,
            p.offsets.size() >= SUPPORT_CAP,
            !files.empty() && max_quality_read_length >= p.sequence.size()};
    };
    for (const auto& p : patterns)
        if (p.units.size() >= SUPPORT_CAP || p.offsets.size() >= SUPPORT_CAP)
            ++result.counts["confidence_capped_patterns"];
    std::map<std::string, std::set<std::string>> repeat_contexts;
    for (std::size_t i = 0; i < pairs.size(); ++i) {
        auto& row = result.rows[i];
        row.repeat = support(row_patterns[i][0]); row.spacer = support(row_patterns[i][1]);
        row.context = support(row_patterns[i][2]);
        row.repeat_score = row.repeat.index(); row.spacer_score = row.spacer.index();
        row.context_score = row.context.index();
        row.context_observability = files.empty() ? "no_reads" :
            max_quality_read_length == 0 ? "no_quality" :
            row.context.observable ? "observable" : "span_exceeds_max_quality_read";
        if (row.context.accepted()) {
            repeat_contexts[canonical_repeats[i]].insert(cores[i]);
            ++result.counts["confidence_context_supported_rows"];
        }
    }
    // Reference copies and multiple categories never accumulate votes. Every
    // score keeps the category, certainty and context count of ONE source array.
    for (const auto& record : records) {
        if (record.confidence < (record.bona_fide ? array_database::MIN_BONA_FIDE_CONFIDENCE :
                array_database::MIN_POSSIBLE_CONFIDENCE)) continue;
        const double weight = record.weight();
        for (const auto& spacer : record.spacers) {
            const auto found = by_spacer.find(spacer);
            if (found != by_spacer.end()) for (auto i : found->second)
                result.rows[i].spacer_database_score = std::max(result.rows[i].spacer_database_score, weight);
        }
        std::set<std::size_t> matched_rows;
        for (const auto& repeat : record.repeats) {
            const auto found = by_repeat.find(repeat);
            if (found == by_repeat.end()) continue;
            for (auto i : found->second) {
                result.rows[i].repeat_database_score = std::max(result.rows[i].repeat_database_score, weight);
                if (result.rows[i].context.accepted() && record.spacers.count(canonical_spacers[i]))
                    matched_rows.insert(i);
            }
        }
        // Collapse duplicated rows and strand copies before leave-one-out counts.
        std::set<std::pair<std::string, std::string>> matched;
        for (auto i : matched_rows) matched.emplace(canonical_spacers[i], cores[i]);
        for (auto i : matched_rows) {
            std::set<std::string> other_spacers, other_cores;
            for (const auto& [spacer, core] : matched)
                if (spacer != canonical_spacers[i] && core != cores[i]) {
                    other_spacers.insert(spacer); other_cores.insert(core);
                }
            const auto count = std::min(other_spacers.size(), other_cores.size());
            const auto minimum = record.bona_fide ? array_database::MIN_BONA_FIDE_SPACERS :
                array_database::MIN_POSSIBLE_SPACERS;
            if (count + 1 < minimum) continue;
            const double score = weight * std::min<std::size_t>(count, 3) / 3.0;
            auto& row = result.rows[i];
            if (score > row.database_score || (score == row.database_score &&
                (row.database_source.empty() || record.id < row.database_source))) {
                row.database_score = score; row.database_source = record.id;
                row.database_other_spacers = other_spacers.size();
                row.database_other_cores = other_cores.size();
            }
        }
    }
    for (std::size_t i = 0; i < pairs.size(); ++i)
        alternatives[cores[i]][identities[i]] = std::max(
            alternatives[cores[i]][identities[i]], result.rows[i].database_score);
    for (std::size_t i = 0; i < pairs.size(); ++i) {
        auto& row = result.rows[i];
        row.repeat_context_core_count = repeat_contexts[canonical_repeats[i]].size();
        if (row.context.observable)
            row.repeat_family_score = std::min<std::size_t>(row.repeat_context_core_count, 4) / 4.0;
        const auto& choices = alternatives.at(cores[i]);
        row.boundary_alternatives = choices.size();
        for (const auto& [identity, score] : choices) if (score > 0) ++row.boundary_database_alternatives;
        if (row.database_score > 0) {
            row.boundary_score = row.database_score / row.boundary_database_alternatives;
            row.boundary_status = row.boundary_database_alternatives == 1 ?
                "catalogue_preferred" : "catalogue_ambiguous";
        } else if (choices.at(identities[i]) > 0) {
            row.boundary_status = "same_cut_without_context_support";
        } else if (row.boundary_database_alternatives > 0) {
            row.boundary_score = 0;
            row.boundary_status = "alternative_to_catalogue_cut";
        } else row.boundary_status = choices.size() > 1 ? "ambiguous_unscored" : "unresolved";
        if (choices.size() > 1) ++result.counts["confidence_ambiguous_rows"];
    }
    return result;
}

inline void write_score(std::ostream& out, double value) {
    if (std::isfinite(value)) out << value; else out << "NA";
}
template<class Arrays>
inline void write_tables(const std::string& directory, const Arrays& arrays, const Result& result) {
    const std::string path = directory + "/element_confidence.tsv";
    std::ofstream out(path);
    if (!out) throw std::runtime_error("Cannot write element confidence: " + path);
    out << "array_id\tentry_index\trepeat_variant\tspacer\trepeat_score\tspacer_score\tcontext_score";
    for (const char* kind : {"repeat", "spacer", "context"})
        for (const char* field : {"sequence_units", "start_offsets", "units_saturated", "offsets_saturated"})
            out << '\t' << kind << '_' << field;
    out << "\tcontext_observability\trepeat_database_score\tspacer_database_score"
        "\trepeat_context_core_count\trepeat_family_score\tdatabase_source\tdatabase_score"
        "\tdatabase_other_spacers\tdatabase_other_cores\tboundary_alternatives"
        "\tboundary_database_alternatives\tboundary_score\tboundary_status\n";
    out << std::setprecision(17);
    std::size_t row_number = 0;
    for (std::size_t a = 0; a < arrays.size(); ++a)
        for (std::size_t e = 0; e < arrays[a].entries.size(); ++e) {
            if (row_number >= result.rows.size()) throw std::runtime_error("Confidence row count mismatch");
            const auto& row = result.rows[row_number++];
            out << "Array_" << a + 1 << '\t' << e + 1 << '\t' << arrays[a].entries[e].first << '\t'
                << arrays[a].entries[e].second << '\t';
            write_score(out, row.repeat_score); out << '\t'; write_score(out, row.spacer_score);
            out << '\t'; write_score(out, row.context_score);
            for (const auto& s : {row.repeat, row.spacer, row.context})
                out << '\t' << s.sequence_units << '\t' << s.start_offsets << '\t'
                    << s.units_saturated << '\t' << s.offsets_saturated;
            out << '\t' << row.context_observability << '\t' << row.repeat_database_score << '\t'
                << row.spacer_database_score << '\t' << row.repeat_context_core_count << '\t';
            write_score(out, row.repeat_family_score);
            out << '\t' << row.database_source << '\t' << row.database_score << '\t'
                << row.database_other_spacers << '\t' << row.database_other_cores << '\t'
                << row.boundary_alternatives << '\t' << row.boundary_database_alternatives << '\t';
            write_score(out, row.boundary_score); out << '\t' << row.boundary_status << '\n';
        }
    if (row_number != result.rows.size()) throw std::runtime_error("Confidence row count mismatch");
    out.close();
    if (!out) throw std::runtime_error("Error writing element confidence: " + path);
}

}} // namespace Zuck::element_confidence
#endif
