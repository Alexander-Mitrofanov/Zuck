#ifndef ZUCK_ARRAY_DATABASE_H
#define ZUCK_ARRAY_DATABASE_H

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <map>
#include <limits>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>
#include "crispr_array/array_database_reads.h"

namespace Zuck { namespace array_database {

// Scores rank evidence; CRISPRidentify certainty scores are not probabilities.
constexpr std::size_t MIN_BONA_FIDE_SPACERS = 2;
constexpr std::size_t MIN_POSSIBLE_SPACERS = 3;
constexpr double MIN_BONA_FIDE_CONFIDENCE = 0.75;
constexpr double MIN_POSSIBLE_CONFIDENCE = 0.5;
constexpr std::size_t MAX_CORE_INTERPRETATIONS = 128;
constexpr std::size_t MIN_REPEAT = 20, MAX_REPEAT = 55;
constexpr std::size_t MIN_SPACER = 20, MAX_SPACER = 70;

inline std::string reverse_complement(const std::string& s) {
    std::string r(s.size(), 'N');
    for (std::size_t i = 0; i < s.size(); ++i) {
        const char b = s[i];
        r[s.size() - 1 - i] = b == 'A' ? 'T' : b == 'T' ? 'A' :
                                b == 'C' ? 'G' : b == 'G' ? 'C' : 'N';
    }
    return r;
}
inline std::string canonical(const std::string& s) {
    return std::min(s, reverse_complement(s));
}
inline bool dna(const std::string& s) {
    return !s.empty() && s.find_first_not_of("ACGT") == std::string::npos;
}
inline std::vector<std::string> split(const std::string& s, char separator) {
    std::vector<std::string> out;
    std::size_t start = 0;
    for (;;) {
        const auto end = s.find(separator, start);
        out.push_back(s.substr(start, end == std::string::npos ? end : end - start));
        if (end == std::string::npos) return out;
        start = end + 1;
    }
}
inline std::string circular_key(const std::string& s) {
    std::string best = s;
    for (std::size_t p = 1; p < s.size(); ++p) {
        const std::string rotated = s.substr(p) + s.substr(0, p);
        if (rotated < best) best = rotated;
    }
    return best;
}
inline std::string signature(const std::set<std::string>& repeats,
                             const std::set<std::string>& spacers) {
    std::string key;
    for (const auto& r : repeats) key += r + ',';
    key += '|';
    for (const auto& s : spacers) key += s + ',';
    return key;
}

struct Record {
    std::string id, category;
    double confidence = 0;
    bool bona_fide = false;
    std::set<std::string> repeats, spacers;
    double weight() const { return (bona_fide ? 1.0 : 0.5) * confidence; }
};

struct Evidence {
    std::string repeat, spacer, source_array, category;
    double confidence = 0, score = 0;
    std::size_t distinct_repeat_matches = 0, distinct_spacer_matches = 0;
    std::size_t observed_circular_cores = 0, cycle_records = 0;
    bool novel_spacer = false, read_validation_applied = false;
    std::size_t read_sequence_units = 0, read_start_offsets = 0;
};

struct Group {
    std::string consensus;
    std::vector<std::pair<std::string, std::string>> entries;
    std::vector<Evidence> evidence;
};

struct Result {
    std::vector<Group> groups;
    std::map<std::string, std::size_t> counts;
};

class Database {
    std::vector<Record> records_;
    std::unordered_map<std::string, std::vector<std::size_t>> repeat_records_;
    std::unordered_map<std::string, std::vector<std::string>> repeat_seeds_;
    std::size_t source_rows_ = 0;

    static std::set<std::string> sequences(const std::string& value,
                                         std::size_t minimum, std::size_t maximum) {
        std::set<std::string> out;
        for (const auto& s : split(value, ',')) {
            if (!dna(s) || s.size() < minimum || s.size() > maximum)
                throw std::runtime_error("Invalid sequence in array database");
            out.insert(canonical(s));
        }
        return out;
    }
    static std::size_t count(const std::string& s) {
        std::size_t consumed = 0;
        if (s.empty() || s.find_first_not_of("0123456789") != std::string::npos)
            throw std::runtime_error("Invalid database count");
        const auto n = std::stoull(s, &consumed);
        if (consumed != s.size() || n > std::numeric_limits<std::size_t>::max())
            throw std::runtime_error("Invalid database count");
        return n;
    }

public:
    explicit Database(const std::string& path, bool include_possible = true) {
        std::ifstream input(path);
        std::string line;
        if (!input || !std::getline(input, line))
            throw std::runtime_error("Cannot read array database: " + path);
        if (!line.empty() && line.back() == '\r') line.pop_back();
        if (line != "array_id\taccession\tcategory\tconfidence\trepeat_count\tspacer_count\trepeats\tspacers")
            throw std::runtime_error("Invalid array database header: " + path);
        std::map<std::string, std::size_t> signatures;
        std::set<std::string> ids;
        std::size_t number = 1;
        while (std::getline(input, line)) {
            ++number;
            if (!line.empty() && line.back() == '\r') line.pop_back();
            try {
                const auto fields = split(line, '\t');
                if (fields.size() != 8 || fields[0].empty() || fields[1].empty())
                    throw std::runtime_error("Expected eight populated database fields");
                if (!ids.insert(fields[0]).second)
                    throw std::runtime_error("Duplicate database array id");
                Record record;
                record.id = fields[0];
                record.category = fields[2];
                record.bona_fide = record.category == "Bona-fide" ||
                                   record.category == "Bona-Fide" || record.category == "Bona_fide";
                if (!record.bona_fide && record.category != "Possible")
                    throw std::runtime_error("Unsupported database category");
                std::size_t consumed = 0;
                record.confidence = std::stod(fields[3], &consumed);
                if (consumed != fields[3].size() || !std::isfinite(record.confidence) ||
                    record.confidence < 0 || record.confidence > 1)
                    throw std::runtime_error("Database confidence must be finite in [0,1]");
                const auto nr = count(fields[4]), ns = count(fields[5]);
                record.repeats = sequences(fields[6], MIN_REPEAT, MAX_REPEAT);
                record.spacers = sequences(fields[7], MIN_SPACER, MAX_SPACER);
                if (nr < 2 || ns < 1)
                    throw std::runtime_error("Database array lacks repeat/spacer geometry");
                ++source_rows_;
                if (!include_possible && !record.bona_fide) continue;
                // An ineligible record cannot suppress an eligible duplicate
                // from another confidence category during deduplication.
                if (record.confidence < (record.bona_fide ? MIN_BONA_FIDE_CONFIDENCE : MIN_POSSIBLE_CONFIDENCE))
                    continue;
                const auto key = signature(record.repeats, record.spacers);
                const auto found = signatures.find(key);
                if (found == signatures.end()) {
                    signatures[key] = records_.size();
                    records_.push_back(std::move(record));
                } else {
                    auto& old = records_[found->second];
                    // Duplicate assemblies do not accumulate votes. Select one
                    // complete provenance record, never mix category and score.
                    if (record.weight() > old.weight() ||
                        (record.weight() == old.weight() && record.id < old.id))
                        old = std::move(record);
                }
            } catch (const std::exception& e) {
                throw std::runtime_error(path + ":" + std::to_string(number) + ": " + e.what());
            }
        }
        if (!input.eof()) throw std::runtime_error("Error reading array database: " + path);
        std::set<std::string> oriented_repeats;
        for (std::size_t i = 0; i < records_.size(); ++i)
            for (const auto& r : records_[i].repeats) {
                repeat_records_[r].push_back(i);
                oriented_repeats.insert(r);
                oriented_repeats.insert(reverse_complement(r));
            }
        for (const auto& r : oriented_repeats) repeat_seeds_[r.substr(0, MIN_REPEAT)].push_back(r);
    }

    const std::vector<Record>& records() const { return records_; }

    Result infer(const std::vector<std::string>& cycles,
                 const std::vector<std::string>& reads = {},
                 bool validate_reads = false, bool novel_spacers = false) const {
        if (novel_spacers && !validate_reads)
            throw std::runtime_error("Novel database spacers require read validation");
        Result result;
        result.counts["database_source_rows"] = source_rows_;
        result.counts["database_distinct_array_signatures"] = records_.size();
        result.counts["database_distinct_repeats"] = repeat_records_.size();
        struct Candidate {
            std::string repeat, spacer;
            Evidence evidence;
            bool eligible = false;
            bool read_supported = true;
            database_reads::Support read_support;
        };
        std::map<std::string, std::size_t> unique_cores;
        for (const auto& core : cycles)
            if (dna(core) && core.size() >= MIN_REPEAT + MIN_SPACER &&
                core.size() <= MAX_REPEAT + MAX_SPACER)
                ++unique_cores[circular_key(core)];
        std::vector<std::size_t> core_multiplicity;
        std::vector<std::vector<std::size_t>> core_candidates;
        std::map<std::pair<std::string, std::string>, std::size_t> by_pair;
        std::vector<Candidate> candidates;
        for (const auto& [core, multiplicity] : unique_cores) {
            core_multiplicity.push_back(multiplicity);
            std::set<std::pair<std::string, std::string>> pairs, novel_pairs;
            const std::string doubled = core + core;
            for (std::size_t pos = 0; pos < core.size(); ++pos) {
                const auto found = repeat_seeds_.find(doubled.substr(pos, MIN_REPEAT));
                if (found == repeat_seeds_.end()) continue;
                for (const auto& repeat : found->second) {
                    if (repeat.size() + MIN_SPACER > core.size() ||
                        core.size() - repeat.size() > MAX_SPACER ||
                        doubled.compare(pos, repeat.size(), repeat) != 0) continue;
                    const auto spacer = doubled.substr(pos + repeat.size(), core.size() - repeat.size());
                    const auto canonical_spacer = canonical(spacer);
                    bool known = false;
                    for (auto record : repeat_records_.at(canonical(repeat)))
                        if (records_[record].spacers.count(canonical_spacer)) { known = true; break; }
                    if (known) pairs.emplace(repeat, spacer);
                    else if (novel_spacers) novel_pairs.emplace(repeat, spacer);
                }
            }
            core_candidates.emplace_back();
            if (pairs.size() > MAX_CORE_INTERPRETATIONS) {
                ++result.counts["database_overloaded_cores"];
                continue;
            }
            if (pairs.size() + novel_pairs.size() <= MAX_CORE_INTERPRETATIONS)
                pairs.insert(novel_pairs.begin(), novel_pairs.end());
            else if (!novel_pairs.empty()) ++result.counts["database_overloaded_novel_cores"];
            for (const auto& pair : pairs) {
                auto entry = by_pair.emplace(pair, candidates.size());
                if (entry.second) {
                    Candidate c;
                    c.repeat = pair.first; c.spacer = pair.second;
                    candidates.push_back(std::move(c));
                }
                const auto id = entry.first->second;
                core_candidates.back().push_back(id);
            }
        }
        result.counts["database_observed_pair_candidates"] = candidates.size();
        if (validate_reads) {
            std::vector<std::pair<std::string, std::string>> pairs;
            for (const auto& candidate : candidates) pairs.emplace_back(candidate.repeat, candidate.spacer);
            const auto validated = database_reads::validate(pairs, reads);
            result.counts.insert(validated.counts.begin(), validated.counts.end());
            for (std::size_t i = 0; i < candidates.size(); ++i) {
                candidates[i].read_support = validated.support[i];
                candidates[i].read_supported = validated.support[i].accepted();
            }
        } else result.counts["database_read_validation_applied"] = 0;
        std::vector<std::vector<std::size_t>> record_candidates(records_.size());
        std::vector<std::vector<std::size_t>> novel_record_candidates(records_.size());
        for (std::size_t i = 0; i < candidates.size(); ++i) {
            if (!candidates[i].read_supported) continue;
            const auto canonical_spacer = canonical(candidates[i].spacer);
            for (auto record : repeat_records_.at(canonical(candidates[i].repeat))) {
                if (records_[record].spacers.count(canonical_spacer))
                    record_candidates[record].push_back(i);
                else if (novel_spacers) novel_record_candidates[record].push_back(i);
            }
        }
        for (std::size_t a = 0; a < records_.size(); ++a) {
            const auto& record = records_[a];
            std::set<std::string> matched_repeats, matched_spacers;
            for (auto c : record_candidates[a]) {
                matched_repeats.insert(canonical(candidates[c].repeat));
                matched_spacers.insert(canonical(candidates[c].spacer));
            }
            if (matched_spacers.size() < (record.bona_fide ? MIN_BONA_FIDE_SPACERS : MIN_POSSIBLE_SPACERS) ||
                record.confidence < (record.bona_fide ? MIN_BONA_FIDE_CONFIDENCE : MIN_POSSIBLE_CONFIDENCE)) continue;
            ++result.counts["database_supported_array_signatures"];
            // Two units of spacer diversity, one of repeat diversity. Saturate
            // at four spacer types and two repeat alleles; duplicate reference
            // copies and repeated graph traversals contribute no votes.
            const double score = record.weight() *
                (2.0 * std::min<std::size_t>(4, matched_spacers.size()) / 4.0 +
                 std::min<std::size_t>(2, matched_repeats.size()) / 2.0);
            auto offer = [&](std::size_t c, bool novel) {
                auto& candidate = candidates[c];
                // An exact same-array spacer match always outranks a novel
                // spacer supported only by other known family members.
                const bool same_kind = novel == candidate.evidence.novel_spacer;
                if (!candidate.eligible || (!novel && candidate.evidence.novel_spacer) ||
                    (same_kind && (score > candidate.evidence.score ||
                     (score == candidate.evidence.score && record.id < candidate.evidence.source_array)))) {
                    candidate.eligible = true;
                    candidate.evidence = {candidate.repeat, candidate.spacer, record.id, record.category,
                        record.confidence, score, matched_repeats.size(), matched_spacers.size(), 0, 0,
                        novel, validate_reads, candidate.read_support.sequence_units,
                        candidate.read_support.start_offsets};
                }
            };
            for (auto c : record_candidates[a]) offer(c, false);
            for (auto c : novel_record_candidates[a]) offer(c, true);
        }
        std::map<std::size_t, std::set<std::size_t>> selected;
        for (std::size_t core = 0; core < core_candidates.size(); ++core) {
            std::size_t best = candidates.size();
            double score = 0;
            std::set<std::pair<std::string, std::string>> tied;
            bool known_eligible = false;
            for (auto c : core_candidates[core])
                if (candidates[c].eligible && !candidates[c].evidence.novel_spacer) known_eligible = true;
            for (auto c : core_candidates[core]) {
                const auto& candidate = candidates[c];
                if (known_eligible && candidate.evidence.novel_spacer) continue;
                if (!candidate.eligible) continue;
                const auto identity = std::make_pair(canonical(candidate.repeat), canonical(candidate.spacer));
                if (candidate.evidence.score > score + 1e-12) {
                    score = candidate.evidence.score;
                    best = c;
                    tied.clear(); tied.insert(identity);
                } else if (std::abs(candidate.evidence.score - score) <= 1e-12) {
                    tied.insert(identity);
                }
            }
            if (best == candidates.size()) continue;
            if (tied.size() > 1) { ++result.counts["database_ambiguous_cores"]; continue; }
            selected[best].insert(core);
        }
        std::map<std::string, Group> groups;
        for (const auto& [id, cores] : selected) {
            auto evidence = candidates[id].evidence;
            evidence.observed_circular_cores = cores.size();
            for (auto core : cores) evidence.cycle_records += core_multiplicity[core];
            auto& group = groups[canonical(evidence.repeat)];
            group.consensus = canonical(evidence.repeat);
            group.entries.emplace_back(evidence.repeat, evidence.spacer);
            group.evidence.push_back(std::move(evidence));
        }
        for (auto& [repeat, group] : groups) {
            std::set<std::string> spacers;
            for (const auto& entry : group.entries) spacers.insert(canonical(entry.second));
            if (spacers.size() < 2) continue;
            result.groups.push_back(std::move(group));
        }
        return result;
    }
};

template<class Arrays>
inline void append_groups(Arrays& arrays, const Result& result,
                          const std::string& evidence_path,
                          std::map<std::string, std::size_t>& counts) {
    counts.insert(result.counts.begin(), result.counts.end());
    counts["database_applied"] = 1;
    counts["database_de_novo_groups"] = arrays.size();
    counts["database_groups_added"] = 0;
    counts["database_rows_added"] = 0;
    std::set<std::pair<std::string, std::string>> old_pairs;
    std::set<std::string> old_spacers, old_group_signatures;
    for (const auto& a : arrays) {
        std::set<std::string> repeats, spacers;
        for (const auto& row : a.entries) {
            old_pairs.insert(row);
            old_spacers.insert(canonical(row.second));
            repeats.insert(canonical(row.first)); spacers.insert(canonical(row.second));
        }
        old_group_signatures.insert(signature(repeats, spacers));
    }
    std::ofstream out(evidence_path);
    if (!out) throw std::runtime_error("Cannot write database evidence: " + evidence_path);
    out << "array_id\tentry_index\trepeat_variant\tspacer\tsource_array\tcategory\tconfidence\tevidence_score"
           "\tdistinct_repeat_matches\tdistinct_spacer_matches\tobserved_circular_cores\tcycle_records"
           "\tnew_oriented_pair\tnew_canonical_spacer\tevidence_kind\tread_validation_applied"
           "\tread_distinct_sequence_units\tread_distinct_start_offsets\tread_span_rule\n";
    out << std::setprecision(17);
    std::set<std::string> new_spacers;
    std::set<std::pair<std::string, std::string>> new_pairs;
    for (const auto& group : result.groups) {
        std::set<std::string> repeats, spacers;
        for (const auto& row : group.entries) {
            repeats.insert(canonical(row.first)); spacers.insert(canonical(row.second));
        }
        if (old_group_signatures.count(signature(repeats, spacers))) continue;
        arrays.push_back({group.consensus, group.entries});
        ++counts["database_groups_added"];
        counts["database_rows_added"] += group.entries.size();
        for (std::size_t i = 0; i < group.evidence.size(); ++i) {
            const auto& e = group.evidence[i];
            const bool new_pair = !old_pairs.count({e.repeat, e.spacer});
            const bool new_spacer = !old_spacers.count(canonical(e.spacer));
            if (new_pair) new_pairs.emplace(e.repeat, e.spacer);
            if (new_spacer) new_spacers.insert(canonical(e.spacer));
            out << "Array_" << arrays.size() << '\t' << i + 1 << '\t' << e.repeat << '\t'
                << e.spacer << '\t' << e.source_array << '\t' << e.category << '\t'
                << e.confidence << '\t' << e.score << '\t' << e.distinct_repeat_matches << '\t'
                << e.distinct_spacer_matches << '\t' << e.observed_circular_cores << '\t'
                << e.cycle_records << '\t' << new_pair << '\t' << new_spacer << '\t'
                << (e.novel_spacer ? "novel_spacer_trusted_repeat" : "known_same_array") << '\t'
                << e.read_validation_applied << '\t' << e.read_sequence_units << '\t'
                << e.read_start_offsets << '\t'
                << (e.read_validation_applied ? "whole_R_S_nextR16_q10_ee1" : "not_applied") << '\n';
        }
    }
    counts["database_new_oriented_pairs"] = new_pairs.size();
    counts["database_new_canonical_spacers"] = new_spacers.size();
    out.close();
    if (!out) throw std::runtime_error("Error writing database evidence: " + evidence_path);
}
}} // namespace Zuck::array_database
#endif
