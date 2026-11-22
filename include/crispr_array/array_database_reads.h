#ifndef ZUCK_ARRAY_DATABASE_READS_H
#define ZUCK_ARRAY_DATABASE_READS_H

#include "crispr_array/read_linked_groups.h"

namespace Zuck { namespace database_reads {

constexpr std::size_t NEXT_REPEAT_ANCHOR = 16;
constexpr int MINIMUM_QUALITY = 10;
constexpr double MAXIMUM_EXPECTED_ERRORS = 1.0;
constexpr std::size_t MINIMUM_SEQUENCE_UNITS = 2;
constexpr std::size_t MINIMUM_START_OFFSETS = 2;

struct Support {
    std::size_t sequence_units = 0, start_offsets = 0;
    bool accepted() const {
        return sequence_units >= MINIMUM_SEQUENCE_UNITS && start_offsets >= MINIMUM_START_OFFSETS;
    }
};
struct Result {
    std::vector<Support> support;
    std::map<std::string, std::size_t> counts;
};

// A circular graph path does not prove its emitted linear adjacency. Require
// the literal full repeat + full spacer + next-repeat prefix in one read.
// Reverse complement the entire span, never its constituent parts separately.
inline Result validate(const std::vector<std::pair<std::string, std::string>>& pairs,
                       const std::vector<std::string>& files) {
    if (files.empty() || files.size() > 2)
        throw std::runtime_error("Database read validation requires one or two read files");
    Result result;
    result.support.resize(pairs.size());
    result.counts = {{"database_read_validation_applied", 1},
        {"database_read_records", 0}, {"database_read_records_without_quality", 0},
        {"database_read_input_units", 0}, {"database_read_distinct_support_sequence_units", 0},
        {"database_read_supported_pair_candidates", 0}};
    read_evidence::ExactMatcher matcher;
    std::vector<std::size_t> lengths;
    for (std::size_t i = 0; i < pairs.size(); ++i) {
        const auto& [repeat, spacer] = pairs[i];
        if (repeat.size() < NEXT_REPEAT_ANCHOR || spacer.empty())
            throw std::runtime_error("Invalid database read span");
        const auto pattern = repeat + spacer + repeat.substr(0, NEXT_REPEAT_ANCHOR);
        lengths.push_back(pattern.size());
        matcher.add(pattern, i);
        const auto reverse = read_linked::reverse_complement(pattern);
        if (reverse != pattern) matcher.add(reverse, i);
    }
    matcher.build();
    std::vector<read_linked::Support> support(pairs.size());
    std::map<std::pair<std::string, std::string>, std::size_t> units;
    auto scan = [&](const klibpp::KSeq& read, std::map<std::size_t, std::set<std::size_t>>& hits) {
        ++result.counts["database_read_records"];
        if (read.qual.empty()) { ++result.counts["database_read_records_without_quality"]; return; }
        read_evidence::QualitySpans quality(read.qual, MINIMUM_QUALITY);
        matcher.scan(read.seq, [&](std::size_t id, std::size_t end) {
            const auto start = end - lengths[id];
            if (quality.passes(start, end, MAXIMUM_EXPECTED_ERRORS))
                hits[id].insert(start);
        });
    };
    read_linked::Reader first(files[0]);
    std::unique_ptr<read_linked::Reader> second;
    if (files.size() == 2) second = std::make_unique<read_linked::Reader>(files[1]);
    klibpp::KSeq a, b;
    while (first.next(a)) {
        if (second && (!second->next(b) || read_linked::read_name(a.name) != read_linked::read_name(b.name)))
            throw std::runtime_error("Database validation read pairs have different names or lengths");
        std::map<std::size_t, std::set<std::size_t>> hits;
        scan(a, hits); if (second) scan(b, hits);
        ++result.counts["database_read_input_units"];
        if (hits.empty()) continue;
        const auto unit = units.emplace(std::make_pair(a.seq, second ? b.seq : ""), units.size());
        for (const auto& [id, offsets] : hits) {
            support[id].sequence_units.insert(unit.first->second);
            support[id].start_offsets.insert(offsets.begin(), offsets.end());
        }
    }
    if (second && second->next(b))
        throw std::runtime_error("Database validation read files contain different record counts");
    result.counts["database_read_distinct_support_sequence_units"] = units.size();
    for (std::size_t i = 0; i < support.size(); ++i) {
        result.support[i] = {support[i].sequence_units.size(), support[i].start_offsets.size()};
        if (result.support[i].accepted()) ++result.counts["database_read_supported_pair_candidates"];
    }
    return result;
}
}} // namespace Zuck::database_reads
#endif
