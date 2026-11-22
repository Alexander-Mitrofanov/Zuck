#ifndef ZUCK_DATABASE_GROUPS_H
#define ZUCK_DATABASE_GROUPS_H

#include "crispr_array/array_database.h"
#include <filesystem>

namespace Zuck { namespace database_groups {

constexpr std::size_t MIN_BONA_FIDE_SPACERS = 3;
constexpr std::size_t MIN_POSSIBLE_SPACERS = 4;
constexpr std::size_t MAX_VIEWS_PER_PARENT = 16;

struct View {
    std::size_t parent = 0;
    std::vector<std::size_t> rows;  // One-based rows in the accepted parent.
    std::set<std::string> spacers;
    std::string source_array, category;
    double confidence = 0, weight = 0;
    std::size_t source_spacers = 0;
};
struct Result {
    std::vector<View> views;
    std::map<std::string, std::size_t> counts;
};

// Reference-assisted content views over accepted, read-corroborated DB rows.
// No source record can add a sequence or remove an original parent row.
template<class Arrays>
Result infer(const Arrays& arrays, std::size_t de_novo_count,
             const std::vector<array_database::Record>& records,
             bool include_possible = false) {
    Result result;
    std::map<std::string, std::vector<std::size_t>> repeat_parents;
    std::vector<std::set<std::string>> parent_spacers(arrays.size());
    for (std::size_t p = de_novo_count; p < arrays.size(); ++p) {
        std::set<std::string> repeats;
        for (const auto& row : arrays[p].entries) {
            repeats.insert(array_database::canonical(row.first));
            parent_spacers[p].insert(array_database::canonical(row.second));
        }
        for (const auto& repeat : repeats) repeat_parents[repeat].push_back(p);
    }
    using Content = std::set<std::string>;
    std::map<std::size_t, std::map<Content, View>> candidates;
    for (const auto& record : records) {
        if (!include_possible && !record.bona_fide) continue;
        const auto floor = record.bona_fide ? array_database::MIN_BONA_FIDE_CONFIDENCE
                                          : array_database::MIN_POSSIBLE_CONFIDENCE;
        if (record.confidence < floor) continue;
        std::set<std::size_t> parents;
        for (const auto& repeat : record.repeats) {
            const auto found = repeat_parents.find(repeat);
            if (found != repeat_parents.end()) parents.insert(found->second.begin(), found->second.end());
        }
        for (auto parent : parents) {
            View view;
            view.parent = parent; view.source_array = record.id; view.category = record.category;
            view.confidence = record.confidence; view.weight = record.weight();
            view.source_spacers = record.spacers.size();
            const auto& entries = arrays[parent].entries;
            for (std::size_t i = 0; i < entries.size(); ++i) {
                const auto& row = entries[i];
                const auto spacer = array_database::canonical(row.second);
                if (record.repeats.count(array_database::canonical(row.first)) && record.spacers.count(spacer)) {
                    view.rows.push_back(i + 1);
                    view.spacers.insert(spacer);
                }
            }
            if (view.spacers.size() < (record.bona_fide ? MIN_BONA_FIDE_SPACERS : MIN_POSSIBLE_SPACERS)) continue;
            ++result.counts["database_group_supported_source_intersections"];
            auto& by_content = candidates[parent];
            const auto old = by_content.find(view.spacers);
            if (old == by_content.end() || view.weight > old->second.weight ||
                (view.weight == old->second.weight && view.source_array < old->second.source_array)) {
                const auto content = view.spacers;
                by_content[content] = std::move(view);
            }
        }
    }
    for (auto& [parent, by_content] : candidates) {
        result.counts["database_group_distinct_intersections"] += by_content.size();
        std::vector<View> ordered;
        for (auto& [content, view] : by_content) ordered.push_back(std::move(view));
        std::sort(ordered.begin(), ordered.end(), [](const View& a, const View& b) {
            if (a.spacers.size() != b.spacers.size()) return a.spacers.size() > b.spacers.size();
            return a.spacers < b.spacers;
        });
        std::vector<View> maximal;
        for (auto& view : ordered) {
            const bool subset = std::any_of(maximal.begin(), maximal.end(), [&](const View& other) {
                return std::includes(other.spacers.begin(), other.spacers.end(), view.spacers.begin(), view.spacers.end());
            });
            if (subset) ++result.counts["database_group_nested_intersections_removed"];
            else maximal.push_back(std::move(view));
        }
        if (maximal.size() > MAX_VIEWS_PER_PARENT) {
            ++result.counts["database_group_ambiguous_parents_abstained"];
            result.counts["database_group_ambiguous_views_abstained"] += maximal.size();
            continue;
        }
        for (auto& view : maximal) {
            if (view.spacers == parent_spacers[parent]) {
                ++result.counts["database_group_full_parent_duplicates_removed"];
                continue;
            }
            result.counts["database_group_memberships"] += view.rows.size();
            result.views.push_back(std::move(view));
        }
    }
    result.counts["database_group_views"] = result.views.size();
    return result;
}

inline void write_tables(const std::string& directory, const Result& result) {
    const auto metadata_path = std::filesystem::path(directory) / "database_groups.tsv";
    const auto members_path = std::filesystem::path(directory) / "database_group_members.tsv";
    std::ofstream metadata(metadata_path), members(members_path);
    if (!metadata || !members) throw std::runtime_error("Cannot write database grouping tables");
    metadata << "group_id\tparent_array_id\tsource_array\tcategory\tconfidence\tmatched_spacers\tsource_spacers\tsource_coverage\n";
    members << "group_id\tparent_array_id\tentry_index\n";
    metadata << std::setprecision(17);
    for (std::size_t i = 0; i < result.views.size(); ++i) {
        const auto& view = result.views[i];
        const auto group_id = "DatabaseGroup_" + std::to_string(i + 1);
        const auto parent_id = "Array_" + std::to_string(view.parent + 1);
        metadata << group_id << '\t' << parent_id << '\t' << view.source_array << '\t' << view.category << '\t'
                 << view.confidence << '\t' << view.spacers.size() << '\t' << view.source_spacers << '\t'
                 << static_cast<double>(view.spacers.size()) / view.source_spacers << '\n';
        for (auto row : view.rows) members << group_id << '\t' << parent_id << '\t' << row << '\n';
    }
    metadata.close(); members.close();
    if (!metadata || !members) throw std::runtime_error("Error writing database grouping tables");
}

}} // namespace Zuck::database_groups
#endif
