#ifndef ZUCK_SHORT_SPACER_FILTER_H
#define ZUCK_SHORT_SPACER_FILTER_H

#include <cstddef>
#include <string>
#include <utility>
#include <vector>

namespace Zuck {

inline constexpr std::size_t SHORT_SPACER_ANCHOR_LENGTH = 23;
inline constexpr std::size_t SHORT_SPACER_MAX_TERMINAL_EXTENSION = 2;

// The 23-mer collapse cannot connect shorter terminal variants. Reject only
// singleton evidence under one observed longest anchor; never chain overlaps,
// change an emitted sequence, or infer that a short spacer is itself invalid.
inline bool has_only_short_terminal_variants(
    const std::vector<std::pair<std::string, std::string>>& entries) {
    if (entries.size() < 2) return false;
    const std::string* anchor = &entries.front().second;
    for (const auto& entry : entries)
        if (entry.second.size() > anchor->size() ||
            (entry.second.size() == anchor->size() && entry.second < *anchor))
            anchor = &entry.second;

    bool has_shorter_variant = false;
    for (const auto& entry : entries) {
        const auto& spacer = entry.second;
        if (spacer == *anchor) continue;
        if (spacer.size() >= SHORT_SPACER_ANCHOR_LENGTH ||
            anchor->size() - spacer.size() > SHORT_SPACER_MAX_TERMINAL_EXTENSION)
            return false;
        const bool prefix = anchor->compare(0, spacer.size(), spacer) == 0;
        const bool suffix = anchor->compare(anchor->size() - spacer.size(),
                                            spacer.size(), spacer) == 0;
        if (!prefix && !suffix) return false;
        has_shorter_variant = true;
    }
    return has_shorter_variant;
}

}  // namespace Zuck

#endif
