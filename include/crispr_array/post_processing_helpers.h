#ifndef ZUCK_POST_PROCESSING_HELPERS_H
#define ZUCK_POST_PROCESSING_HELPERS_H

#include <algorithm>
#include <cctype>
#include <limits>
#include <map>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <unordered_map>
#include <vector>

namespace Zuck::post_processing {

// Replay accepts either one compact sequence or legacy whitespace-separated
// k-mers. Normalize before overlap matching so mixed-case DNA means the same
// thing as uppercase DNA. Non-overlapping tokens keep the legacy concatenation.
inline std::string glue_kmers(const std::string& line) {
    std::istringstream input(line);
    std::string result;
    std::string next;
    while (input >> next) {
        for (char& base : next) {
            base = static_cast<char>(std::toupper(static_cast<unsigned char>(base)));
            if (base != 'A' && base != 'C' && base != 'G' && base != 'T' && base != 'N')
                throw std::invalid_argument("Cycle records must contain only A, C, G, T or N");
        }
        std::size_t overlap = std::min(result.size(), next.size());
        while (overlap > 0 && result.compare(result.size() - overlap, overlap, next, 0, overlap) != 0)
            --overlap;
        result.append(next, overlap, std::string::npos);
    }
    return result;
}

// Exact Levenshtein distance up to limit; return limit + 1 beyond it.
// Only the diagonal band can contain an acceptable alignment. Its bounded
// cells reduce the work for the repeat-family radius (normally four edits).
inline int bounded_edit_distance(std::string_view a, std::string_view b, int limit) {
    if (limit < 0 || limit == std::numeric_limits<int>::max())
        throw std::invalid_argument("Edit-distance limit must be nonnegative and below INT_MAX");
    if (a.size() < b.size()) std::swap(a, b);
    const auto radius = static_cast<std::size_t>(limit);
    const int beyond = limit + 1;
    const auto increment = [beyond](int value) { return value < beyond ? value + 1 : beyond; };
    if (a.size() - b.size() > radius) return beyond;
    if (a == b) return 0;
    if (b.empty()) return static_cast<int>(a.size());

    std::vector<int> previous(b.size() + 1, beyond);
    std::vector<int> current(b.size() + 1, beyond);
    for (std::size_t column = 0; column <= std::min(b.size(), radius); ++column)
        previous[column] = static_cast<int>(column);
    for (std::size_t row = 1; row <= a.size(); ++row) {
        const std::size_t begin = row > radius ? row - radius : 1;
        const std::size_t end = std::min(b.size(), row + radius);
        current[0] = row <= radius ? static_cast<int>(row) : beyond;
        if (begin > 1) current[begin - 1] = beyond;
        int minimum = current[0];
        for (std::size_t column = begin; column <= end; ++column) {
            const int substitution = a[row - 1] == b[column - 1]
                ? previous[column - 1] : increment(previous[column - 1]);
            current[column] = std::min({increment(previous[column]),
                                       increment(current[column - 1]), substitution});
            minimum = std::min(minimum, current[column]);
        }
        if (end < b.size()) current[end + 1] = beyond;
        if (minimum > limit) return beyond;
        previous.swap(current);
    }
    return previous.back();
}

// Preserve the historical rule: any shared exact k-mer joins a transitive
// component. A single owner per k-mer is enough; storing all pairwise links is
// unnecessary. Indices within each component retain their input order.
inline std::vector<std::vector<int>> shared_kmer_components(
    const std::vector<std::string>& alleles, std::size_t kmer_size) {
    if (kmer_size == 0) throw std::invalid_argument("K-mer length must be positive");
    if (alleles.size() > static_cast<std::size_t>(std::numeric_limits<int>::max()))
        throw std::length_error("Too many spacer alleles");
    std::vector<int> parent(alleles.size());
    std::vector<std::size_t> size(alleles.size(), 1);
    std::iota(parent.begin(), parent.end(), 0);
    const auto find = [&](int node) {
        while (parent[node] != node) {
            parent[node] = parent[parent[node]];
            node = parent[node];
        }
        return node;
    };
    std::unordered_map<std::string_view, int> owner;
    for (int index = 0; index < static_cast<int>(alleles.size()); ++index) {
        const std::string_view sequence(alleles[index]);
        if (sequence.size() < kmer_size) continue;
        for (std::size_t position = 0; position <= sequence.size() - kmer_size; ++position) {
            const auto [entry, inserted] = owner.emplace(sequence.substr(position, kmer_size), index);
            if (inserted) continue;
            int left = find(index), right = find(entry->second);
            if (left == right) continue;
            if (size[left] < size[right]) std::swap(left, right);
            parent[right] = left;
            size[left] += size[right];
        }
    }
    std::map<int, std::vector<int>> by_root;
    for (int index = 0; index < static_cast<int>(alleles.size()); ++index)
        by_root[find(index)].push_back(index);
    std::vector<std::vector<int>> components;
    components.reserve(by_root.size());
    for (auto& [root, members] : by_root) components.push_back(std::move(members));
    std::sort(components.begin(), components.end(), [&](const auto& left, const auto& right) {
        return alleles[left.front()] < alleles[right.front()];
    });
    return components;
}

}  // namespace Zuck::post_processing
#endif
