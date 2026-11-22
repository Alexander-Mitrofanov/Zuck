#include "crispr_array/evaluation.h"

#include <algorithm>
#include <numeric>

std::size_t get_levenshtein_distance(const std::string& a, const std::string& b) {
    if (a == b) return 0;
    const auto& shorter = a.size() <= b.size() ? a : b;
    const auto& longer = a.size() <= b.size() ? b : a;
    if (shorter.empty()) return longer.size();

    std::vector<std::size_t> previous(shorter.size() + 1);
    std::vector<std::size_t> current(shorter.size() + 1);
    std::iota(previous.begin(), previous.end(), 0);
    for (std::size_t row = 1; row <= longer.size(); ++row) {
        current[0] = row;
        for (std::size_t column = 1; column <= shorter.size(); ++column) {
            const std::size_t substitution =
                previous[column - 1] + (longer[row - 1] != shorter[column - 1]);
            current[column] = std::min({previous[column] + 1,
                                        current[column - 1] + 1, substitution});
        }
        previous.swap(current);
    }
    return previous.back();
}

float get_string_similarity(const std::string& a, const std::string& b) {
    const auto length = std::max(a.size(), b.size());
    if (length == 0) return 1.0f;
    return 1.0f - static_cast<float>(get_levenshtein_distance(a, b)) /
                      static_cast<float>(length);
}

int get_number_of_duplicate_spacers(
    const std::vector<std::string>& spacers,
    const std::string& expected_sequence) {
    int duplicates = 0;
    for (const auto& spacer : spacers) {
        if (spacer.empty()) continue;
        int occurrences = 0;
        std::size_t position = 0;
        while ((position = expected_sequence.find(spacer, position)) != std::string::npos) {
            ++occurrences;
            ++position;  // Retain overlapping matches.
        }
        if (occurrences > 1) duplicates += occurrences - 1;
    }
    return duplicates;
}

std::string get_most_similar_sequence(
    const std::string& sequence,
    const std::vector<std::string>& choices) {
    const std::string* best = nullptr;
    float best_similarity = -1.0f;
    for (const auto& choice : choices) {
        const float similarity = get_string_similarity(sequence, choice);
        if (similarity > best_similarity) {
            best = &choice;
            best_similarity = similarity;
        }
    }
    return best ? *best : std::string{};
}
