#ifndef EVALUATION_H_
#define EVALUATION_H_

#include <cstddef>
#include <string>
#include <vector>

// Unit-cost Levenshtein distance. Uses O(min(a.size(), b.size())) memory.
std::size_t get_levenshtein_distance(const std::string& a, const std::string& b);

// Normalized edit similarity in [0, 1]; two empty sequences are identical.
float get_string_similarity(const std::string& a, const std::string& b);

// Count additional (including overlapping) occurrences of each nonempty
// spacer in the expected sequence. Repeated input spacers are counted again.
int get_number_of_duplicate_spacers(
    const std::vector<std::string>& spacers,
    const std::string& expected_sequence);

// Return the first most-similar choice without modifying choices.
// Return an empty string when there are no choices.
std::string get_most_similar_sequence(
    const std::string& sequence,
    const std::vector<std::string>& choices);

#endif
