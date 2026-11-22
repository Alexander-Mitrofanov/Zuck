#ifndef ZUCK_OPTION_VALUES_H
#define ZUCK_OPTION_VALUES_H

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

// Shared by command-line and settings-file parsing. Consume the whole value:
// standard numeric conversions alone would silently accept values like "4junk".
namespace Zuck::options {

inline bool boolean(std::string value, const std::string& name) {
    std::transform(value.begin(), value.end(), value.begin(),
                   [](unsigned char c) { return std::tolower(c); });
    if (value == "true" || value == "1" || value == "yes") return true;
    if (value == "false" || value == "0" || value == "no") return false;
    throw std::runtime_error("Invalid " + name + " value: " + value +
                             " (expected true or false)");
}

inline std::uint64_t unsigned_integer(
    const std::string& value, const std::string& name,
    std::uint64_t maximum = std::numeric_limits<std::uint64_t>::max()) {
    if (value.empty() || value.find_first_not_of("0123456789") != std::string::npos)
        throw std::runtime_error("Invalid " + name + ": expected a nonnegative integer");
    try {
        const auto number = std::stoull(value);
        if (number <= maximum) return number;
    } catch (const std::out_of_range&) {
        // Report the option name rather than the standard library function.
    }
    throw std::runtime_error("Invalid " + name + ": integer out of range");
}

inline int nonnegative_integer(const std::string& value, const std::string& name) {
    return static_cast<int>(unsigned_integer(value, name, std::numeric_limits<int>::max()));
}

// Unsuffixed values are GiB; B/K/M/G use binary units. Only decimal notation
// is accepted, making CLI and settings files independent of numeric prefixes.
inline double memory_gib(const std::string& value) {
    std::string digits = value;
    char unit = 'G';
    if (!digits.empty() && std::isalpha(static_cast<unsigned char>(digits.back()))) {
        unit = std::toupper(static_cast<unsigned char>(digits.back()));
        digits.pop_back();
    }
    const auto invalid = [&]() {
        return std::runtime_error("Invalid RAM value: " + value +
                                  " (expected a positive number with optional B/K/M/G suffix)");
    };
    if (digits.empty() || digits.find_first_not_of("0123456789.") != std::string::npos ||
        std::count(digits.begin(), digits.end(), '.') > 1 || digits == ".")
        throw invalid();
    double gib;
    try {
        gib = std::stod(digits);
    } catch (const std::exception&) {
        throw invalid();
    }
    switch (unit) {
        case 'B': gib /= 1024.0 * 1024.0 * 1024.0; break;
        case 'K': gib /= 1024.0 * 1024.0; break;
        case 'M': gib /= 1024.0; break;
        case 'G': break;
        default: throw invalid();
    }
    if (!std::isfinite(gib) || gib <= 0.0) throw invalid();
    return gib;
}

// Settings lists accept whitespace, comma or semicolon separators. Quote a
// filename containing any separator; backslashes remain literal path bytes.
inline std::vector<std::string> paths(const std::string& value) {
    std::vector<std::string> result;
    std::string token;
    char quote = 0;
    for (char c : value) {
        if (quote) {
            if (c == quote) quote = 0;
            else token += c;
        } else if (c == '\'' || c == '"') {
            quote = c;
        } else if (std::isspace(static_cast<unsigned char>(c)) || c == ',' || c == ';') {
            if (!token.empty()) { result.push_back(token); token.clear(); }
        } else {
            token += c;
        }
    }
    if (quote) throw std::runtime_error("Unclosed quote in file list: " + value);
    if (!token.empty()) result.push_back(token);
    return result;
}

}  // namespace Zuck::options
#endif
