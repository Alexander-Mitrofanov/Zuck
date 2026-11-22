#include <stdexcept>
#include <string>
#include <utility>
#include <vector>
#include "crispr_array/short_spacer_filter.h"

using Rows = std::vector<std::pair<std::string, std::string>>;

static void expect(std::initializer_list<std::string> spacers, bool filtered) {
    Rows rows;
    for (const auto& spacer : spacers) rows.emplace_back("ACGT", spacer);
    if (Zuck::has_only_short_terminal_variants(rows) != filtered)
        throw std::runtime_error("Incorrect short terminal singleton classification");
}

int main() {
    const std::string short_spacer = "GACTAGAATTATTTGTTGAG";  // 20 nt
    expect({short_spacer, short_spacer + "A"}, true);
    expect({short_spacer, "AA" + short_spacer}, true);
    expect({short_spacer, short_spacer + "AA"}, true);
    expect({short_spacer, short_spacer + "AAA"}, false);  // outside bound
    expect({short_spacer, "A" + short_spacer + "A"}, false);  // internal
    expect({short_spacer, "CACTAGAATTATTTGTTGAG"}, false);  // real variant
    expect({short_spacer, short_spacer + "A", "TCAATTCGAGCTAGTACGAT"}, false);
    expect({short_spacer + "AAA", short_spacer + "AAAA"}, false);  // >=23
    expect({short_spacer, short_spacer + "AA", short_spacer + "AAAA"}, false);
    expect({short_spacer}, false);
    expect({}, false);
    expect({short_spacer, short_spacer}, false);  // exact dedup is separate
}
