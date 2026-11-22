#include <filesystem>
#include <fstream>
#include <random>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>
#include "crispr_array/directed_haplotype_paths.h"

namespace dp = Zuck::directed_haplotype_paths;
struct Array { std::string consensus; std::vector<std::pair<std::string, std::string>> entries; };

static void require(bool condition, const char* reason) {
    if (!condition) throw std::runtime_error(reason);
}

int main() {
    const auto root = std::filesystem::temp_directory_path() /
        ("zuck-directed-path-test-" + std::to_string(std::random_device{}()));
    std::filesystem::create_directory(root);
    try {
        std::mt19937 random(761923);
        auto dna = [&](int n) {
            std::string value;
            while (n--) value += "ACGT"[random() % 4];
            return value;
        };
        std::string repeat = dna(30);
        if (Zuck::read_linked::reverse_complement(repeat) < repeat)
            repeat = Zuck::read_linked::reverse_complement(repeat);
        Array array{repeat, {}};
        std::vector<std::string> spacers;
        for (int i = 0; i < 5; ++i) {
            spacers.push_back(dna(32));
            array.entries.push_back({repeat, spacers.back()});
        }
        auto junction = [&](int from, int to, int flank) {
            return spacers[from].substr(spacers[from].size() - flank) + repeat +
                spacers[to].substr(0, flank);
        };
        auto edge = [&](std::vector<std::string>& reads, int from, int to,
                        int support = 2) {
            for (int i = 0; i < support; ++i)
                reads.push_back(junction(from, to, 20 + i));
        };
        auto fastq = [&](const std::string& name,
                         const std::vector<std::string>& reads, int quality = 40) {
            const auto path = root / name;
            std::ofstream out(path);
            for (std::size_t i = 0; i < reads.size(); ++i)
                out << '@' << i << '\n' << reads[i] << "\n+\n"
                    << std::string(reads[i].size(), static_cast<char>(33 + quality)) << '\n';
            return path.string();
        };
        auto selected = [&](std::size_t count) {
            std::set<dp::RowKey> rows;
            for (std::size_t i = 1; i <= count; ++i) rows.insert({0, i});
            return rows;
        };
        auto path_members = [](const dp::Result& result) {
            std::set<std::size_t> members;
            for (const auto& group : result.groups) if (group.kind == "path")
                members.insert(group.entries.begin(), group.entries.end());
            return members;
        };

        std::vector<std::string> reads;
        edge(reads, 0, 1); edge(reads, 1, 2);
        const auto full_path = dp::infer(std::vector<Array>{array},
            {fastq("full.fq", reads)}, selected(3));
        require(full_path.groups.size() == 1 &&
                path_members(full_path) == std::set<std::size_t>({1, 2, 3}),
                "complete family path was not emitted");
        require(full_path.counts.at("directed_path_memberships") == 3 &&
                full_path.counts.at("directed_path_accepted_edges") == 2,
                "complete path counts differ");

        const auto partition = dp::infer(std::vector<Array>{array},
            {fastq("partition.fq", reads)}, selected(5));
        require(partition.groups.size() == 3 &&
                partition.counts.at("directed_path_path_groups") == 1 &&
                partition.counts.at("directed_path_singleton_groups") == 2 &&
                partition.counts.at("directed_path_memberships") == 5,
                "residual rows were not partitioned exactly once");

        std::vector<std::string> reverse_reads = reads;
        for (auto& value : reverse_reads)
            value = Zuck::read_linked::reverse_complement(value);
        const auto reverse = dp::infer(std::vector<Array>{array},
            {fastq("reverse.fq", reverse_reads)}, selected(3));
        require(path_members(reverse) == std::set<std::size_t>({1, 2, 3}),
                "reverse-complement reads changed the directed path");

        std::vector<std::string> two; edge(two, 0, 1);
        const auto pair = dp::infer(std::vector<Array>{array},
            {fastq("pair.fq", two)}, selected(2));
        require(pair.groups.size() == 1 && pair.groups[0].kind == "path",
                "supported two-vertex path was rejected");

        std::vector<std::string> weak{junction(0, 1, 20)};
        const auto insufficient = dp::infer(std::vector<Array>{array},
            {fastq("weak.fq", weak)}, selected(2));
        require(insufficient.counts.at("directed_path_path_groups") == 0 &&
                insufficient.counts.at("directed_path_singleton_groups") == 2,
                "one read manufactured a path");

        std::vector<std::string> conflict; edge(conflict, 0, 1); edge(conflict, 1, 0);
        const auto reverse_conflict = dp::infer(std::vector<Array>{array},
            {fastq("conflict.fq", conflict)}, selected(2));
        require(reverse_conflict.counts.at("directed_path_path_groups") == 0 &&
                reverse_conflict.counts.at("directed_path_reverse_conflict_edges") == 2,
                "reverse-direction conflict was not rejected");

        std::vector<std::string> ranked; edge(ranked, 0, 1, 4); edge(ranked, 0, 2, 2);
        const auto unique_best = dp::infer(std::vector<Array>{array},
            {fastq("ranked.fq", ranked)}, selected(3));
        require(path_members(unique_best) == std::set<std::size_t>({1, 2}) &&
                unique_best.counts.at("directed_path_competing_edges") == 1,
                "unique reciprocal best edge was not selected");

        std::vector<std::string> tie; edge(tie, 0, 1); edge(tie, 0, 2);
        const auto tied = dp::infer(std::vector<Array>{array},
            {fastq("tie.fq", tie)}, selected(3));
        require(tied.counts.at("directed_path_path_groups") == 0 &&
                tied.counts.at("directed_path_singleton_groups") == 3,
                "tied branch was cut arbitrarily");

        std::vector<std::string> cycle;
        edge(cycle, 0, 1); edge(cycle, 1, 2); edge(cycle, 2, 0);
        const auto cyclic = dp::infer(std::vector<Array>{array},
            {fastq("cycle.fq", cycle)}, selected(3));
        require(cyclic.counts.at("directed_path_path_groups") == 0 &&
                cyclic.counts.at("directed_path_cycle_edges") == 3,
                "directed cycle was cut arbitrarily");

        auto duplicate = array;
        duplicate.entries[1].second = duplicate.entries[0].second;
        std::vector<std::string> ambiguous_reads;
        for (int flank : {20, 21})
            ambiguous_reads.push_back(spacers[0].substr(spacers[0].size() - flank) +
                                      repeat + spacers[0].substr(0, flank));
        const auto ambiguous = dp::infer(std::vector<Array>{duplicate},
            {fastq("ambiguous.fq", ambiguous_reads)}, selected(2));
        require(ambiguous.counts.at("directed_path_path_groups") == 0,
                "duplicated endpoint was assigned by sequence borrowing");

        const std::string palindrome = "ACGTACGTACGTACGTACGTACGTACGT";
        Array palindromic{palindrome, {{palindrome, spacers[0]}, {palindrome, spacers[1]}}};
        const auto pal = dp::infer(std::vector<Array>{palindromic},
            {fastq("palindrome.fq", {})}, selected(2));
        require(pal.counts.at("directed_path_palindromic_haplotype_rows") == 2 &&
                pal.counts.at("directed_path_singleton_groups") == 2,
                "palindromic haplotype was directionally invented");

        dp::write_tables(root, partition);
        require(std::filesystem::file_size(root / "directed_haplotype_paths.tsv") > 100 &&
                std::filesystem::file_size(root / "directed_haplotype_edges.tsv") > 100,
                "directed partition ledgers were not written");
    } catch (...) {
        std::filesystem::remove_all(root);
        throw;
    }
    std::filesystem::remove_all(root);
}
