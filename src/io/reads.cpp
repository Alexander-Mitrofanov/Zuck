#include "io/reads.h"
#include "crispr_array/read_linked_groups.h"
#include <algorithm>
#include <stdexcept>
#include <string_view>

namespace {
template <typename Consume>
void for_each_sequence(const string& path, Consume consume) {
    Zuck::read_linked::Reader input(path);
    klibpp::KSeq record;
    while (input.next(record)) consume(record.seq);
}
} // namespace

vector<string> extract_sequences_from_fastq_file(const string& fastq_file_path) {
    vector<string> sequences;
    for_each_sequence(fastq_file_path, [&](const string& sequence) {
        sequences.push_back(sequence);
    });
    return sequences;
}

void reverse_pair_ends_sequence(string& sequence) {
    std::reverse(sequence.begin(), sequence.end());

    for (char& base : sequence) {
        switch (base) {
            case 'A': base = 'T'; break;
            case 'T': base = 'A'; break;
            case 'C': base = 'G'; break;
            case 'G': base = 'C'; break;
            case 'a': base = 't'; break;
            case 't': base = 'a'; break;
            case 'c': base = 'g'; break;
            case 'g': base = 'c'; break;
        }
    }
}

namespace {
uint64_t k_mer_to_node_id(const SDBG& sdbg, std::string_view k_mer) {
    if (k_mer.size() != sdbg.k()) return SDBG::kNullID;
    vector<uint8_t> bases(k_mer.size());
    for (size_t i = 0; i < k_mer.size(); ++i) {
        switch (k_mer[i]) {
            case 'A': case 'a': bases[i] = 1; break;
            case 'C': case 'c': bases[i] = 2; break;
            case 'G': case 'g': bases[i] = 3; break;
            case 'T': case 't': bases[i] = 4; break;
            default: return SDBG::kNullID;
        }
    }
    return sdbg.IndexBinarySearch(bases.data());
}
} // namespace

vector<uint64_t> get_read_from_sequence(
    const SDBG& sdbg,
    const unordered_set<uint64_t>& nodes_of_cycles,
    const string& sequence
) {
    const uint32_t K = sdbg.k();

    if (K == 0 || sequence.size() <= 2 * K) {
        return {};
    }

    const std::string_view sequence_view(sequence);
    const auto start_k_mer = sequence_view.substr(0, K);
    const auto end_k_mer = sequence_view.substr(sequence.size() - K, K);

    const uint64_t start_node_id = k_mer_to_node_id(sdbg, start_k_mer);
    const uint64_t end_node_id = k_mer_to_node_id(sdbg, end_k_mer);

    if (nodes_of_cycles.find(start_node_id) == nodes_of_cycles.end()
    && nodes_of_cycles.find(end_node_id) == nodes_of_cycles.end()) {
        return {};
    }

    if (start_node_id == SDBG::kNullID || end_node_id == SDBG::kNullID) return {};
    vector<uint64_t> read = {start_node_id};
    read.reserve(sequence.size() - K + 1);
    for (size_t i = 1; i < sequence.size() - K; ++i) {
        const auto k_mer = sequence_view.substr(i, K);
        const uint64_t node_id = k_mer_to_node_id(sdbg, k_mer);

        if (node_id == SDBG::kNullID) return {};
        read.push_back(node_id);
    }
    read.push_back(end_node_id);

    return read;
}

vector<vector<uint64_t>> get_reads(
    const SDBG& sdbg,
    const string& fastq_file_1,
    const optional<string>& fastq_file_2,
    const vector<vector<uint64_t>>& cycles
) {

    unordered_set<uint64_t> nodes_of_cycles;
    for (const auto& cycle : cycles) {
        for (const auto& node : cycle) {
            nodes_of_cycles.insert(node);
        }
    }

    vector<vector<uint64_t>> reads;
    const auto append_read = [&](const string& sequence) {
        auto read = get_read_from_sequence(sdbg, nodes_of_cycles, sequence);
        if (!read.empty()) reads.push_back(std::move(read));
    };
    for_each_sequence(fastq_file_1, append_read);
    if (fastq_file_2) {
        for_each_sequence(*fastq_file_2, [&](string sequence) {
            reverse_pair_ends_sequence(sequence);
            append_read(sequence);
        });
    }
    return reads;
}
