#include "core/tmp_utils.h"

void trim_string(string& s) {
    s.erase(0, s.find_first_not_of(" \t\n\r"));
    s.erase(s.find_last_not_of(" \t\n\r") + 1);
}

pair<string, optional<string>> get_fastq_files_from_settings(
    const Settings& settings
) {
    const auto files = settings.effective_input_reads();
    if (files.empty() || files.size() > 2)
        throw std::runtime_error("Expected one or two input read files");
    return {files[0], files.size() == 2 ? optional<string>(files[1]) : std::nullopt};
}

vector<vector<uint64_t>> cycles_map_to_cycles(
    const unordered_map<uint64_t, vector<vector<uint64_t>>>& cycles_map
) {
    vector<vector<uint64_t>> cycles;

    for (const auto& [_, inner_cycles] : cycles_map) {
        for (const auto& cycle : inner_cycles) {
            cycles.push_back(cycle);
        }
    }

    return cycles;
}

int get_cycle_count(
    const vector<vector<uint64_t>>& cycles
) {
    return static_cast<int>(cycles.size());
}

void print_sdbg_graph_to_dot_file_convert(const string& lib_file_path) {
    SDBG sdbg;
    sdbg.LoadFromFile(lib_file_path.c_str());


    std::cout << "--------------DOT-FILE------------------" << std::endl;
    std::cout << "digraph TinyGraph {\n";

    for (size_t i = 0; i < sdbg.size(); ++i) {
        std::cout << "    " << i << ";\n";
    }

    for (size_t node = 0; node < sdbg.size(); ++node) {
        int outdegree = sdbg.EdgeOutdegree(node);
        if (outdegree > 0) {
            vector<uint64_t> outgoings(outdegree);
            if (sdbg.OutgoingEdges(node, outgoings.data()) != -1) {
                for (int i = 0; i < outdegree; ++i) {
                    const uint64_t neighbor = outgoings[i];
                    if (!sdbg.IsValidEdge(neighbor)) continue;

                    std::cout << "    " << node << " -> " << neighbor << ";\n";
                }
            }
        }
    }

    std::cout << "}\n";
    std::cout << "--------------DOT-FILE-END--------------" << std::endl;
}


string fetch_node_label(SDBG& sdbg, const size_t& node) {
    if (node >= sdbg.size()) throw std::out_of_range("Path node is outside the graph");
    string label;
    vector<uint8_t> bases(sdbg.k());
    sdbg.GetLabel(node, bases.data());
    for (auto base : bases) label += base >= 1 && base <= 4 ? "ACGT"[base - 1] : 'N';
    return label;
}

tuple<string, vector<string>, string> get_systems(
    SDBG& sdbg,
    vector<vector<uint64_t>>& ordered_cycles
) {
    if (ordered_cycles.empty() || std::any_of(ordered_cycles.begin(), ordered_cycles.end(),
            [](const auto& cycle) { return cycle.empty(); })) return {};
    size_t smallest_cycle_size = ordered_cycles.front().size();
    for (const auto& cycle : ordered_cycles) {
        if (cycle.size() < smallest_cycle_size) {
            smallest_cycle_size = cycle.size();
        }
    }

    /* Repeat node extend to the right */
    size_t extension_to_right = 0;
    for (size_t i = 0; i + 1 < smallest_cycle_size; ++i) {
        unordered_set<char> next_bp;
        for (const auto& cycle : ordered_cycles) {
            const uint64_t current_node = cycle.at(i);
            const string label = fetch_node_label(sdbg, current_node);
            next_bp.insert(label.at(0));
        }

        const bool bp_branched_off = next_bp.size() > 1;
        if (bp_branched_off) {
            unordered_set<char> next_next_bp;
            for (const auto& cycle: ordered_cycles) {
                const uint64_t current_node = cycle.at(i + 1);
                const string label = fetch_node_label(sdbg, current_node);
                next_next_bp.insert(label.at(0));
            }
            
            const bool is_point_mutation = next_next_bp.size() == 1;
            if (!is_point_mutation) {
                extension_to_right = i;
                break;
            }
        }
    }

    /* Repeat node extend to the left */
    size_t extension_to_left = 0;
    for (size_t i = 0; i + 1 < smallest_cycle_size; ++i) {
        unordered_set<char> prev_bp;
        for (const auto& cycle: ordered_cycles) {
            const uint64_t current_node = cycle.at(cycle.size() - i - 1);
            const string label = fetch_node_label(sdbg, current_node);
            prev_bp.insert(label.back());
        }

        const bool bp_branched_off = prev_bp.size() > 1;
        if (bp_branched_off) {
            unordered_set<char> prev_prev_bp;
            for (const auto& cycle: ordered_cycles) {
                const uint64_t inner_current_node = cycle.at(cycle.size() - i - 2);
                const string inner_label = fetch_node_label(sdbg, inner_current_node);
                prev_prev_bp.insert(inner_label.back());
            }
            
            const bool is_point_mutation = prev_prev_bp.size() == 1;
            if (!is_point_mutation) {
                extension_to_left = i;
                break;
            }
        }
    }

    const auto shared_extension = extension_to_left + extension_to_right;
    const size_t repeat_length = shared_extension > sdbg.k() ? shared_extension - sdbg.k() : 0;

    /* Turn cycles into repeats and spacers */
    vector<string> spacers;
    vector<string> repeats;
    for (const auto& cycle : ordered_cycles) {
        string spacer;
        string repeat;

        const size_t offset_repeat = cycle.size() - extension_to_left;

        for (size_t i = 0; i < cycle.size(); ++i) {
            const uint64_t node = cycle.at((offset_repeat + i) % cycle.size());
            const string label = fetch_node_label(sdbg, node);
            
            if (i < repeat_length) {
                repeat += label.back();
            } else {
                spacer += label.back();
            }
        }

        spacers.push_back(spacer);
        repeats.push_back(repeat);
    }

    /* Get consensus/most-common repeat */
    unordered_map<string, int> repeat_count;
    for (const string& repeat : repeats) {
        repeat_count[repeat]++;
    }

    string consensus_repeat = "";
    int max_count = 0;
    for (const auto& [repeat, count] : repeat_count) {
        if (count > max_count || (count == max_count && repeat < consensus_repeat)) {
            max_count = count;
            consensus_repeat = repeat;
        }
    }

    /* Reconstruct sequences */
    string full_sequence = "";
    for (size_t i = 0; i < spacers.size(); ++i) {
        const string repeat = repeats.at(i);
        const string spacer = spacers.at(i);

        if (repeat == consensus_repeat) {
            full_sequence += repeat;
            full_sequence += spacer;
        }
    }

    // This experimental reconstruction assumes a terminal consensus repeat.
    full_sequence += consensus_repeat;

    return std::make_tuple(consensus_repeat, spacers, full_sequence);
}