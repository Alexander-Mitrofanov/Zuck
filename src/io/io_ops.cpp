#include "io/io_ops.h"
#ifdef DEBUG
#include <fstream>
#include <stdexcept>
#include "io/result_bundle.h"
#include "core/option_values.h"

std::unordered_map<uint64_t, std::vector<std::vector<uint64_t>>> io_ops::read_cycles(const std::string& file_path) {
    std::unordered_map<uint64_t, std::vector<std::vector<uint64_t>>> cycles;
    std::ifstream file(file_path);
    if (!file.is_open()) {
        throw std::runtime_error("Cannot open cycle JSON: " + file_path);
    }
    const auto document = json::parse(file);
    if (!document.is_object())
        throw std::runtime_error("Cycle JSON must be an object of path arrays: " + file_path);
    for (const auto& [key, paths] : document.items()) {
        const auto cycle_id = Zuck::options::unsigned_integer(key, "cycle ID");
        if (!paths.is_array())
            throw std::runtime_error("Cycle JSON group must contain an array of paths: " + key);
        std::vector<std::vector<uint64_t>> group;
        for (const auto& path : paths) {
            if (!path.is_array())
                throw std::runtime_error("Cycle JSON path must be an array: " + key);
            std::vector<uint64_t> nodes;
            for (const auto& node : path) {
                if (!node.is_number_integer() ||
                    (!node.is_number_unsigned() && node.get<int64_t>() < 0))
                    throw std::runtime_error("Cycle JSON nodes must be nonnegative integers: " + key);
                nodes.push_back(node.get<uint64_t>());
            }
            group.push_back(std::move(nodes));
        }
        if (!cycles.emplace(cycle_id, std::move(group)).second)
            throw std::runtime_error("Cycle JSON contains duplicate numeric IDs: " + key);
    }
    return cycles;
}

void io_ops::write_cycles(const std::string& file_path, const std::unordered_map<uint64_t, std::vector<std::vector<uint64_t>>>& cycles) {
    json j = json::object();
    for (const auto& [cycle_id, cycle_vec] : cycles) {
        j[std::to_string(cycle_id)] = cycle_vec;
    }
    std::ofstream file(file_path);
    if (!file.is_open()) {
        throw std::runtime_error("Cannot create cycle JSON: " + file_path);
    }
    file << j.dump(4); // Dump with indentation of 4 spaces
    result_bundle::close_checked(file, file_path);
}


void io_ops::write_nodes_gfa(const std::string& file_path, const SDBG& sdbg) {
    std::ofstream file(file_path);
    if (!file.is_open()) {
        throw std::runtime_error("Cannot create graph GFA: " + file_path);
    }
    // Write segments (nodes)
    for (uint64_t edge_id = 0; edge_id < sdbg.size(); ++edge_id) {
        if (!sdbg.IsValidEdge(edge_id)) continue;
        // Get label as sequence (optional, here just node id as string)
        file << "S\t" << edge_id << "\t*\n";
    }
    // Write links (edges)
    for (uint64_t edge_id = 0; edge_id < sdbg.size(); ++edge_id) {
        if (!sdbg.IsValidEdge(edge_id)) continue;
        int outdeg = sdbg.EdgeOutdegree(edge_id);
        if (outdeg > 0) {
            std::vector<uint64_t> outgoings(outdeg);
            sdbg.OutgoingEdges(edge_id, outgoings.data());
            for (int i = 0; i < outdeg; ++i) {
                // Overlap is unknown, use 0M (no overlap info)
                file << "L\t" << edge_id << "\t+\t" << outgoings[i] << "\t+\t0M\n";
            }
        }
    }
    result_bundle::close_checked(file, file_path);
}

#endif
