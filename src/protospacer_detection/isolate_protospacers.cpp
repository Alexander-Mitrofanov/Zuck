#include "protospacer_detection/isolate_protospacers.h"
#include <vector>
#include <sstream>
#include <algorithm>
#include <stdexcept>
#include <charconv>
#include <optional>

namespace {
uint64_t parse_node_id(const std::string& text, const std::string& filename) {
    uint64_t value = 0;
    const auto parsed = std::from_chars(text.data(), text.data() + text.size(), value);
    if (text.empty() || parsed.ec != std::errc{} || parsed.ptr != text.data() + text.size())
        throw std::runtime_error("Invalid node/group identifier in " + filename + ": " + text);
    return value;
}
} // namespace

IsolateProtospacers::IsolateProtospacers(const SDBG& g, const std::string& cyclesFilePath) : graph(g) {
    this->cycleNodes= ReadCycles(cyclesFilePath);
    protospacerNodes = IsolateProtospacersMethod();
}

IsolateProtospacers::IsolateProtospacers(const SDBG& g, const std::map<uint64_t, std::vector<std::vector<uint64_t>>>& repeatToSpacerNodes) : graph(g) {
    for (const auto& pair : repeatToSpacerNodes) {
        uint64_t groupId = pair.first;
        const auto& vecVec = pair.second;
        for (const auto& vec : vecVec) {
            if (!vec.empty()) {
                for (auto node : vec)
                    if (node >= graph.size()) throw std::runtime_error("Cycle node is outside the graph");
                uint64_t cycleId = vec[0];
                std::set<uint64_t> nodes(vec.begin(), vec.end());
                cycleNodes[cycleId].insert(nodes.begin(), nodes.end());
                cycleToGroup[cycleId] = groupId;
            }
        }
    }
    protospacerNodes = IsolateProtospacersMethod();
}

std::map<uint64_t, std::set<uint64_t>> IsolateProtospacers::ReadCycles(const std::string& filePath) {
    std::map<uint64_t, std::set<uint64_t>> cyclesMap;

    std::ifstream file(filePath);
    if (!file.is_open()) {
        throw std::runtime_error("Cannot open cycle nodes: " + filePath);
    }

    std::string line;
    while (std::getline(file, line)) {
        std::istringstream record(line);
        std::string token;
        if (!(record >> token)) continue;
        const auto first_node = parse_node_id(token, filePath);
        if (first_node >= graph.size()) throw std::runtime_error("Cycle node is outside the graph: " + filePath);
        std::set<uint64_t> nodes;
        while (record >> token) {
            const auto node = parse_node_id(token, filePath);
            if (node >= graph.size()) throw std::runtime_error("Cycle node is outside the graph: " + filePath);
            if (node != first_node) nodes.insert(node);
        }
        if (!nodes.empty()) {
            cyclesMap[first_node].insert(nodes.begin(), nodes.end());
            cycleToGroup[first_node] = first_node;
        }
    }
    if (!file.eof()) throw std::runtime_error("Cannot read cycle nodes: " + filePath);

    file.close();
    return cyclesMap;
}

IN_OUT_PAIR_MAP_SET IsolateProtospacers::IsolateProtospacersMethod() {
    std::map<uint64_t, std::set<uint64_t>> incomingOutersMap, outgoingOutersMap;

    for (const auto& pair : cycleNodes) {
        uint64_t k = pair.first;
        const std::set<uint64_t>& cycleNodesSet = pair.second;

        // Find outer incomings: nodes outside that point into the cycle
        std::set<uint64_t> incomingOuters;
        for (uint64_t cycleNode : cycleNodesSet) {
            int indegree = graph.EdgeIndegree(cycleNode);
            if (indegree > 0) {
                std::vector<uint64_t> incomings(indegree);
                graph.IncomingEdges(cycleNode, incomings.data());
                for (uint64_t neighbor : incomings) {
                    if (cycleNodesSet.find(neighbor) == cycleNodesSet.end() && cycleNodes.find(neighbor) == cycleNodes.end()) {
                        incomingOuters.insert(neighbor);
                    }
                }
            }
        }

        // Find outer outgoings: nodes outside that the cycle points to
        std::set<uint64_t> outgoingOuters;
        for (uint64_t cycleNode : cycleNodesSet) {
            int outdegree = graph.EdgeOutdegree(cycleNode);
            if (outdegree > 0) {
                std::vector<uint64_t> outgoings(outdegree);
                graph.OutgoingEdges(cycleNode, outgoings.data());
                for (uint64_t neighbor : outgoings) {
                    if (cycleNodesSet.find(neighbor) == cycleNodesSet.end() && cycleNodes.find(neighbor) == cycleNodes.end()) {
                        outgoingOuters.insert(neighbor);
                    }
                }
            }
        }

        // Only add if both exist
        if (!incomingOuters.empty() && !outgoingOuters.empty()) {
            incomingOutersMap[k] = incomingOuters;
            outgoingOutersMap[k] = outgoingOuters;
        }
    }

    return {incomingOutersMap, outgoingOutersMap};
}

std::pair<std::map<uint64_t, std::set<uint64_t>>, std::map<uint64_t, std::set<uint64_t>>> IsolateProtospacers::FilterSharedGroups(const std::map<uint64_t, std::set<uint64_t>>& inGroup, const std::map<uint64_t, std::set<uint64_t>>& outGroup) {
    std::map<uint64_t, std::set<uint64_t>> possibleInProtospacerNodes, possibleOutProtospacerNodes;
    for (const auto& pair : inGroup) {
        uint64_t k = pair.first;
        const std::set<uint64_t>& values = pair.second;
        if (outGroup.find(k) != outGroup.end()) {
            possibleInProtospacerNodes[k] = values;
        }
    }
    for (const auto& pair : outGroup) {
        uint64_t k = pair.first;
        const std::set<uint64_t>& values = pair.second;
        if (inGroup.find(k) != inGroup.end()) {
            possibleOutProtospacerNodes[k] = values;
        }
    }
    return {possibleInProtospacerNodes, possibleOutProtospacerNodes};
}

void IsolateProtospacers::DepthLimitedSearch(uint64_t currentNode, int depth, std::vector<uint64_t>& path, std::set<uint64_t>& visited, const std::set<uint64_t>& outNodes, const std::set<uint64_t>& cycleNodeSet, int maxDepth, int minDepth, std::vector<std::vector<uint64_t>>& successfulPaths) {
    if (depth > maxDepth) return;
    visited.insert(currentNode);
    path.push_back(currentNode);

    if (outNodes.find(currentNode) != outNodes.end() && depth >= minDepth) {
        successfulPaths.push_back(path);
    } else {
        int outdegree = graph.EdgeOutdegree(currentNode);
        if (outdegree > 0) {
            std::vector<uint64_t> outgoings(outdegree);
            graph.OutgoingEdges(currentNode, outgoings.data());
            for (uint64_t neighbor : outgoings) {
                if (visited.find(neighbor) == visited.end() && (cycleNodeSet.find(neighbor) != cycleNodeSet.end() || outNodes.find(neighbor) != outNodes.end()) && cycleNodes.find(neighbor) == cycleNodes.end()) {
                    DepthLimitedSearch(neighbor, depth + 1, path, visited, outNodes, cycleNodeSet, maxDepth, minDepth, successfulPaths);
                }
            }
        }
    }

    path.pop_back();
    visited.erase(currentNode);
}

// do a dls from each incoming protospacer node
std::map<uint64_t, std::map<uint64_t, std::vector<std::vector<uint64_t>>>> IsolateProtospacers::DepthLimitedPathsFromInToOut(
    const std::map<uint64_t, std::set<uint64_t>>& inGroup,
    const std::map<uint64_t, std::set<uint64_t>>& outGroup,
    int maxDepth,
    int minDepth
) {
    
    if (minDepth < 0 || maxDepth < minDepth)
        throw std::invalid_argument("Invalid protospacer search depth bounds");
    auto [possibleInProtospacerNodes, possibleOutProtospacerNodes] = FilterSharedGroups(inGroup, outGroup);
    
    
    std::map<uint64_t, std::map<uint64_t, std::vector<std::vector<uint64_t>>>> groupedPaths;
    for (const auto& inPair : possibleInProtospacerNodes) {
        uint64_t cycleStartNode = inPair.first;
        const std::set<uint64_t>& inNodes = inPair.second;

        auto outIt = possibleOutProtospacerNodes.find(cycleStartNode);
        if (outIt == possibleOutProtospacerNodes.end()) {
            continue; // No corresponding outgoing protospacer nodes
        }
        const std::set<uint64_t>& outNodes = outIt->second;

        // Get the full cycle node set
        auto cycleIt = this->cycleNodes.find(cycleStartNode);
        if (cycleIt == this->cycleNodes.end()) {
            continue; // No cycle nodes found
        }
        const std::set<uint64_t>& cycleNodeSet = cycleIt->second;

        std::vector<std::vector<uint64_t>> cyclePaths;
        for (uint64_t startNode : inNodes) {
            std::vector<uint64_t> path;
            std::set<uint64_t> visited;
            int cycleMaxDepth = static_cast<int>(std::min<std::size_t>(
                static_cast<std::size_t>(maxDepth), cycleNodeSet.size() + 1));
            DepthLimitedSearch(startNode, 0, path, visited, outNodes, cycleNodeSet, cycleMaxDepth, minDepth, cyclePaths);
        }

        // Filter out subpaths for this cycle
        // Sort by length descending
        std::sort(cyclePaths.begin(), cyclePaths.end(), [](const std::vector<uint64_t>& a, const std::vector<uint64_t>& b) {
            return a.size() != b.size() ? a.size() > b.size() : a < b;
        });
        std::vector<std::vector<uint64_t>> filtered;
        for (const auto& p : cyclePaths) {
            bool isSubpath = false;
            for (const auto& longer : filtered) {
                if (longer.size() > p.size() && std::search(longer.begin(), longer.end(), p.begin(), p.end()) != longer.end()) {
                    isSubpath = true;
                    break;
                }
            }
            if (!isSubpath) {
                filtered.push_back(p);
            }
        }

        // Select node-disjoint paths (greedy: longest first)
        std::vector<std::vector<uint64_t>> disjointPaths;
        std::set<uint64_t> usedNodes;
        for (const auto& p : filtered) {
            bool canAdd = true;
            for (uint64_t node : p) {
                if (usedNodes.count(node)) {
                    canAdd = false;
                    break;
                }
            }
            if (canAdd) {
                disjointPaths.push_back(p);
                for (uint64_t node : p) {
                    usedNodes.insert(node);
                }
            }
        }
        filtered = disjointPaths;

        // Trim each path to [1:-1] for protospacers
        std::vector<std::vector<uint64_t>> trimmedPaths;
        for (const auto& p : filtered) {
            if (p.size() > 2) {
                std::vector<uint64_t> trimmed(p.begin() + 1, p.end() - 1);
                trimmedPaths.push_back(trimmed);
            }
        }

        uint64_t groupId = cycleToGroup[cycleStartNode];
        groupedPaths[groupId][cycleStartNode] = trimmedPaths;
    }

    // Global subpath filtering across all paths
    std::vector<std::vector<uint64_t>> allPaths;
    std::map<std::vector<uint64_t>, std::set<uint64_t>> pathToCycles;
    for (const auto& group : groupedPaths) {
        for (const auto& cycle : group.second) {
            uint64_t cycleId = cycle.first;
            for (const auto& path : cycle.second) {
                allPaths.push_back(path);
                pathToCycles[path].insert(cycleId);
            }
        }
    }

    // Sort by length descending
    std::sort(allPaths.begin(), allPaths.end(), [](const std::vector<uint64_t>& a, const std::vector<uint64_t>& b) {
        return a.size() != b.size() ? a.size() > b.size() : a < b;
    });

    std::vector<std::vector<uint64_t>> globalFiltered;
    for (const auto& p : allPaths) {
        bool isSubpath = false;
        for (const auto& longer : globalFiltered) {
            if (longer.size() > p.size() && std::search(longer.begin(), longer.end(), p.begin(), p.end()) != longer.end()) {
                isSubpath = true;
                break;
            }
        }
        if (!isSubpath) {
            globalFiltered.push_back(p);
        }
    }

    // Unique the global filtered
    std::set<std::vector<uint64_t>> uniqueGlobal(globalFiltered.begin(), globalFiltered.end());
    globalFiltered.assign(uniqueGlobal.begin(), uniqueGlobal.end());

    // Assign each unique path to the cycle with the smallest ID that had it
    std::map<uint64_t, std::map<uint64_t, std::vector<std::vector<uint64_t>>>> newGroupedPaths;
    for (const auto& path : globalFiltered) {
        auto it = pathToCycles.find(path);
        if (it != pathToCycles.end()) {
            uint64_t minCycleId = *std::min_element(it->second.begin(), it->second.end());
            uint64_t groupId = cycleToGroup[minCycleId];
            newGroupedPaths[groupId][minCycleId].push_back(path);
        }
    }

    return newGroupedPaths;
}

void IsolateProtospacers::WritePathsToFile(const std::map<uint64_t, std::map<uint64_t, std::vector<std::vector<uint64_t>>>>& paths, const std::string& filename) {
    std::ofstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Cannot open protospacer paths: " + filename);
    }

    int in_counter = 0;
    for (const auto& group_paths : paths) {
        uint64_t group_id = group_paths.first;
        const auto& cycle_map = group_paths.second;
        file << "Group " << group_id << ":\n";
        for (const auto& cycle_paths : cycle_map) {
            uint64_t cycle_id = cycle_paths.first;
            const auto& path_list = cycle_paths.second;
            if (!path_list.empty()) {
                file << "  Cycle " << cycle_id << ":\n";
                for (const auto& path : path_list) {
                    in_counter += 1;
                    file << in_counter << "    [";
                    for (size_t i = 0; i < path.size(); ++i) {
                        file << path[i];
                        if (i < path.size() - 1) file << " ";
                    }
                    file << "]\n";
                }
            }
        }
    }

    file.close();
    if (!file) throw std::runtime_error("Cannot write protospacer paths: " + filename);
}


std::map<uint64_t, std::map<uint64_t, std::vector<std::vector<uint64_t>>>> IsolateProtospacers::ReadPathsFromFile(const std::string& filename) {
    std::map<uint64_t, std::map<uint64_t, std::vector<std::vector<uint64_t>>>> paths;
    std::ifstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Cannot open protospacer paths: " + filename);
    }

    std::optional<uint64_t> current_group, current_cycle;
    std::string line;
    while (std::getline(file, line)) {
        line.erase(0, line.find_first_not_of(" \t\r"));
        line.erase(line.find_last_not_of(" \t\r") + 1);
        if (line.empty()) continue;
        const auto header = [&](const std::string& prefix) {
            if (line.back() != ':')
                throw std::runtime_error("Missing header colon in " + filename);
            return parse_node_id(line.substr(prefix.size(), line.size() - prefix.size() - 1), filename);
        };
        if (line.rfind("Group ", 0) == 0) {
            current_group = header("Group ");
            current_cycle.reset();
        } else if (line.rfind("Cycle ", 0) == 0) {
            if (!current_group) throw std::runtime_error("Cycle without group in " + filename);
            current_cycle = header("Cycle ");
            if (*current_cycle >= graph.size()) throw std::runtime_error("Cycle node is outside the graph: " + filename);
        } else {
            const auto begin = line.find('[');
            const auto end = line.find(']');
            if (!current_group || !current_cycle || begin == std::string::npos ||
                end != line.size() - 1 || end <= begin)
                throw std::runtime_error("Malformed path record in " + filename);
            std::istringstream nodes(line.substr(begin + 1, end - begin - 1));
            std::string token;
            std::vector<uint64_t> path;
            while (nodes >> token) {
                const auto node = parse_node_id(token, filename);
                if (node >= graph.size()) throw std::runtime_error("Path node is outside the graph: " + filename);
                path.push_back(node);
            }
            if (path.empty()) throw std::runtime_error("Empty path in " + filename);
            paths[*current_group][*current_cycle].push_back(std::move(path));
        }
    }
    if (!file.eof()) throw std::runtime_error("Cannot read protospacer paths: " + filename);
    return paths;
}
