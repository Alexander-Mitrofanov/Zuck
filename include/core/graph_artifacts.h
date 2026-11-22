#ifndef ZUCK_GRAPH_ARTIFACTS_H
#define ZUCK_GRAPH_ARTIFACTS_H

#include <algorithm>
#include <cctype>
#include <string>

namespace Zuck {

// Names owned by MEGAHIT graph construction. Share this list between input
// protection and cleanup; custom graph directories may also contain user data.
inline bool graph_artifact(const std::string& name) {
    const auto numbered = [&](const std::string& prefix) {
        return name.size() > prefix.size() && name.compare(0, prefix.size(), prefix) == 0 &&
            std::all_of(name.begin() + prefix.size(), name.end(),
                        [](unsigned char c) { return std::isdigit(c); });
    };
    return name == "data.lib" || name == "outfile_prefix.bin" ||
        name == "outfile_prefix.lib_info" || name == "graph.sdbg_info" ||
        name == "graph.counting" || numbered("graph.sdbg.") || numbered("graph.mercy_cand.");
}

}  // namespace Zuck
#endif
