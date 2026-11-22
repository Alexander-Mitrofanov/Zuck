#ifndef SDBG_BUILD_H
#define SDBG_BUILD_H

#include <string>
#include "core/settings.h"

// Build MEGAHIT's read library and graph using Zuck's validated settings.
class SDBGBuild {
public:
    explicit SDBGBuild(const Settings& settings);

private:
    std::string WriteLibFile(const std::vector<std::string>& inputs) const;
    void BuildLib() const;
    void BuildSDBG() const;

    const Settings& settings;
};

#endif
