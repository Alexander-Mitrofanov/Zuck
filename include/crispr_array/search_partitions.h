#ifndef ZUCK_SEARCH_PARTITIONS_H
#define ZUCK_SEARCH_PARTITIONS_H

#include <algorithm>
#include <cstddef>
#include <stdexcept>
#include <omp.h>

namespace Zuck {

// Partition identity determines visited history. A worker can execute any
// number of partitions, but resets its scratch state before each one. Changing
// the requested or effective OpenMP team cannot change the start sequence
// within a partition. The fixed partition count remains an inference setting.
template<class Reset, class Visit>
void visit_search_partitions(std::size_t starts, int partitions, int workers,
                             Reset reset, Visit visit) {
    if (partitions < 1 || workers < 1)
        throw std::invalid_argument("Search partitions and workers must be positive");
    if (!starts) return;
    const int count = static_cast<int>(std::min(starts, static_cast<std::size_t>(partitions)));
    const int team = std::min(workers, count);
    const std::size_t width = starts / count, remainder = starts % count;
    #pragma omp parallel for num_threads(team) schedule(dynamic, 1)
    for (int partition = 0; partition < count; ++partition) {
        const int worker = omp_get_thread_num();
        reset(worker);
        const auto p = static_cast<std::size_t>(partition);
        const std::size_t begin = p * width + std::min(p, remainder);
        const std::size_t end = begin + width + (p < remainder);
        for (std::size_t index = begin; index < end; ++index) visit(index, worker);
    }
}

} // namespace Zuck
#endif
