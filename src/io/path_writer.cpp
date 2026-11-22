#include "io/path_writer.h"
#ifdef DEBUG
#include "io/result_bundle.h"

PathWriter::PathWriter(string, SDBG& graph, const vector<uint64_t>& path,
                       string destination, string)
    : sdbg(graph), genome_id(std::move(destination)) {
    std::filesystem::create_directories(genome_id);
    CollectPathsIntoStringStreams(path);
}

string PathWriter::FetchNodeLabel(size_t node) {
    if (node >= sdbg.size()) throw std::out_of_range("Path node is outside the graph");
    vector<uint8_t> bases(sdbg.k());
    sdbg.GetLabel(node, bases.data());
    string label;
    label.reserve(bases.size());
    for (auto base : bases)
        label += base >= 1 && base <= 4 ? "ACGT"[base - 1] : 'N';
    return label;
}

void PathWriter::CollectPathsIntoStringStreams(const vector<uint64_t>& path) {
    const auto cycles_path = std::filesystem::path(genome_id) / "cycles.txt";
    const auto labels_path = std::filesystem::path(genome_id) / "labels.txt";
    ofstream cycles(cycles_path, std::ios::app);
    ofstream labels(labels_path, std::ios::app);
    if (!cycles || !labels)
        throw std::runtime_error("Cannot open diagnostic path output in " + genome_id);
    // Preserve the historical filenames: cycles contains labels, labels IDs.
    for (auto node : path) {
        cycles << FetchNodeLabel(node) << ' ';
        labels << node << ' ';
    }
    cycles << '\n';
    labels << '\n';
    result_bundle::close_checked(cycles, cycles_path);
    result_bundle::close_checked(labels, labels_path);
}
#endif
