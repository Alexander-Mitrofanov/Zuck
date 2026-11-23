// Reference-free, deterministic sample of strand-canonical FASTQ word types.
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <zlib.h>

#include "../libs/megahit/src/sequence/io/kseq.h"
// gzread can return a partial final block together with a gzip framing error.
// Check immediately, before a parser asks zlib for another block.
static int read_chunk(gzFile file, void *buffer, unsigned size) {
    int n = gzread(file, buffer, size);
    int status = Z_OK;
    gzerror(file, &status);
    if (status != Z_OK && status != Z_STREAM_END)
        throw std::runtime_error("Damaged or incomplete gzip stream");
    return n;
}
KSEQ_INIT(gzFile, read_chunk)

namespace {
constexpr int k = 31;
constexpr int minimum_q = 15;
constexpr uint16_t last_bin = 256;
constexpr size_t maximum_types = 200000;
constexpr uint64_t word_mask = (uint64_t{1} << (2 * k)) - 1;

// SplitMix64 finalizer is a permutation. A sampled key retains its identity.
uint64_t mix(uint64_t key) {
    key ^= key >> 30;
    key *= UINT64_C(0xbf58476d1ce4e5b9);
    key ^= key >> 27;
    key *= UINT64_C(0x94d049bb133111eb);
    return key ^ (key >> 31);
}

struct Spectrum {
    std::unordered_map<uint64_t, uint16_t> counts;
    uint64_t mask = 31;
    uint64_t reads = 0, bases = 0, windows = 0;
    double correct_windows = 0;
    std::array<int, 256> codes;
    std::array<double, 94> log_correct{};

    Spectrum() {
        codes.fill(-1);
        const std::string upper = "ACGT", lower = "acgt";
        for (int i = 0; i < 4; ++i) codes[upper[i]] = codes[lower[i]] = i;
        for (int q = minimum_q; q < 94; ++q) {
            log_correct[q] = std::log1p(-std::pow(10.0, -q / 10.0));
        }
    }

    void observe_type(uint64_t canonical_word) {
        const auto key = mix(canonical_word);
        if ((key & mask) != 0) return;
        auto &n = counts[key];
        if (n < last_bin) ++n;
        if (counts.size() > maximum_types) {
            // Every occurrence of any retained key was counted at the old rate.
            mask = (mask << 1) | 1;
            for (auto i = counts.begin(); i != counts.end();) {
                if ((i->first & mask) != 0) i = counts.erase(i);
                else ++i;
            }
        }
    }

    void add_read(const kseq_t &record) {
        const auto length = record.seq.l;
        if (length == 0 || record.qual.l != length)
            throw std::runtime_error("Expected FASTQ sequence with equally long Phred+33 quality");
        ++reads;
        bases += length;
        uint64_t forward = 0, reverse = 0;
        size_t valid = 0;
        double log_probability = 0;
        for (size_t i = 0; i < length; ++i) {
            const int code = codes[static_cast<unsigned char>(record.seq.s[i])];
            const int q = static_cast<unsigned char>(record.qual.s[i]) - 33;
            if (q < 0 || q >= 94) throw std::runtime_error("Invalid FASTQ quality byte");
            if (code < 0 || q < minimum_q) {
                valid = 0;
                forward = reverse = 0;
                log_probability = 0;
                continue;
            }
            forward = ((forward << 2) | uint64_t(code)) & word_mask;
            reverse = (reverse >> 2) | (uint64_t(3 - code) << (2 * (k - 1)));
            log_probability += log_correct[q];
            if (valid >= k) {
                const int previous_q = static_cast<unsigned char>(record.qual.s[i - k]) - 33;
                log_probability -= log_correct[previous_q];
            }
            if (++valid >= k) {
                ++windows;
                correct_windows += std::exp(log_probability);
                observe_type(std::min(forward, reverse));
            }
        }
    }

    void scan(const char *path) {
        std::unique_ptr<gzFile_s, decltype(&gzclose)> file(gzopen(path, "rb"), gzclose);
        if (!file) throw std::runtime_error("Could not open input FASTQ");
        gzbuffer(file.get(), 256 * 1024);
        std::unique_ptr<kseq_t, decltype(&kseq_destroy)> parser(kseq_init(file.get()), kseq_destroy);
        int status;
        while ((status = kseq_read(parser.get())) >= 0) add_read(*parser);
        int zstatus = Z_OK;
        gzerror(file.get(), &zstatus);
        if (status != -1 || (zstatus != Z_OK && zstatus != Z_STREAM_END))
            throw std::runtime_error("Damaged or incomplete FASTQ stream");
        if (gzclose(file.release()) != Z_OK) throw std::runtime_error("Incomplete gzip file");
    }

    void write() const {
        std::array<uint64_t, last_bin + 1> histogram{};
        for (const auto &entry : counts) ++histogram[entry.second];
        std::cout << std::setprecision(14) << "{\"k\":" << k
                  << ",\"reads\":" << reads << ",\"bases\":" << bases
                  << ",\"minimum_phred\":" << minimum_q << ",\"eligible_windows\":" << windows
                  << ",\"correct_window_exposure\":" << (bases ? correct_windows / bases : 0)
                  << ",\"sampling_denominator\":" << mask + 1 << ",\"histogram\":[";
        for (size_t i = 0; i < histogram.size(); ++i) {
            if (i) std::cout << ',';
            std::cout << histogram[i];
        }
        std::cout << "]}\n";
    }
};
}  // namespace

int main(int argc, char **argv) {
    if (argc < 2) {
        std::cerr << "Usage: zuck-kmer-spectrum READS.fastq[.gz] [MORE.fastq[.gz] ...]\n";
        return 2;
    }
    try {
        Spectrum spectrum;
        for (int i = 1; i < argc; ++i) spectrum.scan(argv[i]);
        spectrum.write();
    } catch (const std::exception &e) {
        std::cerr << "zuck-kmer-spectrum: " << e.what() << '\n';
        return 1;
    }
    return 0;
}
