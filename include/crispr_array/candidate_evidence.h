#ifndef ZUCK_CANDIDATE_EVIDENCE_H
#define ZUCK_CANDIDATE_EVIDENCE_H

#include "crispr_array/element_confidence.h"
#include "io/result_bundle.h"

namespace Zuck { namespace candidate_evidence {

// Preserve exact oriented observations before any repeat-family or spacer
// collapse. These are candidates, not additional accepted CRISPR calls.
class Ledger {
    using Pair = std::pair<std::string, std::string>;
    struct Evidence {
        std::size_t traversals = 0;
        std::set<std::string> cores;
    };
    std::map<Pair, Evidence> observations_;
    std::vector<Pair> pairs_;
    element_confidence::Result evidence_;
public:
    void observe(const std::string& repeat, const std::string& spacer,
                 const std::string& core) {
        auto& row = observations_[{repeat, spacer}];
        ++row.traversals;
        row.cores.insert(element_confidence::circular_identity(core));
    }

    void score(const std::vector<std::string>& reads) {
        pairs_.clear();
        for (const auto& item : observations_) pairs_.push_back(item.first);
        // Read evidence is assessed before the primary representative is
        // chosen. Catalogue/family preferences are deliberately not reported:
        // they would depend on this much larger candidate universe.
        evidence_ = element_confidence::annotate(pairs_, reads);
    }

    template<class Arrays, class Lengths>
    void write(const std::filesystem::path& directory, const Arrays& arrays, const Lengths& lengths,
               std::map<std::string, std::size_t>& counts) const {
        std::set<Pair> emitted;
        for (const auto& array : arrays)
            emitted.insert(array.entries.begin(), array.entries.end());
        const auto path = directory / "candidate_alleles.tsv";
        std::ofstream out(path);
        if (!out) throw std::runtime_error("Cannot stage candidate evidence: " + path.string());
        out << "candidate_id\trepeat_variant\tspacer\ttraversal_observations\tdistinct_circular_cores"
               "\temitted_primary_pair\tpasses_length_constraints\trepeat_score\tspacer_score"
               "\tcontext_score\tcontext_sequence_units\tcontext_start_offsets"
               "\tcontext_units_saturated\tcontext_offsets_saturated\tcontext_observability\n";
        std::size_t supported = 0, omitted_supported = 0;
        for (std::size_t i = 0; i < pairs_.size(); ++i) {
            const auto& pair = pairs_[i];
            const auto& observation = observations_.at(pair);
            const auto& row = evidence_.rows.at(i);
            const bool retained = emitted.count(pair) != 0;
            if (row.context.accepted()) { ++supported; if (!retained) ++omitted_supported; }
            out << "Candidate_" << i + 1 << '\t' << pair.first << '\t' << pair.second
                << '\t' << observation.traversals << '\t' << observation.cores.size()
                << '\t' << retained << '\t' << lengths.pair_allowed(pair.first, pair.second) << '\t';
            element_confidence::write_score(out, row.repeat_score); out << '\t';
            element_confidence::write_score(out, row.spacer_score); out << '\t';
            element_confidence::write_score(out, row.context_score);
            out << '\t' << row.context.sequence_units << '\t' << row.context.start_offsets
                << '\t' << row.context.units_saturated << '\t' << row.context.offsets_saturated
                << '\t' << row.context_observability << '\n';
        }
        result_bundle::close_checked(out, path);
        counts["candidate_alleles"] = pairs_.size();
        counts["candidate_context_supported_pairs"] = supported;
        counts["candidate_context_supported_unemitted_pairs"] = omitted_supported;
        counts["candidate_read_records"] = evidence_.counts.at("confidence_read_records");
    }
};

}} // namespace Zuck::candidate_evidence
#endif
