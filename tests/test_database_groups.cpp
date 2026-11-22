#include "crispr_array/database_groups.h"
#include <iostream>

namespace db = Zuck::array_database;
namespace groups = Zuck::database_groups;
struct Array { std::string consensus; std::vector<std::pair<std::string,std::string>> entries; };
void require(bool condition, const char* message) { if (!condition) throw std::runtime_error(message); }
std::string spacer(unsigned n) {
    std::string sequence = "ACGTTGCAAGTCCGATGACT";
    for (int i=0;i<6;++i) { sequence.push_back("ACGT"[n%4]); n/=4; }
    return sequence;
}
db::Record record(std::string id, std::string repeat, std::initializer_list<unsigned> indices,
                  bool bona=true, double confidence=.9) {
    db::Record r; r.id=id; r.category=bona?"Bona-fide":"Possible";
    r.bona_fide=bona; r.confidence=confidence; r.repeats.insert(db::canonical(repeat));
    for (auto i:indices) r.spacers.insert(db::canonical(spacer(i)));
    return r;
}
int main() {
    const std::string repeat = "GTTCACTGCCGTACAGGCAGCTTAGAAA";
    std::vector<Array> arrays{{repeat,{}},{repeat,{}}};
    for (unsigned i=0;i<6;++i) arrays[1].entries.emplace_back(repeat,spacer(i));
    // Every observed orientation/multiplicity belonging to a member is retained.
    arrays[1].entries.emplace_back(db::reverse_complement(repeat),db::reverse_complement(spacer(0)));
    std::vector<db::Record> records{
        record("large",repeat,{0,1,2,3}), record("nested",repeat,{0,1,2}),
        record("other",repeat,{3,4,5}), record("strong_duplicate",repeat,{0,1,2,3},true,.95)};
    auto result=groups::infer(arrays,1,records);
    require(result.views.size()==2,"Two incomparable maximal source contents must remain");
    const auto& first=result.views[0];
    require(first.spacers.size()==4 && first.rows.size()==5,"Retain all observed member rows, including reverse orientation");
    require(first.source_array=="strong_duplicate" && first.confidence==.95,"Dedup chooses one complete stronger evidence record");
    require(result.counts.at("database_group_nested_intersections_removed")==1,"Nested source subsets must not be enumerated");
    require(groups::infer(arrays,arrays.size(),records).views.empty(),"De novo parents cannot be relabelled as verified DB groups");
    auto full=records; full.push_back(record("full",repeat,{0,1,2,3,4,5}));
    require(groups::infer(arrays,1,full).views.empty(),"Full parent evidence removes nested alternatives and is itself redundant");
    require(groups::infer(arrays,1,{record("weak",repeat,{0,1,2},false,.6)},true).views.empty(),"Possible views require four distinct spacers");
    require(groups::infer(arrays,1,{record("possible",repeat,{0,1,2,3},false,.6)},true).views.size()==1,"Four Possible source matches can define a content view");
    require(groups::infer(arrays,1,{record("possible",repeat,{0,1,2,3},false,.6)},false).views.empty(),"Disabled Possible category cannot regain influence through grouping");
    require(groups::infer(arrays,1,{record("low",repeat,{0,1,2,3},true,.74)}).views.empty(),"Low confidence cannot create a source view");
    auto wrong=record("wrong",repeat,{0,1,2,3}); wrong.repeats={db::canonical("AACGTTGCAAGTCCGATGACTAACGTT")};
    require(groups::infer(arrays,1,{wrong}).views.empty(),"Spacer-only source hits cannot group an unrelated repeat");
    std::vector<Array> large{{repeat,{}}};
    for (unsigned i=0;i<35;++i) large[0].entries.emplace_back(repeat,spacer(i));
    std::vector<db::Record> many;
    for (unsigned i=0;i<17;++i) many.push_back(record(std::to_string(i),repeat,{0,2*i+1,2*i+2}));
    const auto ambiguous=groups::infer(large,0,many);
    require(ambiguous.views.empty() && ambiguous.counts.at("database_group_ambiguous_parents_abstained")==1,
            "More than sixteen incomparable source contents must abstain without favorable selection");
    std::cout << "Database grouping maximality, ambiguity, complete membership, confidence and co-occurrence controls passed\n";
}
