#include "crispr_array/boundary_candidates.h"
#include <chrono>
#include <iostream>

namespace bc = Zuck::boundary_candidates;
namespace db = Zuck::array_database;
namespace ec = Zuck::element_confidence;
void require(bool ok, const char* why) { if (!ok) throw std::runtime_error(why); }
int main() {
    const bc::Pair original{"GTTCACTGCCGTACAGGCAGCTTAGAAA", "ACGTGGCATGACTGATCGTACGATCGTACGAA"};
    std::set<bc::Pair> pairs;
    for (int a=-3; a<=3; ++a) for (int b=-3; b<=3; ++b) {
        const auto p = bc::shifted(original,a,b);
        require(p.first.size()==original.first.size()+b-a, "Repeat coordinates must move independently");
        require(p.second.size()==original.second.size()+a-b, "Spacers must include both lengthenings and shortenings");
        require(ec::circular_identity(p.first+p.second)==ec::circular_identity(original.first+original.second),
                "Boundary shifts must conserve the complete observed circle");
        const bc::Pair rc{db::reverse_complement(original.first), db::reverse_complement(original.second)};
        const auto rp = bc::shifted(rc,-b,-a);
        require(rp==bc::Pair{db::reverse_complement(p.first),db::reverse_complement(p.second)},
                "Reverse orientation must map both coordinates correctly");
        pairs.insert(p);
    }
    require(pairs.size()==49,"All 49 combinations are needed, including diagonal shifts");
    require(bc::shifted(original,0,0)==original,"Zero shift must reproduce the original pair");
    const auto root=std::filesystem::temp_directory_path()/
        ("zuck-boundary-"+std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()));
    std::filesystem::create_directory(root);
    const auto reads=root/"reads.fastq";
    const auto span=original.first+original.second+original.first;
    {
        std::ofstream out(reads);
        for (const auto& p : {"AAA"+span+"CC", "T"+span+"GGG"})
            out << '@' << p.size() << '\n' << p << "\n+\n" << std::string(p.size(),'I') << '\n';
    }
    bc::Search search(3);
    search.prepare({original}, {reads.string()}, {});
    const std::vector<bc::Pair> vector_pairs(pairs.begin(),pairs.end());
    const auto direct=ec::annotate(vector_pairs,{reads.string()});
    search.each([&](std::size_t, int, int, const bc::Pair& p, std::size_t context) {
        const auto index=std::find(vector_pairs.begin(),vector_pairs.end(),p)-vector_pairs.begin();
        const auto s=search.support(context), d=direct.rows[index].context;
        require(s.sequence_units==d.sequence_units && s.start_offsets==d.start_offsets && s.index()==d.index(),
                "Shared context optimization must equal literal per-variant scoring");
    });
    bc::Search missing(3); missing.prepare({original},{},{});
    missing.each([&](std::size_t,int,int,const bc::Pair&,std::size_t context) {
        require(!missing.support(context).accepted() && !missing.support(context).observable,
                "Missing reads cannot supply support");
    });
    bc::Search zero(0); zero.prepare({original},{},{});
    require(zero.counts.at("boundary_combinations")==1,"Radius zero must be a one-cut ablation");
    auto record = [&](const std::set<bc::Pair>& seeds) {
        db::Record r; r.id="source"; r.bona_fide=true; r.category="Bona-fide"; r.confidence=0.9;
        for (const auto& p : seeds) {r.repeats.insert(db::canonical(p.first));r.spacers.insert(db::canonical(p.second));}
        return r;
    };
    // Three database spacers from different cuts of ONE observed circle must
    // not masquerade as three independent array anchors.
    const auto alias_record=record({original,bc::shifted(original,0,1),bc::shifted(original,0,-1)});
    bc::Search aliases(3);aliases.prepare({original},{reads.string()},{alias_record});
    require(aliases.select("catalogue-all").empty(),"Alternative cuts cannot corroborate their own circular core");
    const bc::Pair second{original.first,"AGGTCTCGATAGCATGCCGATATCCGTAGTCA"};
    const bc::Pair third{original.first,"TGCATCCGGTAGCTACGATGCGATCCAGTTAG"};
    const std::set<bc::Pair> three{original,second,third};
    {
        std::ofstream out(reads);std::size_t n=0;
        for (const auto& p:three) {
            const auto context=p.first+p.second+p.first;
            for(const auto& s:{"AAA"+context+"CC","T"+context+"GGG"})
                out << '@' << ++n << '\n' << s << "\n+\n" << std::string(s.size(),'I') << '\n';
        }
    }
    bc::Search two(3);two.prepare({original,second},{reads.string()},{record(three)});
    require(!two.select("catalogue-all").empty() && two.select("catalogue").empty(),
            "Default needs two other known anchors; one is only a diagnostic policy");
    bc::Search corroborated(3);corroborated.prepare(three,{reads.string()},{record(three)});
    require(corroborated.select("catalogue")==three,"Three real cores should support the known cuts");
    auto duplicate=record(three);duplicate.id="duplicate";
    bc::Search duplicated(3);duplicated.prepare(three,{reads.string()},{record(three),duplicate});
    require(duplicated.select("catalogue")==three &&
            duplicated.assessment(original).same_array_score==corroborated.assessment(original).same_array_score,
            "Duplicated sources cannot inflate evidence");
    std::vector<db::Record> isolated;
    for(const auto& p:three) isolated.push_back(record({p}));
    bc::Search separate(3);separate.prepare(three,{reads.string()},isolated);
    require(separate.select("catalogue-all").empty(),"Spacers from unrelated source arrays cannot be pooled");
    struct Array {std::string consensus;std::vector<bc::Pair> entries;};
    struct Limits {
        bool pair_allowed(const std::string&,const std::string& s) const {return s.size()>=33;}
        bool repeat_allowed(std::size_t) const {return true;}
    };
    std::vector<Array> arrays;
    std::map<std::string,std::size_t> counts;
    corroborated.append(arrays,Limits{},"catalogue",root/"selected.tsv",counts);
    require(arrays.empty(),"Length filtering leaving fewer than two spacers must suppress the group");
    search.write(root/"boundary_candidates.tsv");
    std::filesystem::remove_all(root);
    std::cout << "Independent boundary shifts, strand covariance, exact context sharing and missingness passed\n";
}
