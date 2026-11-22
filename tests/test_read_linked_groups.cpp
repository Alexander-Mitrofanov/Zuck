#include <filesystem>
#include <fstream>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>
#include "crispr_array/read_linked_groups.h"

namespace rl = Zuck::read_linked;
struct Array { std::string consensus; std::vector<std::pair<std::string, std::string>> entries; };
static void require(bool condition, const char* reason) {
    if (!condition) throw std::runtime_error(reason);
}
int main() {
    const auto root = std::filesystem::temp_directory_path() /
        ("zuck-read-links-test-" + std::to_string(std::random_device{}()));
    std::filesystem::create_directory(root);
    try {
        std::mt19937 random(92231);
        auto dna = [&](int n) { std::string s; while (n--) s += "ACGT"[random()%4]; return s; };
        const std::string repeat = dna(30);
        std::vector<std::string> spacers;
        Array parent{repeat, {}};
        for (int i=0; i<5; ++i) { spacers.push_back(dna(32)); parent.entries.push_back({repeat,spacers.back()}); }
        // Both genuinely observed orientations must remain referenced by a view.
        parent.entries.push_back({rl::reverse_complement(repeat),rl::reverse_complement(spacers[1])});
        const auto original = parent.entries;
        auto junction = [&](int l, int r, int flank) {
            return spacers[l].substr(32-flank)+repeat+spacers[r].substr(0,flank);
        };
        auto edge = [&](std::vector<std::string>& reads,int l,int r) {
            reads.push_back(junction(l,r,20)); reads.push_back(junction(l,r,21));
        };
        auto fastq = [&](const std::string& name,const std::vector<std::string>& reads,int q=40) {
            auto file=root/name;std::ofstream out(file);
            for (std::size_t i=0;i<reads.size();++i)
                out<<'@'<<i<<'\n'<<reads[i]<<"\n+\n"<<std::string(reads[i].size(),char(33+q))<<'\n';
            return file.string();
        };
        std::vector<std::string> reads;edge(reads,0,1);edge(reads,1,2);
        const auto path=fastq("good.fq",reads);
        for (int invalid_case = 0; invalid_case < 5; ++invalid_case) {
            rl::Config config;
            if (invalid_case == 0) config.minimum_repeat = -1;
            if (invalid_case == 1) config.maximum_repeat = config.minimum_repeat - 1;
            if (invalid_case == 2) config.minimum_quality = 94;
            if (invalid_case == 3) config.maximum_expected_errors = std::numeric_limits<double>::quiet_NaN();
            if (invalid_case == 4) config.minimum_sequence_units = 0;
            bool rejected = false;
            try { rl::infer(std::vector<Array>{parent}, {path}, config); }
            catch (const std::runtime_error&) { rejected = true; }
            require(rejected, "Invalid evidence configuration must fail before scanning reads");
        }
        auto result=rl::infer(std::vector<Array>{parent},{path});
        require(result.groups.size()==1,"supported complete path missing");
        require(result.scaffold_groups.size()==1,
                "2/2 supported partial path did not enter scaffold channel");
        require(result.groups[0].entry_indices==std::vector<std::size_t>({1,2,3,6}),"view lost literal parent orientations");
        require(parent.entries==original,"input evidence changed");
        for (const auto& e:result.groups[0].edges)
            require(e.sequence_units==2 && e.start_offsets==2,"support ledger wrong");
        auto strong_edge = [&](std::vector<std::string>& target,int l,int r) {
            target.push_back(junction(l,r,20)); target.push_back(junction(l,r,21));
            target.push_back(junction(l,r,22));
        };
        std::vector<std::string> full_path;
        strong_edge(full_path,0,1);strong_edge(full_path,1,2);
        strong_edge(full_path,2,3);strong_edge(full_path,3,4);
        const auto full=rl::infer(std::vector<Array>{parent},{fastq("full.fq",full_path)});
        require(full.groups.empty() && full.scaffold_groups.size()==1 &&
                full.scaffold_groups[0].entry_indices.size()==parent.entries.size() &&
                full.scaffold_components.size()==1 &&
                full.scaffold_components[0].status=="accepted_full_parent_path",
                "full-parent path was not retained only for scaffolding");
        Array pair{repeat,{{repeat,spacers[0]},{repeat,spacers[1]}}};
        std::vector<std::string> pair_reads;strong_edge(pair_reads,0,1);
        const auto two=rl::infer(std::vector<Array>{pair},{fastq("two.fq",pair_reads)});
        require(two.groups.empty() && two.scaffold_groups.size()==1 &&
                two.scaffold_groups[0].entry_indices==std::vector<std::size_t>({1,2}) &&
                two.scaffold_components[0].status=="accepted_two_node",
                "strong two-node scaffold was not retained");
        // Raw read reverse-complement orientation must give the same content view.
        for (auto& s:reads) s=rl::reverse_complement(s);
        require(rl::infer(std::vector<Array>{parent},{fastq("reverse.fq",reads)}).groups.size()==1,"reverse links lost");
        reads.clear();edge(reads,0,1);edge(reads,1,2);
        // Identical input duplication cannot manufacture support even when mates
        // provide two observed offsets within the same repeated pair.
        const auto duplicate1=fastq("dupe1.fq",{reads[0],reads[0],reads[2],reads[2]});
        const auto duplicate2=fastq("dupe2.fq",{reads[1],reads[1],reads[3],reads[3]});
        require(rl::infer(std::vector<Array>{parent},{duplicate1,duplicate2}).groups.empty(),"duplicates created support");
        require(rl::infer(std::vector<Array>{parent},{fastq("low.fq",reads,9)}).groups.empty(),"low quality accepted");
        require(rl::infer(std::vector<Array>{parent},{fastq("manyerrors.fq",reads,10)}).groups.empty(),"expected-error budget ignored");
        std::string q15;
        for(std::size_t i=0;i<reads.size();++i) {
            std::string quality(reads[i].size(),'I');quality[0]='0'; // Q15, bounded total expected error.
            q15+='@'+std::to_string(i)+'\n'+reads[i]+"\n+\n"+quality+'\n';
        }
        {std::ofstream out(root/"q15.fq");out<<q15;}
        require(rl::infer(std::vector<Array>{parent},{(root/"q15.fq").string()}).groups.size()==1,"safe exact Q15 junction rejected");
        std::vector<std::string> mismatch=reads;
        for(auto& s:mismatch) s[1]=s[1]=='A'?'C':'A';
        require(rl::infer(std::vector<Array>{parent},{fastq("mismatch.fq",mismatch)}).groups.empty(),"available flank mismatch accepted");
        auto ambiguous=parent;
        ambiguous.entries.push_back({repeat,dna(16)+spacers[0].substr(16)});
        require(rl::infer(std::vector<Array>{ambiguous},{path}).groups.empty(),"ambiguous endpoint accepted");
        auto branch=reads;edge(branch,1,3);
        const auto branch_result=rl::infer(
            std::vector<Array>{parent},{fastq("branch.fq",branch)});
        require(branch_result.groups.empty() && branch_result.scaffold_groups.size()==1 &&
                branch_result.scaffold_groups[0].entry_indices==
                    std::vector<std::size_t>({1,2,3,6}) &&
                branch_result.scaffold_components.size()==1 &&
                branch_result.scaffold_components[0].status==
                    "resolved_branched_forest" &&
                branch_result.scaffold_components[0].selected_group_count==1 &&
                branch_result.counts.at(
                    "read_linked_scaffold_forest_rejected_degree_edges")==1,
                "branched component did not produce a logged maximum-evidence path forest");
        auto cycle=reads;edge(cycle,2,0);
        const auto cycle_result=rl::infer(
            std::vector<Array>{parent},{fastq("cycle.fq",cycle)});
        require(cycle_result.groups.empty() && cycle_result.scaffold_groups.size()==1 &&
                cycle_result.scaffold_groups[0].entry_indices==
                    std::vector<std::size_t>({1,2,3,6}) &&
                cycle_result.scaffold_components[0].status=="resolved_cycle_forest" &&
                cycle_result.counts.at(
                    "read_linked_scaffold_forest_rejected_cycle_edges")==1,
                "cyclic component did not produce a logged maximum-evidence path forest");
        {std::ifstream in(path);std::string text((std::istreambuf_iterator<char>(in)),{});gzFile gz=gzopen((root/"good.fq.gz").c_str(),"wb");require(gz!=nullptr,"gzip output open");gzwrite(gz,text.data(),text.size());gzclose(gz);}
        require(rl::infer(std::vector<Array>{parent},{(root/"good.fq.gz").string()}).groups.size()==1,"gzip changed inference");
        {std::ofstream out(root/"reads.fa");for(std::size_t i=0;i<reads.size();++i)out<<'>'<<i<<'\n'<<reads[i]<<'\n';}
        const auto fasta=rl::infer(std::vector<Array>{parent},{(root/"reads.fa").string()});
        require(fasta.groups.empty() && fasta.counts.at("read_linked_records_without_quality")==4,"FASTA quality was invented");
        auto throws=[&](const std::vector<std::string>& files){bool failed=false;try{rl::infer(std::vector<Array>{parent},files);}catch(const std::exception&){failed=true;}require(failed,"invalid input accepted");};
        {std::ofstream out(root/"bad.fq");out<<"@x\nACGT\n+\nII\n";} throws({(root/"bad.fq").string()});
        {std::ofstream out(root/"noquality.fq");out<<"@x\nACGT\n";} throws({(root/"noquality.fq").string()});
        {std::ofstream out(root/"names.fq");out<<"@different\nACGT\n+\nIIII\n";} throws({path,(root/"names.fq").string()});
        throws({path,fastq("short.fq",{reads[0]})}); throws({(root/"absent.fq").string()});
        rl::write_tables(root.string(),result);
        require(std::filesystem::file_size(root/"read_linked_groups.tsv")>50,"membership table missing");
        rl::write_tables(root.string(),full);
        require(std::filesystem::file_size(root/"read_linked_scaffolds.tsv")>50 &&
                std::filesystem::file_size(root/"read_linked_scaffold_components.tsv")>80,
                "scaffold decision tables missing");
    } catch (...) {std::filesystem::remove_all(root);throw;}
    std::filesystem::remove_all(root);
}
