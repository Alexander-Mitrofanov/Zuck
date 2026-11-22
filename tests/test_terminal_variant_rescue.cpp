#include <stdexcept>
#include <string>
#include <vector>
#include "crispr_array/terminal_variant_rescue.h"

static void require(bool value) {
    if (!value) throw std::runtime_error("terminal variant rescue invariant failed");
}
int main() {
    const std::string anchor = "ACGTTAGCGTCAAGGCTAACGTAGCTGACGAT";
    std::vector<std::string> alleles{anchor, anchor.substr(1), anchor.substr(0, anchor.size()-1)};
    std::vector<std::size_t> support{3,3,3};
    std::vector<int> component{0,1,2};
    require((Zuck::co_dominant_terminal_variants(alleles,support,component,0)==std::vector<int>{1,2}));
    support[1]=2;
    require((Zuck::co_dominant_terminal_variants(alleles,support,component,0)==std::vector<int>{2}));
    support={2,2,2};
    require(Zuck::co_dominant_terminal_variants(alleles,support,component,0).empty());
    support={3,3,3};
    alleles[1]=anchor.substr(2);  // two-base terminal change exceeds bound
    alleles[2]=anchor.substr(0,10)+anchor.substr(11);  // internal deletion
    require(Zuck::co_dominant_terminal_variants(alleles,support,component,0).empty());
    alleles={anchor,anchor.substr(1)};support={3,3};component={0,1};
    for (int i=2;i<9;++i) {
        alleles.push_back(anchor+std::string(i,'A'));support.push_back(1);component.push_back(i);
    }
    require(Zuck::co_dominant_terminal_variants(alleles,support,component,0).empty());
    component.pop_back();
    require((Zuck::co_dominant_terminal_variants(alleles,support,component,0)==std::vector<int>{1}));
}
