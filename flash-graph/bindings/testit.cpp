#include "CGraph.h"

int main() {

#if 0
    auto G = fg::CGraph("/mnt/nfs/disa/data/graphs/test_v50_e145.adj",
        "/mnt/nfs/disa/data/graphs/test_v50_e145.index",
        "/home/disa/Research/graphyti/graphyti/src/conf/run_graph.txt");
#else
    auto G = fg::CGraph("test-u.adj", "test-u.index",
        "/home/disa/Research/graphyti/graphyti/src/conf/run_graph.txt");
#endif

#if 0
        auto c = G.coreness();
        printf("The cores are: \n[ ");
        for (auto const& v : c) {
            printf("%lu ", v);
        }
        printf("]\n");
#endif
        printf("My graph\n%s\n", G.to_str().c_str());

        auto d = G.degree();
        printf("The degree is: \n[ ");
        for (auto const& v : d) {
            printf("%u ", v);
        }
        printf("]\n");
}
