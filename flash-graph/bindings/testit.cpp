#include "CGraph.h"

int main() {

    auto G = fg::CGraph("test-d.adj", "test-d.index",
            "/home/disa/Research/graphyti/graphyti/src/conf/run_graph.txt");

        auto c = G.coreness();
        printf("The cores are: \n[ ");
        for (auto const& v : c) {
            printf("%lu ", v);
        }
        printf("]\n");

        auto d = G.degree();
        printf("The degree is: \n[ ");
        for (auto const& v : c) {
            printf("%lu ", v);
        }
        printf("]\n");
}
