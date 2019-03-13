#include "CGraph.h"

int main() {

    auto G = fg::CGraph("test-u.adj", "test-u.index",
            "/home/disa/Research/graphyti/package/src/conf/run_graph.txt");
        auto c = G.coreness();

        printf("The cores are: \n[ ");
        for (auto const& v : c) {
            printf("%lu ", v);
        }
        printf("\n");
}
