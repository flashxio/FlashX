/*
 * Copyright 2014 Open Connectome Project (http://openconnecto.me)
 * Written by Disa Mhembere (disa@jhu.edu)
 *
 * This file is part of FlashGraph.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#include "libgraph-algs/sem_kmeans.h"
#include "libgraph-algs/sem_kmeans_util.h"
#include "libgraph-algs/dist_matrix.h"

static std::vector<prune_cluster::ptr> g_clusters;
constexpr unsigned NCOL = 5;

// Any item with an iterator can be tested for equivalence
template <typename T>
bool all_equal(const T& arg0, const T& arg1) {
    return std::equal(arg0.begin(), arg0.end(), arg1.begin());
}

std::vector<std::vector<double>> test_init_g_clusters(const size_t k=4) {
    BOOST_LOG_TRIVIAL(info) << "Running init g_clusters";
    BOOST_VERIFY(k == 4);

    const std::vector<double> v1 {1, 2, 3, 4, 5};
    const std::vector<double> v2 {6, 7, 8, 9, 10};
    const std::vector<double> v3 {6E-12, -23423.7, .82342342432, 93., 10};
    const std::vector<double> v4 {-.2342, -23.342, -.000003232, -3.234232, 1};

    std::vector<std::vector<double>> means = {v1, v2, v3, v4};

    for (size_t cl = 0; cl < k; cl++) {
        g_clusters.push_back(prune_cluster::create(means[cl])); // ctor & init

        printf("c:%lu =>\n", cl);
        print_vector(g_clusters[cl]->get_mean());
        print_vector(means[cl]);

        BOOST_VERIFY(all_equal(means[cl], g_clusters[cl]->get_mean()));
    }
    printf("Exiting test_init_g_clusters!\n");
    return means;
}

void test_cluster() {
    size_t k = 4;
    std::vector<std::vector<double>> means = test_init_g_clusters(k);
    /* test prev_mean */

    // NOTE: Rotate means clockwise: 0 => 1; 1 => 2; 2 => 3; 3 => 0
    // dist(3,0), (1, 0), (1,2), (2,3)
    const std::vector<double> pdists =
        {ceil(sqrt(721.074)), ceil(sqrt(125)),
            ceil(sqrt(549004845.993)), ceil(sqrt(547586097.288))};

    std::vector<double> tmp;
    std::vector<double> prev = g_clusters.back()->get_mean();
    for (size_t cl = 0; cl < k; cl++) {
        tmp = g_clusters[cl]->get_mean();
        g_clusters[cl]->set_prev_mean();
        g_clusters[cl]->set_mean(prev);

        // Compute dist to prev
        g_clusters[cl]->set_prev_dist(eucl_dist(&((g_clusters[cl]->get_mean())[0]),
                    &((g_clusters[cl]->get_prev_mean())[0]), NCOL) );

        /* test prev_dist */
        BOOST_VERIFY(ceil(g_clusters[cl]->get_prev_dist()) == pdists[cl]);
        prev = tmp;
    }

    /* Test operator [] */
    for (unsigned i = 0; i < g_clusters[0]->size(); i++) {
        BOOST_VERIFY((*(g_clusters[1]))[i]  == means[0][i]);
    }

    /* Test add member */
    class data_seq_it {
        private:
            std::vector<double> data;
            unsigned pos;
        public:
            data_seq_it(std::vector<double>& data) {
                pos = 0;
                this->data = data;
            }
            bool has_next() { return pos < data.size(); }
            double next() { return data[pos++]; }
    };

    // Add 0 to 1
    data_seq_it dsi0(means[0]);
    data_seq_it dsi1(means[1]);

    g_clusters[2]->add_member(dsi0);
    g_clusters[1]->add_member(dsi1);

    BOOST_VERIFY(all_equal(g_clusters[2]->get_mean(),
                g_clusters[1]->get_mean()));
    
    std::cout << "\nMembers in 2 = " <<
        g_clusters[2]->get_num_members() << std::endl;
    std::cout << "Members in 1 = " <<
        g_clusters[1]->get_num_members() << std::endl;

    BOOST_VERIFY(g_clusters[2]->get_num_members() == 1);
    BOOST_VERIFY(g_clusters[1]->get_num_members() == 1);

    // Test remove member
    dsi0 = data_seq_it(means[0]);
    dsi1 = data_seq_it(means[1]);
    g_clusters[1]->remove_member(dsi1);
    g_clusters[1]->remove_member(dsi0);
    std::vector<double> zeros {0,0,0,0,0};

    BOOST_VERIFY(all_equal(g_clusters[1]->get_mean(), zeros));
    BOOST_VERIFY(g_clusters[1]->get_num_members() == -1);
    printf("Exiting test_cluster ==> ");
}

void test_eucl() {
    // Positive
    std::vector<double> v1 {1, 2, 3, 4, 5};
    std::vector<double> v2 {6, 7, 8, 9, 10};
    BOOST_VERIFY(eucl_dist(&v1[0], &v2[0], NCOL) == sqrt(125.0));
    BOOST_VERIFY(eucl_dist(&v2[0], &v1[0], NCOL) == sqrt(125.0));

    // Neg-pos, Pos-neg
    std::vector<double> v3 {6E-12, -23423.7, .82342342432, 93., 10};
    BOOST_VERIFY(ceil(eucl_dist(&v1[0], &v3[0], NCOL)) ==
            ceil(sqrt(548771372.227)));
    BOOST_VERIFY(ceil(eucl_dist(&v3[0], &v1[0], NCOL))
            == ceil(sqrt(548771372.227)));

    // No-op
    std::vector<double> v4 {0, 0, 0, 0, 0};
    BOOST_VERIFY(eucl_dist(&v1[0], &v4[0], NCOL) ==
            eucl_dist(&v4[0], &v1[0], NCOL));
    BOOST_VERIFY(eucl_dist(&v4[0], &v1[0], NCOL) == sqrt(55));

    // Neg-neg
    std::vector<double> v5 {-.2342, -23.342, -.000003232, -3.234232, 1};
    BOOST_VERIFY(ceil(eucl_dist(&v5[0], &v3[0], NCOL))
            == ceil(sqrt(547586097.2884537)));
    BOOST_VERIFY(ceil(eucl_dist(&v3[0], &v5[0], NCOL))
            == ceil(sqrt(547586097.2884537)));

    double arr1[10] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
    double arr2[10] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
    BOOST_VERIFY(eucl_dist(&arr1[0], &arr2[5], NCOL) == sqrt(125));

    printf("Exiting test_eucl ==> ");
}

void test_dist_matrix() {
    constexpr unsigned k = 4;
    test_init_g_clusters();
    dist_matrix::ptr test_mat = dist_matrix::create(k);

    /* Test compute_dist */
    test_mat->compute_dist(g_clusters, k);

    printf("Clusters:\n"); print_clusters<prune_cluster>(g_clusters);
    printf("Cluster distance :\n"); test_mat->print();

    /* Test s_val */
    printf("Printing s_vals:\n");
    for (unsigned i = 0; i < k; i++) {
        BOOST_VERIFY(g_clusters[i]->get_s_val() ==
                test_mat->get_min_dist(i));
    }
    printf("\n");
    printf("Exiting test_dist_matrix ==> ");
}

int main(int argc, char* argv[]) {
    if (argc > 1) {
        if (std::string(argv[1]) == "-e") {
            test_eucl(); std::cout << "Test eucl Success ...\n";
        } else if (std::string(argv[1]) == "-d") {
            test_dist_matrix(); std::cout << "Test distance matrix Success ...\n";
        } else if (std::string(argv[1]) == "-c") {
            test_cluster(); std::cout << "Test cluster Success ...\n";
        } else {
            fprintf(stderr, "Unknown test option '%s'", argv[1]);
        }
    } else { /* Do all tests */
        test_eucl(); std::cout << "Test eucl Success ...\n";
        test_dist_matrix(); std::cout << "Test distance matrix Success ...\n";
        test_cluster(); std::cout << "Test cluster Success ...\n";
    }
    return EXIT_SUCCESS;

}
