/*
 * Copyright 2015 Open Connectome Project (http://openconnecto.me)
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
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY CURRENT_KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#include <signal.h>
#ifdef PROFILER
#include <gperftools/profiler.h>
#endif

#include <vector>
#include <algorithm>

#include "graph_engine.h"
#include "graph_config.h"
#include "FGlib.h"

#define KM_TEST 0

using namespace fg;

namespace {

    typedef safs::page_byte_array::seq_const_iterator<edge_count> data_seq_iterator;
    static unsigned NUM_COLS;
    static unsigned NUM_ROWS;
    static unsigned K;
    static unsigned g_num_changed = 0;
    static struct timeval start, end;

    enum dist_type_t { EUCL, COS }; // Euclidean, Cosine distance

    template <typename T>
        static void print_arr(T* arr, unsigned len) {
            printf("[ ");
            for (unsigned i = 0; i < len; i++) {
                std::cout << arr[i] << " ";
            }
            printf("]\n");
        }

    template <typename T>
        static void print_vector(typename std::vector<T>& v) {
            std::cout << "[";
            typename std::vector<T>::iterator itr = v.begin();
            for (; itr != v.end(); itr++) {
                std::cout << " "<< *itr;
            }
            std::cout <<  " ]\n";
        }

    template <typename T>
        static void print_mat(T* matrix, const unsigned rows, const unsigned cols) {
            for (unsigned row = 0; row < rows; row++) {
                std::cout << "[";
                for (unsigned col = 0; col < cols; col++) {
                    std::cout << " " << matrix[row*cols + col];
                }
                std::cout <<  " ]\n";
            }
        }

    class cluster {
        private:
            std::vector<double> mean; // cluster mean
            unsigned num_members; // cluster assignment count
            bool complete; // Have we already divided by num_members

            void div(const unsigned val) {
                for (unsigned i = 0; i < mean.size(); i++) {
                    mean[i] /= val;
                }
            }

            cluster(const unsigned len) {
                mean.assign(len, 0);
                num_members = 0;
                complete = false;
            }

        public:
            typedef typename std::shared_ptr<cluster> ptr;

            static ptr create(const unsigned len) {
                return ptr(new cluster(len));
            }

            void init(const unsigned len) {
                mean.assign(len, 0);
                num_members = 0;
            }

            void assign(const double val) {
                mean.assign(mean.size(), val);
            }

            void clear() {
                this->assign(0);
                this->num_members = 0;
            }

            const std::vector<double>& get_mean() const {
                return mean;
            }

            const unsigned get_num_members() const {
                return num_members;
            }

            cluster& operator=(const cluster& other) {
                this->mean = other.get_mean();
                this->num_members = other.get_num_members();
                return *this;
            }
    };
}


namespace fg
{
    void compute_sem_kmeans(FG_graph::ptr fg, const size_t k, const std::string init,
            const unsigned MAX_ITERS, const double tolerance)
    {
#ifdef PROFILER
        ProfilerStart("/home/disa/FlashGraph/flash-graph/matrix/kmeans.perf");
#endif
        K = k;
        std::vector<cluster::ptr> clusters; // cluster means/centers
        std::vector<unsigned> cluster_assignments; // Which cluster a sample is in
        NUM_ROWS = 0; // FIXME: Get from graph or header
        NUM_COLS = 0; // FIXME: Store the maximum dimension of the matrix in header

        // Check k
        if (K > NUM_ROWS || K < 2 || K == (unsigned)-1) {
            BOOST_LOG_TRIVIAL(fatal)
                << "'k' must be between 2 and the number of rows in the matrix" <<
                "k = " << K;
            exit(-1);
        }

        // Check Initialization
        if (init.compare("random") && init.compare("kmeanspp") &&
                init.compare("forgy")) {
            BOOST_LOG_TRIVIAL(fatal)
                << "[ERROR]: param init must be one of: 'random', 'forgy', 'kmeanspp'.It is '"
                << init << "'";
            exit(-1);
        }

        gettimeofday(&start , NULL);
        /*** Begin VarInit of data structures ***/
        cluster_assignments.assign(NUM_ROWS, -1);
        for (size_t cl = 0; cl < k; cl++)
            clusters.push_back(cluster::create(NUM_COLS));

        /*** End VarInit ***/

        if (init == "random") {
            BOOST_LOG_TRIVIAL(info) << "Init is '"<< init <<"'";
            // random_partition_init(cluster_assignments); // TODO
            // We must now update cluster centers before we begin
            //M_step();
        }
        if (init == "forgy") {
            BOOST_LOG_TRIVIAL(info) << "Init is '"<< init <<"'";
            //forgy_init(); // TODO
        }
        else {
            if (init == "kmeanspp") {
            BOOST_LOG_TRIVIAL(info) << "Init is '"<< init <<"'";
                //kmeanspp_init(matrix, clusters); // TODO
            }
        }

        BOOST_LOG_TRIVIAL(info) << "SEM-K||means starting ...";

        bool converged = false;

        std::string str_iters = MAX_ITERS == std::numeric_limits<unsigned>::max() ?
            "until convergence ...":
            std::to_string(MAX_ITERS) + " iterations ...";
        BOOST_LOG_TRIVIAL(info) << "Computing " << str_iters;
        unsigned iter = 1;
        while (iter < MAX_ITERS) {
            BOOST_LOG_TRIVIAL(info) << "E-step Iteration " << iter <<
                " . Computing cluster assignments ...";
            // E_step();

            if (g_num_changed == 0 || ((g_num_changed/(double)NUM_ROWS)) <= tolerance) {
                converged = true;
                break;
            } else {
                g_num_changed = 0;
            }

#if KM_TEST
            BOOST_LOG_TRIVIAL(info) << "M-step Updating cluster means ...";
#endif
            // M_step();
            iter++;
        }

        gettimeofday(&end, NULL);
        BOOST_LOG_TRIVIAL(info) << "\n\nAlgorithmic time taken = " <<
            time_diff(start, end) << " sec\n";

#ifdef PROFILER
        ProfilerStop();
#endif
        BOOST_LOG_TRIVIAL(info) << "\n******************************************\n";

        if (converged) {
            BOOST_LOG_TRIVIAL(info) <<
                "K-means converged in " << iter << " iterations";
        } else {
            BOOST_LOG_TRIVIAL(warning) << "[Warning]: K-means failed to converge in "
                << iter << " iterations";
        }
        BOOST_LOG_TRIVIAL(info) << "\n******************************************\n";

        //return;
    }
}
