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
#ifndef __KMEANS_COORDINATOR__
#define __KMEANS_COORDINATOR__

#include <vector>
#include <unordered_map>
#include <memory>

#include "kmeans_thread.h"
#include "libgraph-algs/sem_kmeans_util.h"

#ifdef PROFILER
#include <gperftools/profiler.h>
#endif

namespace {
    typedef std::vector<kmeans_thread::ptr>::iterator thread_iter;
#if 0
    double* g_data; // TEST
#endif

    class kmeans_coordinator {
        private:
            unsigned nthreads, nnodes;
            size_t nrow, ncol;
            std::vector<kmeans_thread::ptr> threads;
            std::string fn; // File on disk
            unsigned* cluster_assignments;
            unsigned* cluster_assignment_counts;
            unsigned k;
            clusters::ptr cltrs;
            init_type_t _init_t;
            dist_type_t _dist_t;
            double tolerance;
            unsigned max_iters;
            unsigned num_changed; // total # of samples that changed clstr in an iter

            // Metadata
            std::vector<unsigned> thd_max_row_idx; // max index stored within each threads partition

            kmeans_coordinator(const std::string fn, const size_t nrow,
                    const size_t ncol, const unsigned k, const unsigned max_iters,
                    const unsigned nnodes, const unsigned nthreads,
                    const double* centers, const init_type_t it,
                    const double tolerance, const dist_type_t dt) {
#if 0
                // TEST //
                std:: cout << "Reading " << nrow <<
                    "rows, " << ncol << "cols\n";
                bin_reader<double> b(fn, nrow, ncol);
                //g_data.resize(nrow*ncol);
                g_data = new double[nrow*ncol];
                b.read(g_data);
#endif

                this->fn = fn;
                this->nrow = nrow;
                this->ncol = ncol;
                this->k = k;
                BOOST_ASSERT_MSG(k >= 1, "[FATAL]: 'k' must be >= 1");
                this->max_iters = max_iters;
                this->nnodes = nnodes;
                this->nthreads = nthreads;
                if (nthreads >  (unsigned)get_num_omp_threads()) {
                    BOOST_LOG_TRIVIAL(warning) << "[WARNING]: Exceeded system"
                        " #virtual cores of: " << get_num_omp_threads();
                }
                this->_init_t = it;
                this->tolerance = tolerance;
                this->_dist_t = dt;
                num_changed = 0;

                BOOST_VERIFY(cluster_assignments = new unsigned [nrow]);
                BOOST_VERIFY(cluster_assignment_counts = new unsigned [k]);

                std::fill(&cluster_assignments[0],
                        (&cluster_assignments[0])+nrow, -1);
                std::fill(&cluster_assignment_counts[0],
                        (&cluster_assignment_counts[0])+k, 0);

                cltrs = clusters::create(k, ncol);
                if (centers) {
                    cltrs->set_mean(centers);
                    if (_init_t != NONE) {
                        BOOST_LOG_TRIVIAL(warning) << "[WARNING]: Init method " <<
                            "ignored because centers provided!";
                    } else {
                        BOOST_LOG_TRIVIAL(info) << "Init-ed centers ...";
                    }
                }
                // NUMA node affinity binding policy is round-robin
                unsigned thds_row = nrow / nthreads;
                for (unsigned thd_id = 0; thd_id < nthreads; thd_id++) {
                    std::pair<size_t, unsigned> tup = get_offset_len_tup(thd_id);
                    thd_max_row_idx.push_back((thd_id*thds_row) + tup.second);
                    threads.push_back(kmeans_thread::create((thd_id % nnodes),
                                thd_id, tup.first, tup.second, thds_row,
                                ncol, cltrs, cluster_assignments, fn));
                }
            }

        public:
            typedef std::shared_ptr<kmeans_coordinator> ptr;

            static ptr create(const std::string fn, const size_t nrow,
                    const size_t ncol, const unsigned k, const unsigned max_iters,
                    const unsigned nnodes, const unsigned nthreads,
                    const double* centers=NULL, const std::string init="kmeanspp",
                    const double tolerance=-1, const std::string dist_type="eucl") {

                init_type_t _init_t;
                if (init == "random")
                    _init_t = RANDOM;
                else if (init == "forgy")
                    _init_t = FORGY;
                else if (init == "kmeanspp")
                    _init_t = PLUSPLUS;
                else if (init == "none")
                    _init_t = NONE;
                else {
                    BOOST_LOG_TRIVIAL(fatal) << "[ERROR]: param init must be one of:"
                       " [random | forgy | kmeanspp]. It is '" << init << "'";
                    exit(-1);
                }

                dist_type_t _dist_t;
                if (dist_type == "eucl")
                    _dist_t = EUCL;
                else if (dist_type == "cos")
                    _dist_t = COS;
                else {
                    BOOST_LOG_TRIVIAL(fatal) << "[ERROR]: param dist_type must be one of:"
                       " 'eucl', 'cos'.It is '" << dist_type << "'";
                    exit(-1);
                }
#if KM_TEST
                printf("kmeans coordinator => NUMA nodes: %u, nthreads: %u, "
                        "nrow: %lu, ncol: %lu, init: '%s', dist_t: '%s', fn: '%s'"
                        "\n\n", nnodes, nthreads, nrow, ncol, init.c_str(),
                        dist_type.c_str(), fn.c_str());
#endif
                return ptr(new kmeans_coordinator(fn, nrow, ncol, k, max_iters,
                            nnodes, nthreads, centers, _init_t, tolerance, _dist_t));
            }

            std::pair<size_t, unsigned> get_offset_len_tup(const unsigned thd_id) {
                unsigned rows_per_thread = nrow / nthreads;
                size_t start_offset = (thd_id*rows_per_thread*ncol);

                if (thd_id == nthreads - 1)
                    rows_per_thread += nrow % nthreads;
                return std::pair<size_t, unsigned>(start_offset, rows_per_thread);
            }

            // Pass file handle to threads to read & numa alloc
            void numa_alloc_data();
            void join_threads();
            void create_thread_map();
            void run_init();
            void run_kmeans();
            void run_min_tri_kmeans();
            void random_partition_init();
            void forgy_init();
            void update_clusters();
            void kmeanspp_init();
            void start_threads(const thread_state_t state);
            void set_thread_clust_idx(const unsigned clust_idx);
            double reduction_on_cuml_sum();
            void set_thd_dist_v_ptr(double* v);

            const double* get_thd_data(const unsigned row_id) const;

            ~kmeans_coordinator() {
                std::vector<kmeans_thread::ptr>::iterator it = threads.begin();
                for (; it != threads.end(); ++it)
                    (*it)->destroy_numa_mem();
                delete [] cluster_assignments;
                delete [] cluster_assignment_counts;

#if 0
                delete [] g_data;
#endif
            }
    };

    typedef std::pair<unsigned, unsigned> thd_row_tup;

    // Pass file handle to threads to read & numa alloc
    void kmeans_coordinator::numa_alloc_data() {
        for (thread_iter it = threads.begin(); it != threads.end(); ++it)
            (*it)->start(ALLOC_DATA);
        join_threads();
    }

    void kmeans_coordinator::join_threads() {
        for (thread_iter it = threads.begin(); it != threads.end(); ++it)
            (*it)->join();
    }

    // <Thread, within-thread-row-id>
    const double* kmeans_coordinator::get_thd_data(const unsigned row_id) const {
        // TODO: Cheapen
        unsigned parent_thd = std::upper_bound(thd_max_row_idx.begin(),
                thd_max_row_idx.end(), row_id) - thd_max_row_idx.begin();
        unsigned rows_per_thread = nrow/nthreads; // All but the last thread

#if 0
        printf("Global row %u, row in parent thd: %u --> %u\n", row_id,
                parent_thd, (row_id-(parent_thd*rows_per_thread)));
#endif

        return &((threads[parent_thd]->get_local_data())
                [(row_id-(parent_thd*rows_per_thread))*ncol]);
    }

    void kmeans_coordinator::random_partition_init() {
        printf("In random init testing data ...\n");
        for (unsigned row = 0; row < nrow; row++) {
            unsigned asgnd_clust = random() % k; // 0...k

            const double* dp = get_thd_data(row);
#if 0
            if (!(eq_all(dp, &g_data[row*ncol], ncol))) {
                printf("Row: %u, Correct data: ", row);
                print_arr<double>(&g_data[row*ncol], ncol);
                printf("Retrived data: ");
                print_arr<double>(dp, ncol);
                assert(0);
            }
#endif
            cltrs->add_member(dp, asgnd_clust);
            cluster_assignments[row] = asgnd_clust;
        }

        for (unsigned cl = 0; cl < k; cl++)
            cltrs->finalize(cl);

#if VERBOSE
        printf("After rand paritions cluster_asgns: ");
        print_arr<unsigned>(cluster_assignments, nrow);
#endif
    }

    void kmeans_coordinator::forgy_init() {
        BOOST_LOG_TRIVIAL(info) << "Forgy init start";
        for (unsigned clust_idx = 0; clust_idx < k; clust_idx++) { // 0...k
            unsigned rand_idx = random() % (nrow - 1); // 0...(nrow-1)
            cltrs->set_mean(get_thd_data(rand_idx), clust_idx);
        }
        BOOST_LOG_TRIVIAL(info) << "Forgy init end";
    }

    void kmeans_coordinator::run_init() {
        switch(_init_t) {
            case RANDOM:
                random_partition_init();
                break;
            case FORGY:
                forgy_init();
                break;
            case PLUSPLUS:
                kmeanspp_init();
                break;
            case NONE:
                break;
            default:
                fprintf(stderr, "[FATAL]: Unknow initialization type\n");
                exit(EXIT_FAILURE);
        }
    }

    void kmeans_coordinator::update_clusters() {
        num_changed = 0; // Always reset here since there's no pruning
        cltrs->clear();

        // Serial aggreate of OMP_MAX_THREADS vectors
        for (thread_iter it = threads.begin(); it != threads.end(); ++it) {
            // Updated the changed cluster count
            num_changed += (*it)->get_num_changed();
            // Summation for cluster centers

#if VERBOSE
            printf("Thread %ld clusters:\n", (it-threads.begin()));
            ((*it)->get_local_clusters())->print_means();
#endif

            cltrs->peq((*it)->get_local_clusters());
        }

        unsigned chk_nmemb = 0;
        for (unsigned clust_idx = 0; clust_idx < k; clust_idx++) {
            cltrs->finalize(clust_idx);
            cluster_assignment_counts[clust_idx] =
                cltrs->get_num_members(clust_idx);
            chk_nmemb += cluster_assignment_counts[clust_idx];
        }
        if (chk_nmemb != nrow)
            printf("chk_nmemb = %u\n", chk_nmemb);

        BOOST_VERIFY(chk_nmemb == nrow);
        BOOST_VERIFY(num_changed <= nrow);

#if KM_TEST
        BOOST_LOG_TRIVIAL(info) << "Global number of changes: " << num_changed;
#endif
    }

    void kmeans_coordinator::run_min_tri_kmeans() {
        assert(0); // TODO: min-tri-kmeans logic here
    }

    double kmeans_coordinator::reduction_on_cuml_sum() {
        double tot = 0;
        for (thread_iter it = threads.begin(); it != threads.end(); ++it)
            tot += (*it)->get_culm_dist();
        return tot;
    }

    void kmeans_coordinator::set_thread_clust_idx(const unsigned clust_idx) {
        for (thread_iter it = threads.begin(); it != threads.end(); ++it)
            (*it)->set_clust_idx(clust_idx);
    }

    void kmeans_coordinator::set_thd_dist_v_ptr(double* v) {
        for (thread_iter it = threads.begin(); it != threads.end(); ++it)
            (*it)->set_dist_v_ptr(v);
    }

    void kmeans_coordinator::start_threads(const thread_state_t state) {
        for (thread_iter it = threads.begin(); it != threads.end(); ++it) {
            (*it)->start(state);
        }
    }

    void kmeans_coordinator::kmeanspp_init() {
        struct timeval start, end;
        gettimeofday(&start , NULL);

        std::vector<double> dist_v;
        dist_v.assign(nrow, std::numeric_limits<double>::max());
        set_thd_dist_v_ptr(&dist_v[0]);

		// Choose c1 uniformly at random
		unsigned selected_idx = random() % nrow; // 0...(nrow-1)

        cltrs->set_mean(get_thd_data(selected_idx), 0);
		dist_v[selected_idx] = 0.0;
#if KM_TEST
		BOOST_LOG_TRIVIAL(info) << "\nChoosing "
            << selected_idx << " as center K = 0";
#endif
		unsigned clust_idx = 0; // The number of clusters assigned

		// Choose next center c_i with weighted prob
		while ((clust_idx + 1) < k) {
            set_thread_clust_idx(clust_idx); // Set the current cluster index
            start_threads(KMSPP_INIT); // Run || distance comp to clust_idx
            join_threads();
			double cuml_dist = reduction_on_cuml_sum(); // Sum the per thread cumulative dists

			cuml_dist = (cuml_dist * ((double)random())) / (RAND_MAX - 1.0);
			clust_idx++;

			for (size_t row = 0; row < nrow; row++) {
				cuml_dist -= dist_v[row];
				if (cuml_dist <= 0) {
#if KM_TEST
					BOOST_LOG_TRIVIAL(info) << "Choosing "
                        << row << " as center k = " << clust_idx;
#endif
                    cltrs->set_mean(get_thd_data(row), clust_idx);
					break;
				}
			}
			BOOST_VERIFY(cuml_dist <= 0);
		}

#if VERBOSE
		BOOST_LOG_TRIVIAL(info) << "\nCluster centers after kmeans++";
        clusters->print_means();
#endif
        gettimeofday(&end, NULL);
        BOOST_LOG_TRIVIAL(info) << "\n\nInitialization time: " <<
            time_diff(start, end) << " sec\n";
    }

    /**
     * Main driver for kmeans
     */
    void kmeans_coordinator::run_kmeans() {
#ifdef PROFILER
		ProfilerStart("matrix/kmeans_coordinator.perf");
#endif
        numa_alloc_data(); // Move data
        struct timeval start, end;
        gettimeofday(&start , NULL);
        run_init(); // Initialize clusters

        // Run kmeans loop
        size_t iter = 1;
        bool converged = false;
        while (iter <= max_iters) {
            BOOST_LOG_TRIVIAL(info) << "E-step Iteration: " << iter;
            start_threads(EM);
            join_threads();
            update_clusters();

            printf("Cluster assignment counts: ");
            print_arr(cluster_assignment_counts, k);

            if (num_changed == 0 ||
                    ((num_changed/(double)nrow)) <= tolerance) {
                converged = true;
                break;
            }
            iter++;
        }
#ifdef PROFILER
		ProfilerStop();
#endif
        gettimeofday(&end, NULL);
        BOOST_LOG_TRIVIAL(info) << "\n\nAlgorithmic time taken = " <<
            time_diff(start, end) << " sec\n";

        BOOST_LOG_TRIVIAL(info) << "\n******************************************\n";
        if (converged) {
            BOOST_LOG_TRIVIAL(info) <<
                "K-means converged in " << iter << " iterations";
        } else {
            BOOST_LOG_TRIVIAL(warning) << "[Warning]: K-means failed to converge in "
                << iter << " iterations";
        }
        BOOST_LOG_TRIVIAL(info) << "\n******************************************\n";
    }
}
#endif
