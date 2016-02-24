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
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY CURRENT_KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#include "kmeans.h"
#define KM_TEST 1
#define VERBOSE 0

#ifdef PROFILER
#include <gperftools/profiler.h>
#endif

namespace {
	static unsigned NUM_COLS;
	static size_t K;
	static unsigned NUM_ROWS;
	short OMP_MAX_THREADS;
	static unsigned g_num_changed = 0;
    static const unsigned INVALID_CLUSTER_ID = std::numeric_limits<unsigned>::max();
	static struct timeval start, end;
    static init_type_t g_init_type;

	/**
	 * \brief This initializes clusters by randomly choosing sample
	 *		membership in a cluster.
	 * See: http://en.wikipedia.org/wiki/K-means_clustering#Initialization_methods
	 *	\param cluster_assignments Which cluster each sample falls into.
	 */
	static void random_partition_init(unsigned* cluster_assignments,
            const double* matrix, clusters::ptr clusters) {
		BOOST_LOG_TRIVIAL(info) << "Random init start";

		// #pragma omp parallel for firstprivate(cluster_assignments, K) shared(cluster_assignments)
		for (unsigned row = 0; row < NUM_ROWS; row++) {
			size_t asgnd_clust = random() % K; // 0...K
            clusters->add_member(&matrix[row*NUM_COLS], asgnd_clust);
			cluster_assignments[row] = asgnd_clust;
		}

		// NOTE: M-Step called in compute func to update cluster counts & centers
#if VERBOSE
		printf("After rand paritions cluster_asgns: "); print_arr(cluster_assignments, NUM_ROWS);
#endif
		BOOST_LOG_TRIVIAL(info) << "Random init end\n";
	}

	/**
	 * \brief Forgy init takes `K` random samples from the matrix
	 *		and uses them as cluster centers.
	 * \param matrix the flattened matrix who's rows are being clustered.
	 * \param clusters The cluster centers (means) flattened matrix.
	 */
	static void forgy_init(const double* matrix, clusters::ptr clusters) {

		BOOST_LOG_TRIVIAL(info) << "Forgy init start";

		for (unsigned clust_idx = 0; clust_idx < K; clust_idx++) { // 0...K
			unsigned rand_idx = random() % (NUM_ROWS - 1); // 0...(n-1)
            clusters->set_mean(&matrix[rand_idx*NUM_COLS], clust_idx);
		}

		BOOST_LOG_TRIVIAL(info) << "Forgy init end";
	}

	/**
	 * \brief A parallel version of the kmeans++ initialization alg.
	 *  See: http://ilpubs.stanford.edu:8090/778/1/2006-13.pdf for algorithm
	 */
	static void kmeanspp_init(const double* matrix, clusters::ptr clusters,
            unsigned* cluster_assignments, std::vector<double>& dist_v) {

		// Choose c1 uniiformly at random
		unsigned selected_idx = random() % NUM_ROWS; // 0...(NUM_ROWS-1)

        clusters->set_mean(&matrix[selected_idx*NUM_COLS], 0);
		dist_v[selected_idx] = 0.0;

#if KM_TEST
		BOOST_LOG_TRIVIAL(info) << "\nChoosing "
            << selected_idx << " as center K = 0";
#endif

		unsigned clust_idx = 0; // The number of clusters assigned

		// Choose next center c_i with weighted prob
		while ((clust_idx + 1) < K) {
			double cum_dist = 0;
#pragma omp parallel for reduction(+:cum_dist) shared (dist_v)
			for (size_t row = 0; row < NUM_ROWS; row++) {
                double dist = dist_comp_raw(&matrix[row*NUM_COLS],
                            &((clusters->get_means())[clust_idx*NUM_COLS]), NUM_COLS);

				if (dist < dist_v[row]) { // Found a closer cluster than before
					dist_v[row] = dist;
                    cluster_assignments[row] = clust_idx;
				}
				cum_dist += dist_v[row];
			}

			cum_dist = (cum_dist * ((double)random())) / (RAND_MAX - 1.0);
			clust_idx++;

			for (size_t i=0; i < NUM_ROWS; i++) {
				cum_dist -= dist_v[i];
				if (cum_dist <= 0) {
#if KM_TEST
					BOOST_LOG_TRIVIAL(info) << "Choosing "
                        << i << " as center K = " << clust_idx;
#endif
                    clusters->set_mean(&(matrix[i*NUM_COLS]), clust_idx);
					break;
				}
			}
			assert (cum_dist <= 0);
		}

#if VERBOSE
		BOOST_LOG_TRIVIAL(info) << "\nCluster centers after kmeans++"; clusters->print_means();
#endif
	}


	/**
	 * \brief Update the cluster assignments while recomputing distance matrix.
	 * \param matrix The flattened matrix who's rows are being clustered.
	 * \param clusters The cluster centers (means) flattened matrix.
	 *	\param cluster_assignments Which cluster each sample falls into.
	 */
	static void EM_step(const double* matrix, clusters::ptr cls,
            unsigned* cluster_assignments, unsigned* cluster_assignment_counts,
            const bool prune_init=false) {

		std::vector<clusters::ptr> pt_cl(OMP_MAX_THREADS);
        // Per thread changed cluster count. OMP_MAX_THREADS
		std::vector<size_t> pt_num_change(OMP_MAX_THREADS);

        for (int i = 0; i < OMP_MAX_THREADS; i++)
            pt_cl[i] = clusters::create(K, NUM_COLS);

#pragma omp parallel for firstprivate(matrix, pt_cl)\
        shared(cluster_assignments) schedule(static)
		for (unsigned row = 0; row < NUM_ROWS; row++) {

            size_t asgnd_clust = INVALID_CLUSTER_ID;
            double best, dist;
            dist = best = std::numeric_limits<double>::max();

            for (unsigned clust_idx = 0; clust_idx < K; clust_idx++) {
                dist = dist_comp_raw(&matrix[row*NUM_COLS],
                        &(cls->get_means()[clust_idx*NUM_COLS]), NUM_COLS);

                if (dist < best) {
                    best = dist;
                    asgnd_clust = clust_idx;
                }
            }

            BOOST_VERIFY(asgnd_clust != INVALID_CLUSTER_ID);

            if (asgnd_clust != cluster_assignments[row]) {
                pt_num_change[omp_get_thread_num()]++;
            }
            cluster_assignments[row] = asgnd_clust;
            pt_cl[omp_get_thread_num()]->add_member(&matrix[row*NUM_COLS], asgnd_clust);
            // Accumulate for local copies
		}

#if VERBOSE
		BOOST_LOG_TRIVIAL(info) << "Clearing cluster assignment counts";
		BOOST_LOG_TRIVIAL(info) << "Clearing cluster centers ...";
#endif
        cls->clear();

		// Serial aggreate of OMP_MAX_THREADS vectors
        // TODO: Pool these
		for (int thd = 0; thd < OMP_MAX_THREADS; thd++) {
			// Updated the changed cluster count
			g_num_changed += pt_num_change[thd];
			// Summation for cluster centers
            cls->peq(pt_cl[thd]);
		}

        unsigned chk_nmemb = 0;
        for (unsigned clust_idx = 0; clust_idx < K; clust_idx++) {
            cls->finalize(clust_idx);
            cluster_assignment_counts[clust_idx] = cls->get_num_members(clust_idx);
            chk_nmemb += cluster_assignment_counts[clust_idx];
        }
        BOOST_VERIFY(chk_nmemb == NUM_ROWS);

#if KM_TEST
		BOOST_LOG_TRIVIAL(info) << "Global number of changes: " << g_num_changed;
#endif
	}
}

namespace fg
{
	unsigned compute_kmeans(const double* matrix, double* clusters_ptr,
			unsigned* cluster_assignments, unsigned* cluster_assignment_counts,
			const unsigned num_rows, const unsigned num_cols, const unsigned k,
            const unsigned MAX_ITERS, const int max_threads, const std::string init,
            const double tolerance, const std::string dist_type)
	{
#ifdef PROFILER
		ProfilerStart("matrix/kmeans.perf");
#endif
		NUM_COLS = num_cols;
		K = k;
		NUM_ROWS = num_rows;
		assert(max_threads > 0);

		OMP_MAX_THREADS = std::min(max_threads, get_num_omp_threads());
		omp_set_num_threads(OMP_MAX_THREADS);
		BOOST_LOG_TRIVIAL(info) << "Running on " << OMP_MAX_THREADS << " threads!";

        // Check k
		if (K > NUM_ROWS || K < 2 || K == (unsigned)-1) {
			BOOST_LOG_TRIVIAL(fatal)
				<< "'k' must be between 2 and the number of rows in the matrix" <<
				"k = " << K;
			exit(-1);
		}

		gettimeofday(&start , NULL);
		/*** Begin VarInit of data structures ***/
        std::fill(&cluster_assignments[0], (&cluster_assignments[0])+NUM_ROWS, -1);
        std::fill(&cluster_assignment_counts[0], (&cluster_assignment_counts[0])+K, 0);

        clusters::ptr clusters = clusters::create(K, NUM_COLS);

        if (init == "none")
            clusters->set_mean(clusters_ptr);

        // For pruning
        std::vector<double> dist_v;
        dist_v.assign(NUM_ROWS, std::numeric_limits<double>::max());

		/*** End VarInit ***/
        BOOST_LOG_TRIVIAL(info) << "Dist_type is " << dist_type;
        if (dist_type == "eucl") {
            g_dist_type = EUCL;
        } else if (dist_type == "cos") {
            g_dist_type = COS;
        } else {
            BOOST_LOG_TRIVIAL(fatal)
                << "[ERROR]: param dist_type must be one of: 'eucl', 'cos'.It is '"
                << dist_type << "'";
            exit(-1);
        }

        if (init == "random") {
            random_partition_init(cluster_assignments, matrix, clusters);
            g_init_type = RANDOM;
        } else if (init == "forgy") {
            forgy_init(matrix, clusters);
            g_init_type = FORGY;
        } else if (init == "kmeanspp") {
            kmeanspp_init(matrix, clusters, cluster_assignments, dist_v);
            g_init_type = PLUSPLUS;
        } else if (init == "none") {
            g_init_type = NONE;
        } else {
            BOOST_LOG_TRIVIAL(fatal)
                << "[ERROR]: param init must be one of: "
                "'random', 'forgy', 'kmeanspp'.It is '"
                << init << "'";
            exit(-1);
        }

        if (g_init_type == NONE || g_init_type == FORGY
                || g_init_type == RANDOM) {
			EM_step(matrix, clusters, cluster_assignments,
                    cluster_assignment_counts, true);
        }

        BOOST_LOG_TRIVIAL(info) << "Init is '" << init << "'";
		BOOST_LOG_TRIVIAL(info) << "Matrix K-means starting ...";

		bool converged = false;
		std::string str_iters = MAX_ITERS == std::numeric_limits<unsigned>::max() ?
			"until convergence ...":
			std::to_string(MAX_ITERS) + " iterations ...";
		BOOST_LOG_TRIVIAL(info) << "Computing " << str_iters;
		unsigned iter = 1;

		while (iter < MAX_ITERS) {
			// Hold cluster assignment counter
			BOOST_LOG_TRIVIAL(info) << "E-step Iteration " << iter <<
				". Computing cluster assignments ...";
			EM_step(matrix, clusters, cluster_assignments,
                    cluster_assignment_counts);
#if KM_TEST
            printf("Cluster assignment counts: ");
            print_arr(cluster_assignment_counts, K);
#endif

			if (g_num_changed == 0 || ((g_num_changed/(double)NUM_ROWS))
                    <= tolerance) {
				converged = true;
				break;
			} else {
				g_num_changed = 0;
			}
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

		return iter;
	}
}
