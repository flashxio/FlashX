/*
 * Copyright 2014 Open Connectome Project (http://openconnecto.me)
 * Written by Disa Mhembere (disa@jhu.edu)
 *
 * This file is part of FlashMatrix.
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
    static dist_type_t g_dist_type;

	static struct timeval start, end;

    // TODO: Doc
    static double const eucl_dist(const unsigned row, const unsigned clust_idx,
            const double* matrix, const double* clusters) {
        double dist = 0;
        for (unsigned col = 0; col < NUM_COLS; col++) {
            double diff = matrix[(row*NUM_COLS) + col] -
                clusters[clust_idx*NUM_COLS + col];
            dist += diff * diff;
        }
        return sqrt(dist); // TODO: rm sqrt
    }

    // TODO: Doc
    double const cos_dist(const unsigned row, const unsigned clust_idx,
            const double* matrix, const double* clusters) {
        double numr, ldenom, rdenom;
        numr = ldenom = rdenom = 0;

        for (unsigned col = 0; col < NUM_COLS; col++) {
            double a = matrix[(row*NUM_COLS) + col];
            double b = clusters[clust_idx*NUM_COLS + col];

            numr += a*b;
            ldenom += a*a;
            rdenom += b*b;
        }
        return  1 - (numr / ((sqrt(ldenom)*sqrt(rdenom))));
    }

    static size_t get_best_cluster(const unsigned row, const double* matrix, const double* clusters) {
        size_t asgnd_clust = INVALID_CLUSTER_ID;
        double best = std::numeric_limits<double>::max();
        double dist = -1;

        for (unsigned clust_idx = 0; clust_idx < K; clust_idx++) {
            if (g_dist_type == EUCL)
                dist = eucl_dist(row, clust_idx, matrix, clusters);
            else if (g_dist_type == COS)
                dist = cos_dist(row, clust_idx, matrix, clusters);
            else
                BOOST_ASSERT_MSG(false, "Unknown distance metric!");

            if (dist < best) {
                best = dist;
                asgnd_clust = clust_idx;
            }
        }
        return asgnd_clust;
    }

	/**
	 * \brief This initializes clusters by randomly choosing sample
	 *		membership in a cluster.
	 * See: http://en.wikipedia.org/wiki/K-means_clustering#Initialization_methods
     * \param cluster_assignments Which cluster each sample falls into.
	 */
	static void random_partition_init(unsigned* cluster_assignments) {
		BOOST_LOG_TRIVIAL(info) << "Random init start";

		// #pragma omp parallel for firstprivate(cluster_assignments, K) shared(cluster_assignments)
		for (unsigned row = 0; row < NUM_ROWS; row++) {
			size_t assigned_cluster = random() % K; // 0...K
			cluster_assignments[row] = assigned_cluster;
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
	static void forgy_init(const double* matrix, double* clusters) {

		BOOST_LOG_TRIVIAL(info) << "Forgy init start";

		for (unsigned clust_idx = 0; clust_idx < K; clust_idx++) { // 0...K
			unsigned rand_idx = random() % (NUM_ROWS - 1); // 0...(n-1)
            std::copy(&matrix[rand_idx*NUM_COLS],
                    &matrix[(rand_idx*NUM_COLS)+NUM_COLS],
                    &clusters[clust_idx*NUM_COLS]);
		}

		BOOST_LOG_TRIVIAL(info) << "Forgy init end";
	}

	/**
	 * \brief A parallel version of the kmeans++ initialization alg.
	 *  See: http://ilpubs.stanford.edu:8090/778/1/2006-13.pdf for algorithm
	 */
	static void kmeanspp_init(const double* matrix, double* clusters) {
		// Init distance array
        std::vector<double> dist_v = std::vector<double>(NUM_ROWS);
        std::fill(dist_v.begin(), dist_v.end(), std::numeric_limits<double>::max());

		// Choose c1 uniiformly at random
		unsigned rand_idx = random() % NUM_ROWS; // 0...(NUM_ROWS-1)

        std::copy(&matrix[rand_idx*NUM_COLS],
                &matrix[(rand_idx*NUM_COLS)+NUM_COLS], &clusters[0]);
		dist_v[rand_idx] = 0.0;

#if KM_TEST
		BOOST_LOG_TRIVIAL(info) << "\nChoosing " << rand_idx << " as center K = 0";
#endif

		unsigned clust_idx = 0; // The number of clusters assigned

		// Choose next center c_i with weighted prob
		while ((clust_idx + 1) < K) {
			double cum_dist = 0;
#pragma omp parallel for reduction(+:cum_dist) shared (dist_v)
			for (size_t row = 0; row < NUM_ROWS; row++) {

                double dist = -1;
                if (g_dist_type == EUCL)
                    // Per thread instance
                    dist = eucl_dist(row, clust_idx, matrix, clusters);
                else if (g_dist_type == COS)
                    dist = cos_dist(row, clust_idx, matrix, clusters);
                else
                    BOOST_ASSERT_MSG(false, "Unknown distance metric!");

				if (dist < dist_v[row]) { // Found a closer cluster than before
					dist_v[row] = dist;
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
                    std::copy(&(matrix[i*NUM_COLS]), &(matrix[i*(NUM_COLS)+NUM_COLS]),
                            &(clusters[clust_idx*NUM_COLS]));
					break;
				}
			}
			assert (cum_dist <= 0);
		}

#if VERBOSE
		BOOST_LOG_TRIVIAL(info) << "\n\nOriginal matrix: "; print_mat(matrix, NUM_ROWS, NUM_COLS);
		BOOST_LOG_TRIVIAL(info) << "\nCluster centers after kmeans++"; print_mat(clusters, K, NUM_COLS);
#endif
	}

	/**
	 * \brief Update the cluster assignments while recomputing distance matrix.
	 * \param matrix The flattened matrix who's rows are being clustered.
	 * \param clusters The cluster centers (means) flattened matrix.
	 *	\param cluster_assignments Which cluster each sample falls into.
	 */
	static void E_step(const double* matrix, double* clusters,
			unsigned* cluster_assignments, unsigned* cluster_assignment_counts) {
		// Create per thread vectors
		std::vector<std::vector<unsigned>> pt_cl_as_cnt(OMP_MAX_THREADS); // K * OMP_MAX_THREADS
		std::vector<std::vector<double>> pt_cl(OMP_MAX_THREADS); // K * NUM_COLS * OMP_MAX_THREADS
		std::vector<size_t> pt_num_change(OMP_MAX_THREADS); // Per thread changed cluster count. OMP_MAX_THREADS

		for (int i = 0; i < OMP_MAX_THREADS; i++) {
			pt_cl_as_cnt[i].resize(K); // C++ default is 0 value in all
			pt_cl[i].resize(K*NUM_COLS);
		}

#pragma omp parallel for firstprivate(matrix, clusters) shared(pt_cl_as_cnt, cluster_assignments)\
		schedule(static)
		for (unsigned row = 0; row < NUM_ROWS; row++) {
			size_t asgnd_clust = get_best_cluster(row, matrix, clusters);

            assert(asgnd_clust != INVALID_CLUSTER_ID); // TODO: RM

			if (asgnd_clust != cluster_assignments[row]) {
				pt_num_change[omp_get_thread_num()]++;
			}
			cluster_assignments[row] = asgnd_clust;

			pt_cl_as_cnt[omp_get_thread_num()][asgnd_clust]++; // Add to local copy
			// Accumulate for local copies
#if KM_TEST
			assert(omp_get_thread_num() <= OMP_MAX_THREADS);
#endif
			for (unsigned col = 0; col < NUM_COLS; col++) {
				pt_cl[omp_get_thread_num()][asgnd_clust*NUM_COLS + col] =
					pt_cl[omp_get_thread_num()][asgnd_clust*NUM_COLS + col] + matrix[row*NUM_COLS + col];
			}
		}

#if VERBOSE
		BOOST_LOG_TRIVIAL(info) << "Clearing cluster assignment counts";
		BOOST_LOG_TRIVIAL(info) << "Clearing cluster centers ...";
#endif
        std::fill(&cluster_assignment_counts[0], (&cluster_assignment_counts[0])+K, 0);
        std::fill(&clusters[0], (&clusters[0])+(K*NUM_COLS), 0);

		// Serial aggreate of OMP_MAX_THREADS vectors
		for (int i = 0; i < OMP_MAX_THREADS; i++) {
			// Updated the changed cluster count
			g_num_changed += pt_num_change[i];

			// Counts
			for (unsigned j = 0; j < K; j++) {
				cluster_assignment_counts[j] += pt_cl_as_cnt[i][j];
			}

			// Summation for cluster centers
#pragma omp parallel for firstprivate(pt_cl) shared(clusters) schedule(static)
			for (unsigned row = 0; row < K; row++) { /* ClusterID */
				for (unsigned col = 0; col < NUM_COLS; col++) { /* Features */
					clusters[row*NUM_COLS+col] = pt_cl[i][row*NUM_COLS + col] + clusters[row*NUM_COLS + col];
				}
			}
		}

#if KM_TEST
		printf("Cluster assignment counts: "); print_arr(cluster_assignment_counts, K);
		//printf("Cluster assignments:\n"); print_arr(cluster_assignments, NUM_ROWS);
		BOOST_LOG_TRIVIAL(info) << "Global number of changes: " << g_num_changed;
#endif
	}

	/**
	 * \brief Update the cluster means
	 * \param matrix The matrix who's rows are being clustered.
	 * \param clusters The cluster centers (means).
	 * \param cluster_assignment_counts How many members each cluster has.
	 * \param cluster_assignments Which cluster each sample falls into.
	 */
	static void M_step(const double* matrix, double* clusters, unsigned*
			cluster_assignment_counts, const unsigned* cluster_assignments, bool init=false) {

#if KM_TEST
		BOOST_LOG_TRIVIAL(info) << "M_step start";
#endif

		if (init) {
			BOOST_LOG_TRIVIAL(info) << "Clearing cluster centers ...";
            std::fill(&clusters[0], (&clusters[0])+(K*NUM_COLS), 0);

			BOOST_LOG_TRIVIAL(info) << "Clearing cluster assignment counts";
            std::fill(&cluster_assignment_counts[0], (&cluster_assignment_counts[0])+K, 0);

            // FIXME: Could be loop over threads not rows if we replicate clusters
			for (unsigned row = 0; row < NUM_ROWS; row++) {

				unsigned asgnd_clust = cluster_assignments[row];
				cluster_assignment_counts[asgnd_clust]++;

				for (unsigned col = 0; col < NUM_COLS; col++) {
					clusters[asgnd_clust*NUM_COLS + col] =
						clusters[asgnd_clust*NUM_COLS + col] + matrix[row*NUM_COLS + col];
				}
			}
		}

		// Take the mean of all added means
#if KM_TEST
		BOOST_LOG_TRIVIAL(info) << " Div by in place ...";
		//double* clusters = __builtin_assume_aligned(clusters);
#endif
		for (unsigned clust_idx = 0; clust_idx < K; clust_idx++) {
			if (cluster_assignment_counts[clust_idx] > 0) { // Avoid div by 0
				// #pragma omp parallel for firstprivate(cluster_assignment_counts, NUM_COLS) shared(clusters)
				for (unsigned col = 0; col < NUM_COLS; col++) {
					clusters[clust_idx*NUM_COLS + col] =
						clusters[clust_idx*NUM_COLS + col] / cluster_assignment_counts[clust_idx];
				}
			}
		}

# if 0
		printf("Cluster centers: \n"); print_mat(clusters, K, NUM_COLS);
#endif
	}
}

namespace fg
{
	unsigned compute_kmeans(const double* matrix, double* clusters,
			unsigned* cluster_assignments, unsigned* cluster_assignment_counts,
			const unsigned num_rows, const unsigned num_cols, const size_t k, const unsigned MAX_ITERS,
			const int max_threads, const std::string init, const double tolerance,
            const std::string dist_type)
	{
#ifdef PROFILER
		ProfilerStart("/home/disa/FlashGraph/flash-graph/matrix/kmeans.perf");
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
			exit(EXIT_FAILURE);
		}

		gettimeofday(&start , NULL);
		/*** Begin VarInit of data structures ***/
        std::fill(&cluster_assignments[0], (&cluster_assignments[0])+NUM_ROWS, -1);
        std::fill(&cluster_assignment_counts[0], (&cluster_assignment_counts[0])+K, 0);

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
            exit(EXIT_FAILURE);
        }

        if (init == "random") {
            random_partition_init(cluster_assignments);
            // We must now update cluster centers before we begin
            M_step(matrix, clusters, cluster_assignment_counts,
                    cluster_assignments, true);
        } else if (init == "forgy") {
            forgy_init(matrix, clusters);
        } else if (init == "kmeanspp") {
            kmeanspp_init(matrix, clusters);
        } else if (init == "none") {
            ;
        } else {
            BOOST_LOG_TRIVIAL(fatal)
                << "[ERROR]: param init must be one of: "
                "'random', 'forgy', 'kmeanspp'.It is '"
                << init << "'";
            exit(EXIT_FAILURE);
        }

        BOOST_LOG_TRIVIAL(info) << "Init is '" << init << "'";
		BOOST_LOG_TRIVIAL(info) << "Matrix K-means starting ...";
#if 0
        std::string fn = "centers"+std::to_string(K)+"-"+std::to_string(NUM_COLS)+".bin";
        FILE* f = fopen(fn.c_str(), "wb");
        fwrite(clusters, sizeof(double)*K*NUM_COLS, 1, f);
        fclose(f);
#endif
		bool converged = false;

		std::string str_iters =  MAX_ITERS == std::numeric_limits<unsigned>::max() ?
			"until convergence ...":
			std::to_string(MAX_ITERS) + " iterations ...";
		BOOST_LOG_TRIVIAL(info) << "Computing " << str_iters;
		unsigned iter = 1;

		while (iter < MAX_ITERS) {
			// Hold cluster assignment counter

			BOOST_LOG_TRIVIAL(info) << "E-step Iteration " << iter <<
				" . Computing cluster assignments ...";
			E_step(matrix, clusters, cluster_assignments, cluster_assignment_counts);

			if (g_num_changed == 0 || ((g_num_changed/(double)NUM_ROWS)) <= tolerance) {
				converged = true;
				break;
			} else {
				g_num_changed = 0;
			}

#if KM_TEST
			BOOST_LOG_TRIVIAL(info) << "M-step Updating cluster means ...";
#endif
			M_step(matrix, clusters, cluster_assignment_counts, cluster_assignments);
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
