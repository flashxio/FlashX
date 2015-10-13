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
#define KM_TEST 0
#define KM_ENABLE_DIST_TYPE 0

#ifdef PROFILER
#include <gperftools/profiler.h>
#endif

namespace {
	static unsigned NEV;
	static size_t K;
	static unsigned NUM_ROWS;
	short OMP_MAX_THREADS;
	static unsigned g_num_changed = 0;
    static const unsigned INVALID_CLUSTER_ID = std::numeric_limits<unsigned>::max();

	static struct timeval start, end;

	/**
	 * \brief Print an arry of some length `len`.
	 *	 \param len The length of the array.
	 */
	template <typename T>
		static void print_arr(T* arr, unsigned len) {
			printf("[ ");
			for (unsigned i = 0; i < len; i++) {
				std::cout << arr[i] << " ";
			}
			printf("]\n");
		}

	/*
	 * \Internal
	 * \brief Simple helper used to print a vector.
	 * \param v The vector to print.
	 */
	template <typename T>
		static void print_vector(typename std::vector<T>& v) {
			std::cout << "[";

			typename std::vector<T>::iterator itr = v.begin();
			for (; itr != v.end(); itr++) {
				std::cout << " "<< *itr;
			}
			std::cout <<  " ]\n";
		}

	/*\Internal
	 * \brief print a col wise matrix of type double / double.
	 * Used for testing only.
	 * \param matrix The col wise matrix.
	 * \param rows The number of rows in the mat
	 * \param cols The number of cols in the mat
	 */
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

    // TODO: Doc
    static double const eucl_dist(const unsigned row, const unsigned clust_idx,
                                const double* matrix, const double* clusters) {
        double dist = 0;
        for (unsigned col = 0; col < NEV; col++) {
            double diff = matrix[(row*NEV) + col] - clusters[clust_idx*NEV + col];
            dist += diff * diff;
        }
        return  dist;
    }
    

    /**
      * TODO: Doc
      */
    static double const cos_dist(const unsigned row, const unsigned clust_idx,
                                const double* matrix, const double* clusters) {
        double numr, ldenom, rdenom;
        numr = ldenom = rdenom = 0;

        for (unsigned col = 0; col < NEV; col++) {
            double a = matrix[(row*NEV) + col];
            double b = clusters[clust_idx*NEV + col];

            numr += a*b;
            ldenom += a*a;
            rdenom += b*b;
        }
        return  1 - (numr / ((sqrt(ldenom)*sqrt(rdenom))));
    }

    // TODO: Doc
#if KM_ENABLE_DIST_TYPE
    static size_t get_best_cluster(const unsigned row, const double* matrix, const double* clusters, 
            double (const *dist_func) (const unsigned, const unsigned, const double*, const double*)) {
#else
    static size_t get_best_cluster(const unsigned row, const double* matrix, const double* clusters) {
#endif

        size_t asgnd_clust = INVALID_CLUSTER_ID;
        double best = std::numeric_limits<double>::max();
        for (unsigned clust_idx = 0; clust_idx < K; clust_idx++) {
#if KM_ENABLE_DIST_TYPE
            double dist = dist_func(row, clust_idx, matrix, clusters);
#else
            double dist = eucl_dist(row, clust_idx, matrix, clusters);
#endif
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
	 *	\param cluster_assignments Which cluster each sample falls into.
	 */
	static void random_partition_init(unsigned* cluster_assignments) {
		BOOST_LOG_TRIVIAL(info) << "Random init start";

		// #pragma omp parallel for firstprivate(cluster_assignments, K) shared(cluster_assignments)
		for (unsigned vid = 0; vid < NUM_ROWS; vid++) {
			size_t assigned_cluster = random() % K; // 0...K
			cluster_assignments[vid] = assigned_cluster;
		}

		// NOTE: M-Step is called in compute func to update cluster counts & centers
#if KM_TEST
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
			memcpy(&clusters[clust_idx*NEV], &matrix[rand_idx*NEV], sizeof(clusters[0])*NEV);
		}

		BOOST_LOG_TRIVIAL(info) << "Forgy init end";
	}

	/**
	 * \brief A parallel version of the kmeans++ initialization alg.
	 *  See: http://ilpubs.stanford.edu:8090/778/1/2006-13.pdf for algorithm
	 */
#if KM_ENABLE_DIST_TYPE
	static void kmeanspp_init(const double* matrix, double* clusters, 
            double (const *dist_func)(const unsigned, const unsigned, const double*, const double*)) {
#else
	static void kmeanspp_init(const double* matrix, double* clusters) {
#endif

		// Init distance array
		double* dist_v = new double [NUM_ROWS*sizeof(double)];
		for (size_t i = 0; i < NUM_ROWS; i++)
			dist_v[i] = std::numeric_limits<double>::max();

		// Choose c1 uniiformly at random
		unsigned rand_idx = random() % (NUM_ROWS - 1); // 0...(n-1)

		memcpy(&clusters[0], &matrix[rand_idx*NEV], sizeof(clusters[0])*NEV);
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

#if KM_ENABLE_DIST_TYPE
				double dist = dist_func(row, clust_idx, matrix, clusters); // Per thread instance
#else
				double dist = eucl_dist(row, clust_idx, matrix, clusters); // Per thread instance
#endif

				if (dist < dist_v[row]) { // Found a closer cluster than before
					dist_v[row] = dist;
				}
				cum_dist += dist_v[row];
			}

			cum_dist = (cum_dist * ((double)random())) / (RAND_MAX - 1.);
			clust_idx++;

			for (size_t i=0; i < NUM_ROWS; i++) {
				cum_dist -= dist_v[i];
				if (cum_dist <= 0) {
#if KM_TEST
					BOOST_LOG_TRIVIAL(info) << "Choosing " << i << " as center K = " << clust_idx;
#endif
					memcpy(&(clusters[clust_idx*NEV]), &(matrix[i*NEV]), sizeof(clusters[0])*NEV);
					break;
				}
			}
			assert (cum_dist <= 0);
		}
		delete [] dist_v;

#if KM_TEST
		BOOST_LOG_TRIVIAL(info) << "\n\nOriginal matrix: "; print_mat(matrix, NUM_ROWS, NEV);
		BOOST_LOG_TRIVIAL(info) << "\nCluster centers after kmeans++"; print_mat(clusters, K, NEV);
#endif
	}

	/**
	 * \brief Update the cluster assignments while recomputing distance matrix.
	 * \param matrix The flattened matrix who's rows are being clustered.
	 * \param clusters The cluster centers (means) flattened matrix.
	 *	\param cluster_assignments Which cluster each sample falls into.
	 */
#if KM_ENABLE_DIST_TYPE
	static void E_step(const double* matrix, double* clusters,
			unsigned* cluster_assignments, unsigned* cluster_assignment_counts,
            const double (const *dist_func)(const unsigned, const unsigned, const double*, const double*)) {
#else
	static void E_step(const double* matrix, double* clusters,
			unsigned* cluster_assignments, unsigned* cluster_assignment_counts) {
#endif


#if KM_TEST
		gettimeofday(&start, NULL);
#endif

		// Create per thread vectors
		std::vector<std::vector<unsigned>> pt_cl_as_cnt; // K * OMP_MAX_THREADS
		std::vector<std::vector<double>> pt_cl; // K * nev * OMP_MAX_THREADS
		std::vector<size_t> pt_num_change; // Per thread changed cluster count. OMP_MAX_THREADS

		pt_cl_as_cnt.resize(OMP_MAX_THREADS);
		pt_cl.resize(OMP_MAX_THREADS);
		pt_num_change.resize(OMP_MAX_THREADS);

		for (int i = 0; i < OMP_MAX_THREADS; i++) {
			pt_cl_as_cnt[i].resize(K); // C++ default is 0 value in all
			pt_cl[i].resize(K*NEV);
		}

#pragma omp parallel for firstprivate(matrix, clusters) shared(pt_cl_as_cnt, cluster_assignments)\
		schedule(static)
		for (unsigned row = 0; row < NUM_ROWS; row++) {
#if KM_ENABLE_DIST_TYPE
			size_t asgnd_clust = get_best_cluster(row, matrix, clusters, dist_func);
#else
			size_t asgnd_clust = get_best_cluster(row, matrix, clusters);
#endif

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
			for (unsigned col = 0; col < NEV; col++) {
				pt_cl[omp_get_thread_num()][asgnd_clust*NEV + col] =
					pt_cl[omp_get_thread_num()][asgnd_clust*NEV + col] + matrix[row*NEV + col];
			}
		}

#if KM_TEST
		BOOST_LOG_TRIVIAL(info) << "Clearing cluster assignment counts";
		BOOST_LOG_TRIVIAL(info) << "Clearing cluster centers ...";
#endif
		memset(cluster_assignment_counts, 0, sizeof(cluster_assignment_counts[0])*K);
		memset(clusters, 0, sizeof(clusters[0])*K*NEV);

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
				for (unsigned col = 0; col < NEV; col++) { /* Features */
					clusters[row*NEV+col] = pt_cl[i][row*NEV + col] + clusters[row*NEV + col];
				}
			}
		}

#if KM_TEST
		gettimeofday(&end, NULL);
		BOOST_LOG_TRIVIAL(info) << "E-step time taken = " << time_diff(start, end) << " sec";
#endif

#if KM_TEST
		printf("Cluster assignment counts: "); print_arr(cluster_assignment_counts, K);
		printf("Cluster assignments:\n"); print_arr(cluster_assignments, NUM_ROWS);
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
			memset(clusters, 0, sizeof(clusters[0])*K*NEV);

			BOOST_LOG_TRIVIAL(info) << "Clearing cluster assignment counts";
			memset(cluster_assignment_counts, 0, sizeof(cluster_assignment_counts[0])*K);

			for (unsigned vid = 0; vid < NUM_ROWS; vid++) {

				unsigned asgnd_clust = cluster_assignments[vid];
				cluster_assignment_counts[asgnd_clust]++;

				for (unsigned col = 0; col < NEV; col++) {
					clusters[asgnd_clust*NEV + col] =
						clusters[asgnd_clust*NEV + col] + matrix[vid*NEV + col];
				}
			}
		}

		// Take the mean of all added means
#if KM_TEST
		BOOST_LOG_TRIVIAL(info) << " Div by in place ...";
		double* clusters = __builtin_assume_aligned(clusters);
#endif
		for (unsigned clust_idx = 0; clust_idx < K; clust_idx++) {
			if (cluster_assignment_counts[clust_idx] > 0) { // Avoid div by 0
				// #pragma omp parallel for firstprivate(cluster_assignment_counts, NEV) shared(clusters)
				for (unsigned col = 0; col < NEV; col++) {
					clusters[clust_idx*NEV + col] =
						clusters[clust_idx*NEV + col] / cluster_assignment_counts[clust_idx];
				}
			}
		}

# if KM_TEST
		printf("Cluster centers: \n"); print_mat(clusters, K, NEV);
#endif
	}
}

namespace fg
{
	unsigned compute_kmeans(const double* matrix, double* clusters,
			unsigned* cluster_assignments, unsigned* cluster_assignment_counts,
			const unsigned num_rows, const unsigned nev, const size_t k, const unsigned MAX_ITERS,
			const int max_threads, const std::string init, const double tolerance,
            const std::string dist_type)
	{
#ifdef PROFILER
		ProfilerStart("/home/disa/FlashGraph/flash-graph/matrix/kmeans.perf");
#endif
		NEV = nev;
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
		memset(cluster_assignments, -1,
				sizeof(cluster_assignments[0])*NUM_ROWS);
		memset(cluster_assignment_counts, 0, sizeof(cluster_assignment_counts[0])*K);

		/*** End VarInit ***/
#if KM_ENABLE_DIST_TYPE
        const double (const *dist_func)(const unsigned, const unsigned, const double*, const double*);
        if (dist_type == "eucl") {
            BOOST_LOG_TRIVIAL(info) << "Dist_type is " << dist_type;
            dist_func = &eucl_dist;
        } else if (dist_type == "cos") {
            BOOST_LOG_TRIVIAL(info) << "Dist_type is " << dist_type;
            dist_func = &cos_dist;
        } else {
            BOOST_LOG_TRIVIAL(fatal)
                << "[ERROR]: param dist_type must be one of: 'eucl', 'cos'.It is '"
                << dist_type << "'";
            exit(-1);
        }
#endif

        if (init == "random") {
            BOOST_LOG_TRIVIAL(info) << "Init is '" << init << "'";
            random_partition_init(cluster_assignments);
            // We must now update cluster centers before we begin
            M_step(matrix, clusters, cluster_assignment_counts, cluster_assignments, true);
        } else if (init == "forgy") {
            BOOST_LOG_TRIVIAL(info) << "Init is '" << init << "'";
            forgy_init(matrix, clusters);
        } else if (init == "kmeanspp") {
            BOOST_LOG_TRIVIAL(info) << "Init is '" << init << "'";
#if KM_ENABLE_DIST_TYPE
            kmeanspp_init(matrix, clusters, dist_func);
#else
            kmeanspp_init(matrix, clusters);
#endif
        } else {
            BOOST_LOG_TRIVIAL(fatal)
                << "[ERROR]: param init must be one of: 'random', 'forgy', 'kmeanspp'.It is '"
                << init << "'";
            exit(-1);
        }

		BOOST_LOG_TRIVIAL(info) << "Matrix K-means starting ...";

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
#if KM_ENABLE_DIST_TYPE
			E_step(matrix, clusters, cluster_assignments, cluster_assignment_counts, dist_func);
#else
			E_step(matrix, clusters, cluster_assignments, cluster_assignment_counts);
#endif

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
