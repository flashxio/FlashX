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

bool updated = false;
static vsize_t NEV;
static size_t K;
static vsize_t NUM_ROWS;

static struct timeval start, end;

template <typename T>
void print_arr(T* arr, vsize_t len) {
    printf("[ ");
    for (vsize_t i = 0; i < len; i++) {
        std::cout << arr[i] << " ";
    }
    printf("]\n");
}

template <typename T>
void print_vector(typename std::vector<T>& v) 
{
	std::cout << "[";

	typename std::vector<T>::iterator itr = v.begin();
	for (; itr != v.end(); itr++) {
		std::cout << " "<< *itr;
	}
	std::cout <<  " ]\n";
}

ev_float_t dist_sq(ev_float_t arg1, ev_float_t arg2) 
{
	ev_float_t diff = arg1 - arg2;
	return diff*diff;
}

void random_partition_init(vsize_t* cluster_assignments) 
{
	BOOST_LOG_TRIVIAL(info) << "Random init start";

// #pragma omp parallel for firstprivate(cluster_assignments, K) shared(cluster_assignments)
	for (vsize_t vid = 0; vid < NUM_ROWS; vid++) {
		size_t assigned_cluster = random() % K; // 0...K
		cluster_assignments[vid] = assigned_cluster;
	}

	// NOTE: M-Step is called in compute func to update cluster counts & centers 
#if 0
	printf("After rand paritions cluster_asgns: "); print_arr(cluster_assignments, NUM_ROWS);
#endif
	BOOST_LOG_TRIVIAL(info) << "Random init end\n";
}

void forgy_init(const ev_float_t* matrix, ev_float_t* clusters) {

	BOOST_LOG_TRIVIAL(info) << "Forgy init start";

// #pragma omp parallel for firstprivate(matrix) shared (clusters)
	for (vsize_t clust_idx = 0; clust_idx < K; clust_idx++) { // 0...K
		vsize_t rand_idx = random() % (NUM_ROWS - 1); // 0...(n-1)
		memcpy(&clusters[clust_idx*NEV], &matrix[rand_idx*NEV], sizeof(clusters[0])*NEV);
	}
	BOOST_LOG_TRIVIAL(info) << "Forgy init end";
}

void kmeanspp_init()
{
	BOOST_LOG_TRIVIAL(fatal) << "kmeanspp not yet implemented";
	exit(-1);
}

template <typename T>
void print_mat(T* matrix, const vsize_t rows, const vsize_t cols) {
	for (vsize_t row = 0; row < rows; row++) {
		std::cout << "[";
		for (vsize_t col = 0; col < cols; col++) {
			std::cout << " " << matrix[row*cols + col];
		}
		std::cout <<  " ]\n";
	}	
}

void E_step(const ev_float_t* matrix, ev_float_t* clusters, vsize_t* cluster_assignments)
{
#pragma omp parallel for firstprivate(matrix, clusters) shared(cluster_assignments)
    for (vsize_t row=0; row < NUM_ROWS; row++) {
        size_t asgnd_clust = std::numeric_limits<size_t>::max();
        ev_float_t best = std::numeric_limits<ev_float_t>::max();
        for (size_t clust_idx = 0; clust_idx < K; clust_idx++) {
            ev_float_t sum = 0;
            for (vsize_t col = 0; col < NEV; col++) {
                sum += dist_sq(matrix[(row*NEV) + col], clusters[clust_idx*NEV + col]);
            }
            if (sum < best) {
                best = sum;
                asgnd_clust = clust_idx;
            }
        }

		if (asgnd_clust != cluster_assignments[row]) { 
			if (!updated) updated = true;
		}
        cluster_assignments[row] = asgnd_clust;
    }

#if 0
    printf("Cluster assignments:\n"); print_arr(cluster_assignments, NUM_ROWS);
#endif
}

void M_step(const ev_float_t* matrix, ev_float_t* clusters, 
		vsize_t* cluster_assignment_counts, const vsize_t* cluster_assignments)
{
	BOOST_LOG_TRIVIAL(info) << "M_step start";

    gettimeofday(&start, NULL);
	BOOST_LOG_TRIVIAL(info) << "Clearing cluster centers ...";
	memset(clusters, 0, sizeof(clusters[0])*K*NEV);

	BOOST_LOG_TRIVIAL(info) << "Clearing cluster assignments";
	memset(cluster_assignment_counts, 0, sizeof(cluster_assignment_counts[0])*K);

	// TODO: Parallelize me somehow
	for (vsize_t vid = 0; vid < NUM_ROWS; vid++) {

		// Update cluster counts
		vsize_t asgnd_clust = cluster_assignments[vid];
		cluster_assignment_counts[asgnd_clust]++;	

		for (vsize_t col = 0; col < NEV; col++) {
			clusters[asgnd_clust*NEV + col] = clusters[asgnd_clust*NEV + col] + matrix[vid*NEV + col];	
		}
	}

#if 0
	printf("Cluster assignment counts: "); print_arr(cluster_assignment_counts, K);
	printf("Cluster centers: \n"); print_mat(clusters, K, NEV);
#endif

	// Take the mean of all added means
	BOOST_LOG_TRIVIAL(info) << " Div by in place ...";
//#pragma omp parallel for firstprivate(cluster_assignment_counts, K, NEV) shared(clusters)
	for (vsize_t clust_idx = 0; clust_idx < K; clust_idx++) {
		if (cluster_assignment_counts[clust_idx] > 0) { // Avoid div by 0
#pragma omp parallel for firstprivate(cluster_assignment_counts, NEV) shared(clusters)
			for (vsize_t col = 0; col < NEV; col++) {
				clusters[clust_idx*NEV + col] = 
					clusters[clust_idx*NEV + col] / cluster_assignment_counts[clust_idx];
			}
		}
	}
    gettimeofday(&end, NULL);
    BOOST_LOG_TRIVIAL(info) << "M-step time taken = %.3f\n", time_diff(start, end);

# if 0
	printf("Cluster centers: \n"); print_mat(clusters, K, NEV);
#endif
}

vsize_t compute_kmeans(const ev_float_t* matrix, ev_float_t* clusters, 
		vsize_t* cluster_assignments, vsize_t* cluster_assignment_counts,
		const vsize_t num_rows, const vsize_t nev, const size_t k, const vsize_t MAX_ITERS, 
		const std::string init)
{
	NEV = nev;
	K = k;
	NUM_ROWS = num_rows;

	if (K > NUM_ROWS || K < 2) { 
		BOOST_LOG_TRIVIAL(fatal)
			<< "'k' must be between 2 and the number of rows in the matrix";
		exit(-1);
	}

	if (init.compare("random") != 0 && init.compare("kmeanspp") != 0 &&
			init.compare("forgy") != 0) {
		BOOST_LOG_TRIVIAL(fatal)
			<< "[ERROR]: param init must be one of: 'random', 'kmeanspp'.It is '"
			<< init << "'";
		exit(-1);
	}

#ifdef PROFILER
	ProfilerStart(PROF_FILE_LOC);
#endif

	/*** Begin VarInit of data structures ***/
	memset(cluster_assignments, INVALID_VERTEX_ID, 
			sizeof(cluster_assignments[0])*NUM_ROWS);
	memset(cluster_assignment_counts, 0, sizeof(cluster_assignment_counts[0])*K);

	/*** End VarInit ***/

	if (init == "random") {
		BOOST_LOG_TRIVIAL(info) << "Init is random_partition \n";
		random_partition_init(cluster_assignments);
		// We must now update cluster centers before we begin
		M_step(matrix, clusters, cluster_assignment_counts, cluster_assignments);
	}
	if (init == "forgy") {
		BOOST_LOG_TRIVIAL(info) << "Init is forgy";
		forgy_init(matrix, clusters);
	}
	else {
		if (init == "kmeanspp") {
			kmeanspp_init(); // TODO: kmeanspp
		}
	}
	BOOST_LOG_TRIVIAL(info) << "Matrix K-means starting ...";

	bool converged = false;

	std::string str_iters =  MAX_ITERS == std::numeric_limits<vsize_t>::max() ? "until convergence ...":
		std::to_string(MAX_ITERS) + " iterations ...";
	BOOST_LOG_TRIVIAL(info) << "Computing " << str_iters; 
	vsize_t iter = 1;
	while (iter < MAX_ITERS) {
		// Hold cluster assignment counter

		BOOST_LOG_TRIVIAL(info) << "E-step Iteration " << iter << " . Computing cluster assignments ...";
		E_step(matrix, clusters, cluster_assignments);

		if (!updated) {
			converged = true;
			break;
		} else { updated = false; }

		BOOST_LOG_TRIVIAL(info) << "M-step Updating cluster means ...";

		M_step(matrix, clusters, cluster_assignment_counts, cluster_assignments);
		iter++;
	}

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

	BOOST_LOG_TRIVIAL(info) << "\nCluster counts: ";
	return iter;
}
