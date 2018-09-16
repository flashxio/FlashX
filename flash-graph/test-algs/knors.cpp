/*
 * Copyright 2017 Open Connectome Project (http://openconnecto.me)
 * Written by Disa Mhembere (disa@jhu.edu
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

#include <signal.h>
#ifdef PROFILER
#include <gperftools/profiler.h>
#endif

#include "FGlib.h"
#include "ts_graph.h"
#include "../libsafs/log.h"
#include "matrix/FG_sparse_matrix.h"
#include "libgraph-algs/sem_kmeans.h"

using namespace fg;

void print_usage();

void int_handler(int sig_num) {
#ifdef PROFILER
	if (!graph_conf.get_prof_file().empty())
		ProfilerStop();
#endif
	exit(0);
}

void run_sem_kmeans(FG_graph::ptr graph, const unsigned k,
        int argc, char *argv[]) {
	int opt;
	int num_opts = 0;

    unsigned max_iters = std::numeric_limits<unsigned>::max();
    std::string init = "kmeanspp";
    double tolerance = -1;
    std::vector<double>* centers = NULL;
    std::string init_centers_fn = "";
    std::string outdir = "";
    double cache_size_gb = 0;
    unsigned rc_update_start_interval = 5;
    bool no_prune = false;

	while ((opt = getopt(argc, argv, "i:t:l:C:r:I:Po:")) != -1) {
		num_opts++;
		switch (opt) {
			case 'i':
                max_iters = atol(optarg);
				break;
			case 't':
                init = optarg;
				break;
            case 'l':
                tolerance = atof(optarg);
                num_opts++;
                break;
            case 'C':
                init_centers_fn = optarg;
                num_opts++;
                break;
            case 'r':
                cache_size_gb = atof(optarg);
                num_opts++;
                break;
            case 'I':
                rc_update_start_interval = atoi(optarg);
                num_opts++;
                break;
            case 'P':
                no_prune = true;
                num_opts++;
                break;
            case 'o':
                outdir = std::string(optarg);
                num_opts++;
                break;
			default:
				print_usage();
				abort();
		}
	}

    if (!init_centers_fn.empty()) {
        BOOST_LOG_TRIVIAL(info) << "\nReading centers from disk at loc '" << init_centers_fn
            << "' ...";
        BOOST_VERIFY(centers = new std::vector<double>(k*graph->get_dim()));
        kpmbase::bin_io<double> br(init_centers_fn, k, graph->get_dim());
        br.read(centers);
    }

    kpmbase::kmeans_t ret;

    if (no_prune) {
         compute_sem_kmeans(graph, k, init, max_iters, tolerance,
                ret, graph->get_nsamples(), graph->get_dim(), centers);
    } else  {
#if 0 // No Full Elkans algorithm
        compute_triangle_sem_kmeans(graph, k, init,
                max_iters, tolerance, ret,
                graph->get_nsamples(), graph->get_dim(), centers);
#else
        compute_min_triangle_sem_kmeans(graph, k, init,
                max_iters, tolerance, ret, graph->get_nsamples(),
                graph->get_dim(), centers,
                cache_size_gb, rc_update_start_interval);
#endif
    }

    if (!outdir.empty()) {
        printf("\nWriting output to '%s'\n", outdir.c_str());
        ret.write(outdir);
    }

    if (centers) { delete centers; }
}

void print_usage() {
	fprintf(stderr,
			"knors config-file data-file nsamples dim k [alg-options]\n");
	fprintf(stderr, "-t: init type [random, forgy, kmeanspp]\n");
	fprintf(stderr, "-i: max number of iterations\n");
	fprintf(stderr, "-C: File with initial clusters in same format as data\n");
	fprintf(stderr, "-l: convergence tolerance (defualt: -1 = no changes)\n");
    fprintf(stderr, "-P DO NOT use the minimal triangle inequality (~Elkan's alg)\n");
	fprintf(stderr, "-r: size of the row cache in gb\n");
	fprintf(stderr, "-I: row cache update interval\n");
    fprintf(stderr, "-o Write output to an output directory of this name\n");
	fprintf(stderr, "\n");
	graph_conf.print_help();
	safs::params.print_help();
}

int main(int argc, char *argv[])
{
	if (argc < 6) {
		print_usage();
		exit(EXIT_FAILURE);
	}

	argv++;
	argc--;

	std::string conf_file = argv[0];
	std::string graph_file = argv[1];
	vsize_t nsamples = atol(argv[2]);
	vsize_t dim = atol(argv[3]);
	unsigned k = atol(argv[4]);

	// We should increase by 3 instead of 4. getopt() ignores the first
	// argument in the list.
	argv += 4;
	argc -= 4;

	config_map::ptr configs = config_map::create(conf_file);
	if (configs == NULL)
		configs = config_map::ptr();
	signal(SIGINT, int_handler);

    graph_engine::init_flash_graph(configs);

    {
        FG_graph::ptr graph =
            FG_graph::create(graph_file, "_", configs, nsamples, dim);
        run_sem_kmeans(graph, k, argc, argv);
    }

    graph_engine::destroy_flash_graph();
}