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
#include <numeric>

#include "sem_kmeans.h"
#include "../../../../libkcommon/clusters.hpp"

using namespace fg;

namespace {

    static size_t g_io_reqs = 0;

    static kbase::clusters::ptr g_clusters; // cluster means/centers

    static unsigned NUM_ROWS;
    static unsigned g_num_changed = 0;
    static struct timeval start, end;
    static kbase::init_t g_init; // May have to use
    static unsigned  g_kmspp_cluster_idx; // Used for kmeans++ init
    static unsigned g_kmspp_next_cluster; // Sample row selected as next cluster
    static kmspp_stage_t g_kmspp_stage; // Either adding a mean / computing dist
    static kbase::stage_t g_stage; // What phase of the algo we're in
    static unsigned g_iter;
    static std::vector<vertex_id_t> all_vertices;
    static barrier::ptr iter_barrier;
    static unsigned g_max_iters;
    static std::vector<size_t> g_num_members_v;
    static double g_tolerance;
    static bool g_converged = false;
    static std::default_random_engine generator;
    static std::uniform_real_distribution<double> ur_distribution(0.0, 1.0);

    static unsigned kmeanspp_get_next_cluster_id(graph_engine::ptr mat);
    void update_clusters(graph_engine::ptr mat,
            std::vector<size_t>& g_num_members_v);
    class kmeans_vertex: public base_kmeans_vertex
    {
        public:
        kmeans_vertex(vertex_id_t id): base_kmeans_vertex(id) { }

        void run(vertex_program &prog) {
            vertex_id_t id = prog.get_vertex_id(*this);
            switch (g_stage) {
                case kbase::stage_t::INIT:
                    request_vertices(&id, 1);
                    break;
                case kbase::stage_t::ESTEP:
                    if (!g_converged && (g_iter < g_max_iters)) {
                        prog.activate_vertices(&id, 1); // Activate for next iter
                        request_vertices(&id, 1);
                    } else {
                        return;
                    }
                    break;
                default:
                    BOOST_ASSERT_MSG(0, "Unknown g_stage!");
            }
        }

        void run(vertex_program& prog, const page_vertex &vertex) {
            switch (g_stage) {
                case kbase::stage_t::INIT:
                    run_init(prog, vertex, g_init);
                    break;
                case kbase::stage_t::ESTEP:
                    run_distance(prog, vertex);
                    break;
                default:
                    BOOST_ASSERT_MSG(0, "Unknown g_stage!");
            }
        }

        void run_on_message(vertex_program& prog, const vertex_message& msg) { }
        void run_init(vertex_program& prog,
                const page_vertex &vertex, kbase::init_t init);
        void run_distance(vertex_program& prog, const page_vertex &vertex);
        double dist_comp(const page_vertex &vertex, const double* mean);
    };

    class kmeans_vertex_program:
        public base_kmeans_vertex_program<kmeans_vertex, kbase::clusters> {

        private:
            graph_engine::ptr mat;

        public:
        typedef std::shared_ptr<kmeans_vertex_program> ptr;

        kmeans_vertex_program(graph_engine::ptr mat) {
            this->mat = mat;
        }

        static ptr cast2(vertex_program::ptr prog) {
            return std::static_pointer_cast<kmeans_vertex_program, vertex_program>(prog);
        }

        void remove_member(const unsigned id, data_seq_iter& count_it) {
            get_pt_clusters()->remove_member(count_it, id);
        }

        void swap_membership(data_seq_iter& count_it, const unsigned from_id,
                const unsigned to_id) {
            get_pt_clusters()->swap_membership(count_it, from_id, to_id);
        }

        // Per partition
        virtual void run_on_iteration_end() override {
            if (iter_barrier->ping()) { // Make sure all partitions are complete 1st
                if (!g_converged) {
                    BOOST_LOG_TRIVIAL(info) << "E-step Iteration " << g_iter <<
                        " . Computing cluster assignments ...";
                    g_io_reqs += NUM_ROWS;

                    BOOST_LOG_TRIVIAL(info) << "Main: M-step Updating "
                        << "cluster means ...";
                    update_clusters(mat, g_num_members_v);

                    BOOST_LOG_TRIVIAL(info) << "Printing cluster counts ...";
                    kbase::print_vector(g_num_members_v);

                    BOOST_LOG_TRIVIAL(info) << "** Samples changes cluster: "
                        << g_num_changed << " **\n";

                    if ((g_num_changed == 0 ||
                            ((g_num_changed/(double)NUM_ROWS)) <= g_tolerance)
                            && g_iter < g_max_iters) {
                        g_converged = true;
                    } else {
                        g_num_changed = 0;
                        g_iter++;
                    }
                }
            }
        }
    };

    class kmeans_vertex_program_creater: public vertex_program_creater
    {
        graph_engine::ptr mat;

        public:
        kmeans_vertex_program_creater (graph_engine::ptr mat) {
            this->mat = mat;
        }

        vertex_program::ptr create() const {
            return vertex_program::ptr(new kmeans_vertex_program(mat));
        }
    };

    /* Used in kmeans++ initialization */
    class kmeanspp_vertex_program : public vertex_program_impl<kmeans_vertex>
    {
        double pt_cuml_sum;
        graph_engine::ptr mat;

        public:
        typedef std::shared_ptr<kmeanspp_vertex_program> ptr;

        kmeanspp_vertex_program(graph_engine::ptr mat) {
            pt_cuml_sum = 0.0;
            this->mat = mat;
        }

        static ptr cast2(vertex_program::ptr prog) {
            return std::static_pointer_cast<kmeanspp_vertex_program,
                   vertex_program>(prog);
        }

        void pt_cuml_sum_peq (const double val) {
            pt_cuml_sum += val;
        }

        const double get_pt_cuml_sum() const {
            return pt_cuml_sum;
        }

        void reset() { pt_cuml_sum = 0.0; }

        // Per partition
        virtual void run_on_iteration_end() override {
            if (iter_barrier->ping()) { // Make sure all partitions are complete 1st
                if (g_kmspp_stage == DIST) {
                    // Last thread calls to pick a new cluster index
#if VERBOSE
                    BOOST_LOG_TRIVIAL(info) << "Printing clusters "
                        << "after sample set_mean ...";
                    g_clusters->print_means();
#endif
                    if (g_kmspp_cluster_idx+1 < K) { // Still have work to do?
                        g_io_reqs += (NUM_ROWS + 1);
                        g_kmspp_next_cluster = kmeanspp_get_next_cluster_id(mat);
                        // Activate the new vertex
                        activate_vertices(&g_kmspp_next_cluster, 1);
                    }
                    g_kmspp_stage = ADDMEAN;
                } else {
                    g_kmspp_stage = DIST;
                }
            }
        }
    };

    class kmeanspp_vertex_program_creater: public vertex_program_creater
    {
        graph_engine::ptr mat;

        public:
        kmeanspp_vertex_program_creater (graph_engine::ptr mat) {
            this->mat = mat;
        }

        vertex_program::ptr create() const {
            return vertex_program::ptr(new kmeanspp_vertex_program(mat));
        }
    };

    double kmeans_vertex::dist_comp(const page_vertex &vertex, const double* mean) {
        data_seq_iter count_it =
            ((const page_row&)vertex).get_data_seq_it<double>();

        double dist = 0;
        double diff;
        vertex_id_t nid = 0;

        while(count_it.has_next()) {
            double e = count_it.next();
            diff = e - mean[nid++];
            dist += diff*diff;
        }
        BOOST_VERIFY(nid == NUM_COLS);
        return sqrt(dist);
    }

    void kmeans_vertex::run_init(vertex_program& prog,
            const page_vertex &vertex, kbase::init_t init) {
        switch (g_init) {
            case kbase::init_t::RANDOM:
                {
                    kmeans_vertex_program& vprog = (kmeans_vertex_program&) prog;
                    unsigned new_cluster_id = random() % K;
#if VERBOSE
                    printf("Random init: v%u assigned to cluster: c%x\n",
                            prog.get_vertex_id(*this), new_cluster_id);
#endif
                    set_cluster_id(new_cluster_id);
                    data_seq_iter count_it = ((const page_row&)vertex).
                        get_data_seq_it<double>();
                    vprog.add_member(get_cluster_id(), count_it);
                }
                break;
            case kbase::init_t::FORGY:
                {
                    vertex_id_t my_id = prog.get_vertex_id(*this);
#if KM_TEST
                    printf("Forgy init: v%u setting cluster: c%x\n", my_id, g_init_hash[my_id]);
#endif
                    data_seq_iter count_it = ((const page_row&)vertex).
                        get_data_seq_it<double>();
                    g_clusters->set_mean(count_it, g_init_hash[my_id]);
                }
                break;
            case kbase::init_t::PLUSPLUS:
                {
                    data_seq_iter count_it = ((const page_row&)vertex).
                        get_data_seq_it<double>();

                    vertex_id_t my_id = prog.get_vertex_id(*this);
                    if (g_kmspp_stage == ADDMEAN) {
#if KM_TEST
                        printf("kms++ v%u making itself c%u\n", my_id, g_kmspp_cluster_idx);
#endif
                        set_cluster_id(g_kmspp_cluster_idx);
                        g_kmspp_distance[my_id] = 0;
                        g_clusters->add_member(count_it, g_kmspp_cluster_idx);
                        // Activate all
                        prog.activate_vertices(&all_vertices[0], NUM_ROWS);
                    } else if (g_kmspp_stage == DIST) {
                        double _dist = dist_comp(vertex,
                                &(g_clusters->get_means()[g_kmspp_cluster_idx*NUM_COLS]));

                        if (_dist < g_kmspp_distance[my_id]) {
                            g_kmspp_distance[my_id] = _dist;
                            set_cluster_id(g_kmspp_cluster_idx);
                        }
                        ((kmeanspp_vertex_program&)prog).
                            pt_cuml_sum_peq(g_kmspp_distance[my_id]);
                    } else {
                        BOOST_ASSERT_MSG(0, "Unknown g_kmspp_stage type");
                    }
                }
                break;
            default:
                BOOST_ASSERT_MSG(0, "Unknown g_init type");
        }
    }

    void kmeans_vertex::run_distance(vertex_program& prog, const page_vertex &vertex) {
        kmeans_vertex_program& vprog = (kmeans_vertex_program&) prog;
        unsigned old_cluster_id = get_cluster_id();

        double best = std::numeric_limits<double>::max();
        for (unsigned cl = 0; cl < K; cl++) {
            double udist = dist_comp(vertex, &(g_clusters->get_means()[cl*NUM_COLS]));
            if (udist < best) {
                best = udist;
                set_cluster_id(cl);
            }
        }

#if KM_TEST
        BOOST_VERIFY(get_cluster_id() >= 0 && get_cluster_id() < K);
#endif
        data_seq_iter count_it = ((const page_row&)vertex).get_data_seq_it<double>();

        vprog.add_member(get_cluster_id(), count_it);
        if (old_cluster_id != get_cluster_id()) {
            vprog.pt_changed_pp(); // Add a vertex to the count of changed ones
        }
    }

        static FG_vector<unsigned>::ptr get_membership(graph_engine::ptr mat) {
            FG_vector<unsigned>::ptr vec = FG_vector<unsigned>::create(mat);
            mat->query_on_all(vertex_query::ptr(new save_query<unsigned, kmeans_vertex>(vec)));
            return vec;
        }

        void update_clusters(graph_engine::ptr mat,
                std::vector<size_t>& g_num_members_v) {
            g_clusters->clear();

            std::vector<vertex_program::ptr> kms_clust_progs;
            mat->get_vertex_programs(kms_clust_progs);

            for (unsigned thd = 0; thd < kms_clust_progs.size(); thd++) {
                kmeans_vertex_program::ptr kms_prog =
                    kmeans_vertex_program::cast2(kms_clust_progs[thd]);
                kbase::clusters::ptr pt_clusters = kms_prog->get_pt_clusters();
                g_num_changed += kms_prog->get_pt_changed();

                BOOST_VERIFY(g_num_changed <= NUM_ROWS);
                /* Merge the per-thread clusters */
                g_clusters->peq(pt_clusters);

                kms_prog->reset();
            }

            for (unsigned cl = 0; cl < K; cl++) {
                g_clusters->finalize(cl);
                g_num_members_v[cl] = g_clusters->get_num_members(cl);
            }
#if KM_TEST
            long t_members = 0;
            for (unsigned cl = 0; cl < K; cl++) {
                t_members += g_clusters->get_num_members(cl);
                if (t_members > (long) NUM_ROWS) {
                    BOOST_LOG_TRIVIAL(error) << "[FATAL]: Too many members cluster: "
                        << cl << "/" << K << " at members = " << t_members;
                    BOOST_VERIFY(false);
                }
            }
#endif
        }

        /* During kmeans++ we select a new cluster each iteration
           This step get the next sample selected as a cluster center
           */
        static unsigned kmeanspp_get_next_cluster_id(graph_engine::ptr mat) {
#if KM_TEST
            BOOST_LOG_TRIVIAL(info) << "Assigning new cluster ...";
#endif
            std::vector<vertex_program::ptr> kmspp_progs;
            mat->get_vertex_programs(kmspp_progs);

            double cuml_sum = 0;
            BOOST_FOREACH(vertex_program::ptr vprog, kmspp_progs) {
                kmeanspp_vertex_program::ptr kmspp_prog =
                    kmeanspp_vertex_program::cast2(vprog);
                cuml_sum += kmspp_prog->get_pt_cuml_sum();

                // NOTE: We need to reset the cumulative sums
                kmspp_prog->reset();
            }

            cuml_sum = (cuml_sum * ur_distribution(generator))/(RAND_MAX - 1.0);
            BOOST_ASSERT_MSG(cuml_sum != 0, "Cumulative sum == 0!");

            g_kmspp_cluster_idx++;

            for (unsigned row = 0; row < NUM_ROWS; row++) {
#if VERBOSE
                BOOST_LOG_TRIVIAL(info) << "cuml_sum = " << cuml_sum;
#endif
                cuml_sum -= g_kmspp_distance[row];
                if (cuml_sum <= 0) {
#if KM_TEST
                    BOOST_LOG_TRIVIAL(info) << "Choosing v:" << row
                        << " as center K = " << g_kmspp_cluster_idx;
#endif
                    return row;
                }
            }
            BOOST_ASSERT_MSG(false, "Cumulative sum of distances was > than distances!");
            exit(EXIT_FAILURE);
        }

        static inline bool fexists(const std::string& name) {
            struct stat buffer;
            return (stat (name.c_str(), &buffer) == 0);
        }
    }

    namespace fg
    {
        void compute_sem_kmeans(FG_graph::ptr fg,
                const unsigned k, const std::string init,
                const unsigned max_iters, const double tolerance,
                kbase::kmeans_t& ret, const unsigned num_rows,
                const unsigned num_cols, std::vector<double>* centers) {
#ifdef PROFILER
            ProfilerStart("libgraph-algs/min_tri_sem_kmeans.perf");
#endif
            K = k;
            g_max_iters = max_iters;
            g_tolerance = tolerance;

            // Check Initialization
            if ((NULL == centers) && init.compare("random") && init.compare("kmeanspp") &&
                    init.compare("forgy")) {
                BOOST_LOG_TRIVIAL(fatal)
                    << "[ERROR]: init must be one of: 'random', 'forgy', 'kmeanspp'.It is '"
                    << init << "'";
                exit(EXIT_FAILURE);
            }

            graph_index::ptr index = NUMA_graph_index<kmeans_vertex>::create(
                    fg->get_graph_header());
            graph_engine::ptr mat = fg->create_engine(index);

            NUM_ROWS = num_rows;
            NUM_COLS = num_cols;

            // Iteration barrier configured to require nthreads to exit
            unsigned nthreads = atoi(fg->get_configs()->get_option("threads").c_str());
            iter_barrier = barrier::create(nthreads);

            // Check k
            if (K > NUM_ROWS || K < 2 || K == (unsigned)-1) {
                BOOST_LOG_TRIVIAL(fatal)
                    << "'k' must be between 2 and the number of rows in the matrix " <<
                    "k = " << K;
                exit(EXIT_FAILURE);
            }

            BOOST_VERIFY(num_cols > 0);

            BOOST_LOG_TRIVIAL(info) << "Matrix has rows = " << NUM_ROWS << ", cols = " <<
                NUM_COLS;
            gettimeofday(&start , NULL);

            /*** Begin VarInit of data structures ***/
            g_clusters = kbase::clusters::create(K, NUM_COLS);
            if (centers)
                g_clusters->set_mean(*centers); // (*centers) // &(*(centers)[0])

            g_num_members_v.assign(K, 0);

            /*** End VarInit ***/

            if (!centers) {
                g_stage = kbase::stage_t::INIT;

                if (init == "random") {
                    BOOST_LOG_TRIVIAL(info) << "Running init: '"<< init <<"' ...";
                    g_init = kbase::init_t::RANDOM;

                    mat->start_all(vertex_initializer::ptr(),
                            vertex_program_creater::ptr(
                                new kmeans_vertex_program_creater(mat)));
                    mat->wait4complete();
                    g_io_reqs += NUM_ROWS;

                    update_clusters(mat, g_num_members_v);
                }
                if (init == "forgy") {
                    BOOST_LOG_TRIVIAL(info) << "Deterministic Init is: '"<< init <<"'";
                    g_init = kbase::init_t::FORGY;
                    std::uniform_int_distribution<vertex_id_t>
                            distribution(0, NUM_ROWS-1);

                    // Select K in range NUM_ROWS
                    std::vector<vertex_id_t> init_ids; // Used to start engine
                    for (unsigned cl = 0; cl < K; cl++) {
                        vertex_id_t id = distribution(generator);
                        g_init_hash[id] = cl; // <vertex_id, cluster_id>
                        init_ids.push_back(id);
                    }
                    mat->start(&init_ids.front(), K);
                    mat->wait4complete();
                    g_io_reqs += 1;

                } else if (init == "kmeanspp") {
                    BOOST_LOG_TRIVIAL(info) << "Init is '"<< init <<"'";
                    // FIXME: Wasteful
                    all_vertices.resize(NUM_ROWS);
                    std::iota(all_vertices.begin(), all_vertices.end(), 0);
                    g_init = kbase::init_t::PLUSPLUS;
                    // Init g_kmspp_distance to max distance
                    g_kmspp_distance.assign(NUM_ROWS, std::numeric_limits<double>::max());

                    g_kmspp_cluster_idx = 0;

                    std::uniform_int_distribution<vertex_id_t>
                        distribution(0, NUM_ROWS-1);
                    g_kmspp_next_cluster = distribution(generator); // 0 - (NUM_ROWS - 1)

#if KM_TEST
                    BOOST_LOG_TRIVIAL(info) << "Assigning v:" << g_kmspp_next_cluster
                        << " as first cluster";
#endif
                    g_kmspp_distance[g_kmspp_next_cluster] = 0;

                    g_kmspp_stage = ADDMEAN;
                    mat->start(&g_kmspp_next_cluster, 1,
                            vertex_initializer::ptr(),
                            vertex_program_creater::ptr(
                                new kmeanspp_vertex_program_creater(mat)));
                    mat->wait4complete();
                }
            } else
                g_clusters->print_means();

            g_stage = kbase::stage_t::ESTEP;
            BOOST_LOG_TRIVIAL(info) << "knors No Pruning starting ...";
            std::string str_iters = g_max_iters == std::numeric_limits<unsigned>::max() ?
                "until convergence ...":
                std::to_string(g_max_iters) + " iterations ...";
            BOOST_LOG_TRIVIAL(info) << "Computing " << str_iters;
            g_iter = 0;

            if (max_iters > 0) {
                mat->start_all(vertex_initializer::ptr(),
                        vertex_program_creater::ptr(
                            new kmeans_vertex_program_creater(mat)));
                mat->wait4complete();
            }

            gettimeofday(&end, NULL);
            BOOST_LOG_TRIVIAL(info) << "\n\nAlgorithmic time taken = " <<
                time_diff(start, end) << " sec\n";

#ifdef PROFILER
            ProfilerStop();
#endif
            BOOST_LOG_TRIVIAL(info) << "\n******************************************\n";
            printf("Total # of IO requests: %lu\nTotal bytes requested: %lu\n\n",
                    g_io_reqs, (g_io_reqs*(sizeof(double))*NUM_COLS));

            if (g_converged) {
                BOOST_LOG_TRIVIAL(info) <<
                    "K-means converged in " << ++g_iter << " iterations";
            } else {
                BOOST_LOG_TRIVIAL(warning) << "[Warning]: K-means failed to converge in "
                    << g_max_iters << " iterations";
            }
            BOOST_LOG_TRIVIAL(info) << "\n******************************************\n";

            kbase::print_vector(g_num_members_v);

            std::vector<unsigned> mv(NUM_ROWS);
            get_membership(mat)->copy_to<unsigned>(&mv[0], NUM_ROWS);

            ret.set_params(NUM_ROWS, NUM_COLS, g_iter, K);
            ret.set_computed(&mv[0], &g_num_members_v[0],
                    g_clusters->get_means());
        }
    }