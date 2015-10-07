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

#include "sem_kmeans.h"
#include <signal.h>
#ifdef PROFILER
#include <gperftools/profiler.h>
#endif

#define PRUNE 1
#define KM_TEST 0
#define VERBOSE 0
#define INVALID_CLUST_ID -1

using namespace fg;

namespace {
    typedef std::pair<double, double> distpair;
    static unsigned NUM_COLS;
    static unsigned NUM_ROWS;
    static unsigned K;
    static unsigned g_num_changed = 0;
    static struct timeval start, end;
    static std::map<vertex_id_t, unsigned> g_init_hash; // Used for forgy init
    static unsigned  g_kmspp_cluster_idx; // Used for kmeans++ init
    static unsigned g_kmspp_next_cluster; // Sample row selected as the next cluster
    static std::vector<double> g_kmspp_distance; // Used for kmeans++ init
    static unsigned g_iter;
    static bool g_even_iter;
#if PRUNE
    static double g_COMP_THRESH;
    static std::vector<double> g_comp_thresh_v; // To prune computation
#endif

    enum dist_type_t { EUCL, COS }; // Euclidean, Cosine distance

    template <typename T>
        static void print_vector(typename std::vector<T> v) {
            std::cout << "[";
            typename std::vector<T>::iterator itr = v.begin();
            for (; itr != v.end(); itr++) {
                std::cout << " "<< *itr;
            }
            std::cout <<  " ]\n";
        }

    enum init_type_t { RANDOM, FORGY, PLUSPLUS } g_init; // May have to use
    enum kmspp_stage_t { ADDMEAN, DIST } g_kmspp_stage; // Either adding a mean / computing dist
    enum kms_stage_t { INIT, ESTEP } g_stage; // What phase of the algo we're in

    static std::vector<cluster::ptr> g_clusters; // cluster means/centers

    // Helpers //
    static void print_clusters(std::vector<cluster::ptr>& clusters) {
        for (std::vector<cluster::ptr>::iterator it = clusters.begin();
                it != clusters.end(); ++it) {
            print_vector<double>((*it)->get_mean());
        }
        std::cout << "\n";
    }

    static void print_sample(vertex_id_t my_id, data_seq_iterator& count_it, edge_seq_iterator& id_it) {
        std::vector<std::string> v;
        while (count_it.has_next()) {
            edge_count e = count_it.next();
            vertex_id_t nid = id_it.next();
            char buffer [50];
            assert(sprintf(buffer, "%u:%i",nid, e.get_count()));
            v.push_back(std::string(buffer));
        }
        printf("V%u's vector: ", my_id); print_vector<std::string>(v);
    }
    // End helpers //

    class kmeans_vertex: public compute_vertex
    {
        unsigned cluster_id;
        double dist;

        public:
        kmeans_vertex(vertex_id_t id):
            compute_vertex(id) {
                dist = std::numeric_limits<double>::max(); // Start @ max
            }

        unsigned get_result() const {
            return cluster_id;
        }

        const vsize_t get_cluster_id() const {
            return cluster_id;
        }

        void run(vertex_program &prog) {
#if PRUNE
            if (g_even_iter) {
                if (this->dist < g_comp_thresh_v[cluster_id]) {
#if VERBOSE
                    printf("No work for v:%u with dist=%.3f < comp_thresh=%.3f\n",
                            prog.get_vertex_id(*this), this->dist,  g_comp_thresh_v[cluster_id]);
#endif
                    return; // No work for you!
                }
            }
#endif
            vertex_id_t id = prog.get_vertex_id(*this);
            request_vertices(&id, 1);
        }

        void run(vertex_program& prog, const page_vertex &vertex) {
            switch (g_stage) {
                case INIT:
                    run_init(prog, vertex, g_init);
                    break;
                case ESTEP:
                    run_distance(prog, vertex);
                    break;
                default:
                    assert(0);
            }
        }

        // Set a cluster to have the same mean as this sample
        void set_as_mean(const page_vertex &vertex, vertex_id_t my_id, unsigned to_cluster_id) {
            edge_seq_iterator id_it = vertex.get_neigh_seq_it(OUT_EDGE);
            data_seq_iterator count_it = ((const page_directed_vertex&)vertex).
                get_data_seq_it<edge_count>(OUT_EDGE);

            // Build the setter vector that we assign to a cluster center
            std::vector<double> setter;
            setter.assign(NUM_COLS, 0);
            while (id_it.has_next()) {
                vertex_id_t nid = id_it.next();
                edge_count e = count_it.next();
                setter[nid] = (double) e.get_count();
            }
            g_clusters[to_cluster_id]->set_mean(setter);
        }

        void run_on_message(vertex_program& prog, const vertex_message& msg) { }
        void run_init(vertex_program& prog, const page_vertex &vertex, init_type_t init);
        void run_distance(vertex_program& prog, const page_vertex &vertex);
        double get_distance(unsigned cl, edge_seq_iterator& id_it, data_seq_iterator& count_it);
    };

    /* Used in per thread cluster formation */
    class kmeans_vertex_program : public vertex_program_impl<kmeans_vertex>
    {
        std::vector<cluster::ptr> pt_clusters;
        unsigned pt_changed;
        std::vector<distpair> bounds; // std::pair<(double) min,(double) max> - 1 per cluster

        public:
        typedef std::shared_ptr<kmeans_vertex_program> ptr;

        //TODO: Optimize only add cluster when a vertex joins it
        kmeans_vertex_program() {
            for (unsigned thd = 0; thd < K; thd++) {
                pt_clusters.push_back(cluster::create(NUM_COLS));
                pt_changed = 0;
                bounds.push_back(distpair(std::numeric_limits<double>::max(),
                            std::numeric_limits<double>::min()));
            }
        }

        static ptr cast2(vertex_program::ptr prog) {
            return std::static_pointer_cast<kmeans_vertex_program, vertex_program>(prog);
        }

        const std::vector<cluster::ptr>& get_pt_clusters() {
            return pt_clusters;
        }

        void add_member(const unsigned id, edge_seq_iterator& id_it, data_seq_iterator& count_it) {
            pt_clusters[id]->add_member(id_it, count_it);
        }

        void update_bounds(const unsigned cluster_id, const double val) {
            if (val < bounds[cluster_id].first)
                bounds[cluster_id].first = val;
            if (val > bounds[cluster_id].second)
                bounds[cluster_id].second = val;
        }

        const unsigned get_pt_changed() {
            return pt_changed;
        }

       void pt_changed_pp() {
           pt_changed++;
       }

       std::vector<distpair>& get_bounds() {
           return bounds;
       }
    };

    class kmeans_vertex_program_creater: public vertex_program_creater
    {
        public:
            vertex_program::ptr create() const {
                return vertex_program::ptr(new kmeans_vertex_program());
            }
    };

    /* Used in kmeans++ initialization */
    class kmeanspp_vertex_program : public vertex_program_impl<kmeans_vertex>
    {
        double pt_cuml_sum;

        public:
        typedef std::shared_ptr<kmeanspp_vertex_program> ptr;

        kmeanspp_vertex_program() {
            pt_cuml_sum = 0;
        }

        static ptr cast2(vertex_program::ptr prog) {
            return std::static_pointer_cast<kmeanspp_vertex_program, vertex_program>(prog);
        }

        void pt_cuml_sum_peq (double val) {
            pt_cuml_sum += val;
        }

        const double get_pt_cuml_sum() const {
            return pt_cuml_sum;
        }
    };

    class kmeanspp_vertex_program_creater: public vertex_program_creater
    {
        public:
            vertex_program::ptr create() const {
                return vertex_program::ptr(new kmeanspp_vertex_program());
            }
    };

    void kmeans_vertex::run_init(vertex_program& prog, const page_vertex &vertex, init_type_t init) {
        switch (g_init) {
            case RANDOM:
                {
                    unsigned new_cluster_id = random() % K;
                    kmeans_vertex_program& vprog = (kmeans_vertex_program&) prog;
#if KM_TEST
                    printf("Random init: v%u assigned to cluster: c%x\n", prog.get_vertex_id(*this), new_cluster_id);
#endif
                    this->cluster_id = new_cluster_id;
                    edge_seq_iterator id_it = vertex.get_neigh_seq_it(OUT_EDGE); // TODO: Make sure OUT_EDGE has data
                    data_seq_iterator count_it = ((const page_directed_vertex&)vertex).
                        get_data_seq_it<edge_count>(OUT_EDGE); //TODO: Make sure we have a directed graph
                    vprog.add_member(cluster_id, id_it, count_it);
                }
                break;
            case FORGY:
                {
                    vertex_id_t my_id = prog.get_vertex_id(*this);
#if KM_TEST
                    printf("Forgy init: v%u setting cluster: c%x\n", my_id, g_init_hash[my_id]);
#endif
                    set_as_mean(vertex, my_id, g_init_hash[my_id]);
                }
                break;
            case PLUSPLUS:
                {
                    edge_seq_iterator id_it = vertex.get_neigh_seq_it(OUT_EDGE);
                    data_seq_iterator count_it = ((const page_directed_vertex&)vertex).
                        get_data_seq_it<edge_count>(OUT_EDGE);

                    if (g_kmspp_stage == ADDMEAN) {
#if KM_TEST
                        vertex_id_t my_id = prog.get_vertex_id(*this);
                        printf("kms++ v%u making itself c%u\n", my_id, g_kmspp_cluster_idx);
#endif
                        g_clusters[g_kmspp_cluster_idx]->add_member(id_it, count_it);
                    } else {
                        // TODO: Test putting if (my_id != g_kmspp_next_cluster) test
                        vertex_id_t my_id = prog.get_vertex_id(*this);
                        double dist = get_distance(g_kmspp_cluster_idx, id_it, count_it);
                        if (dist < g_kmspp_distance[my_id]) {
#if VERBOSE
                            printf("kms++ v%u updating dist from: %.3f to %.3f\n", my_id, g_kmspp_distance[my_id], dist);
#endif
                            g_kmspp_distance[my_id] = dist;
                        }
                    }
                }
                break;
            default:
                assert(0);
        }
    }

    double kmeans_vertex::get_distance(unsigned cl, edge_seq_iterator& id_it,
         data_seq_iterator& count_it) {
        double dist = 0;
        double diff;
        while(id_it.has_next()) {
            vertex_id_t nid = id_it.next();
            edge_count e = count_it.next();
            diff = e.get_count() - (*g_clusters[cl])[nid]; // TODO: Do we need to take the abs value here?
            dist += diff*diff;
        }
        return dist;
    }

    void kmeans_vertex::run_distance(vertex_program& prog, const page_vertex &vertex) {
        kmeans_vertex_program& vprog = (kmeans_vertex_program&) prog;
        double best = std::numeric_limits<double>::max();
        unsigned new_cluster_id = INVALID_CLUST_ID;

        for (unsigned cl = 0; cl < K; cl++) {
            // TODO: Better access pattern than getting a new iterator every time
            edge_seq_iterator id_it = vertex.get_neigh_seq_it(OUT_EDGE);
            data_seq_iterator count_it =
                ((const page_directed_vertex&)vertex).get_data_seq_it<edge_count>(OUT_EDGE);

            double dist = get_distance(cl, id_it, count_it);
            if (dist < best) { // Get the distance to cluster `cl'
                new_cluster_id = cl;
                best = dist;
            }
        }

        assert(new_cluster_id >= 0 && new_cluster_id < K);
        if (new_cluster_id != this->cluster_id) {
#if VERBOSE
            vertex_id_t my_id = prog.get_vertex_id(*this);
            edge_seq_iterator id_it = vertex.get_neigh_seq_it(OUT_EDGE);
            data_seq_iterator count_it =
                ((const page_directed_vertex&)vertex).get_data_seq_it<edge_count>(OUT_EDGE);
            printf("Vertex%u changed membership from c%u to c%u with best-dist: %.4f\n",
                    my_id, cluster_id, new_cluster_id, best);
            print_sample(my_id, count_it, id_it);
#endif
            vprog.pt_changed_pp(); // Add a vertex to the count of changed ones
        }
        this->cluster_id = new_cluster_id;

        edge_seq_iterator id_it = vertex.get_neigh_seq_it(OUT_EDGE);
        data_seq_iterator count_it =
            ((const page_directed_vertex&)vertex).get_data_seq_it<edge_count>(OUT_EDGE);
        vprog.add_member(cluster_id, id_it, count_it);
        vprog.update_bounds(cluster_id, best);
        this->dist = best;
    }

    static FG_vector<unsigned>::ptr get_membership(graph_engine::ptr mat) {
        FG_vector<unsigned>::ptr vec = FG_vector<unsigned>::create(mat);
        mat->query_on_all(vertex_query::ptr(new save_query<unsigned, kmeans_vertex>(vec)));
        return vec;
    }

    static void clear_clusters() {
        for (unsigned thd = 0; thd < g_clusters.size(); thd++) {
            g_clusters[thd]->clear();
        }
    }

    static void update_clusters(graph_engine::ptr mat, std::vector<unsigned>& num_members_v) {
        clear_clusters();
        std::vector<vertex_program::ptr> kms_clust_progs;
        mat->get_vertex_programs(kms_clust_progs);
#if PRUNE
        std::vector<distpair> bounds;
        std::vector<distpair> tmp_g_bounds;
#endif

        for (unsigned thd = 0; thd < kms_clust_progs.size(); thd++) {
            kmeans_vertex_program::ptr kms_prog = kmeans_vertex_program::cast2(kms_clust_progs[thd]);
            std::vector<cluster::ptr> pt_clusters = kms_prog->get_pt_clusters();
#if PRUNE
            if (!g_even_iter) {
                bounds = kms_prog->get_bounds();
                if (thd == 0) { tmp_g_bounds = bounds; } // Avoid loc to init tmp_g_bounds
            }
#endif
            g_num_changed += kms_prog->get_pt_changed();
            /* Merge the per-thread clusters */
            for (unsigned cl = 0; cl < K; cl++) {
                *(g_clusters[cl]) += *(pt_clusters[cl]);
#if PRUNE
                if (!g_even_iter) {
                    if (bounds[cl].first < tmp_g_bounds[cl].first) {
                        tmp_g_bounds[cl].first = bounds[cl].first;
                    }
                    if (bounds[cl].second > tmp_g_bounds[cl].second) {
                        tmp_g_bounds[cl].second = bounds[cl].second;
                    }
                }
#endif
                if (thd == kms_clust_progs.size()-1) {
                    g_clusters[cl]->finalize();
                    num_members_v[cl] = g_clusters[cl]->get_num_members();
                }
            }
        }
#if PRUNE
        if (!g_even_iter) {
            for (unsigned cl = 0; cl < K; cl++) {
                g_comp_thresh_v[cl] = (fabs(tmp_g_bounds[cl].first -
                            tmp_g_bounds[cl].second) * g_COMP_THRESH) + tmp_g_bounds[cl].first;
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
            kmeanspp_vertex_program::ptr kmspp_prog = kmeanspp_vertex_program::cast2(vprog);
            cuml_sum += kmspp_prog->get_pt_cuml_sum();
        }
        cuml_sum = (cuml_sum * ((double)random())) / (RAND_MAX-1.0);

        g_kmspp_cluster_idx++;

        for (unsigned row = 0; row < NUM_ROWS; row++) {
            cuml_sum -= g_kmspp_distance[row];
            if (cuml_sum <= 0) {
#if KM_TEST
                BOOST_LOG_TRIVIAL(info) << "Choosing v:" << row << " as center K = " << g_kmspp_cluster_idx;
#endif
                return row;
            }
        }
        exit(EXIT_FAILURE);
    }

    // Return all the cluster means only
    static void get_means(std::vector<std::vector<double>>& means) {
        for (std::vector<cluster::ptr>::iterator it = g_clusters.begin();
                it != g_clusters.end(); ++it) {
           means.push_back((*it)->get_mean());
        }
    }
}

namespace fg
{
    sem_kmeans_ret::ptr compute_sem_kmeans(FG_graph::ptr fg, const size_t k, const std::string init,
            const unsigned max_iters, const double tolerance, const double comp_thresh) {
#ifdef PROFILER
        ProfilerStart("/home/disa/FlashGraph/flash-graph/libgraph-algs/sem_kmeans.perf");
#endif
        graph_index::ptr index = NUMA_graph_index<kmeans_vertex>::create(
                fg->get_graph_header());
        graph_engine::ptr mat = fg->create_engine(index);

        K = k;
        NUM_ROWS = mat->get_max_vertex_id() + 1;
        NUM_COLS = NUM_ROWS;

#if KM_TEST
        BOOST_LOG_TRIVIAL(info) << "We have rows = " << NUM_ROWS << ", cols = " <<
           NUM_COLS;
#endif

        // Check k
        if (K > NUM_ROWS || K < 2 || K == (unsigned)-1) {
            BOOST_LOG_TRIVIAL(fatal)
                << "'k' must be between 2 and the number of rows in the matrix" <<
                "k = " << K;
            exit(EXIT_FAILURE);
        }

        // Check Initialization
        if (init.compare("random") && init.compare("kmeanspp") &&
                init.compare("forgy")) {
            BOOST_LOG_TRIVIAL(fatal)
                << "[ERROR]: param init must be one of: 'random', 'forgy', 'kmeanspp'.It is '"
                << init << "'";
            exit(EXIT_FAILURE);
        }

        gettimeofday(&start , NULL);
        /*** Begin VarInit of data structures ***/
        FG_vector<unsigned>::ptr cluster_assignments; // Which cluster a sample is in
        for (size_t cl = 0; cl < k; cl++)
            g_clusters.push_back(cluster::create(NUM_COLS));

#if PRUNE
        g_comp_thresh_v.assign(K, 0); // Init g_COMP_THRESH for each cluster to 0
        g_COMP_THRESH = comp_thresh;
#endif
        std::vector<unsigned> num_members_v;
        num_members_v.resize(K);
        /*** End VarInit ***/
        g_stage = INIT;

        if (init == "random") {
            BOOST_LOG_TRIVIAL(info) << "Running init: '"<< init <<"' ...";
            g_init = RANDOM;
            mat->start_all(vertex_initializer::ptr(),
                    vertex_program_creater::ptr(new kmeans_vertex_program_creater()));
            mat->wait4complete();

            update_clusters(mat, num_members_v);
        }
        if (init == "forgy") {
            BOOST_LOG_TRIVIAL(info) << "Deterministic Init is: '"<< init <<"'";
            g_init = FORGY;

            // Select K in range NUM_ROWS
            std::vector<vertex_id_t> init_ids; // Used to start engine
            for (unsigned cl = 0; cl < K; cl++) {
                vertex_id_t id = random() % NUM_ROWS;
                g_init_hash[id] = cl; // <vertex_id, cluster_id>
                init_ids.push_back(id);
            }
            mat->start(&init_ids.front(), K);
            mat->wait4complete();

        } else if (init == "kmeanspp") {
            BOOST_LOG_TRIVIAL(info) << "Init is '"<< init <<"'";
            g_init = PLUSPLUS;

            // Init g_kmspp_distance to max distance
            g_kmspp_distance.assign(NUM_ROWS, std::numeric_limits<double>::max());

            g_kmspp_cluster_idx = 0;
            g_kmspp_next_cluster = random() % NUM_ROWS;
            g_kmspp_distance[g_kmspp_next_cluster] = 0;

            // Fire up K engines with 2 iters/engine
            while (true) {
                // Start 1 vertex which will activate all
                g_kmspp_stage = ADDMEAN;
                mat->start(&g_kmspp_next_cluster, 1, vertex_initializer::ptr(),
                    vertex_program_creater::ptr(new kmeanspp_vertex_program_creater()));
                mat->wait4complete();
#if KM_TEST
                BOOST_LOG_TRIVIAL(info) << "Printing clusters after sample set_mean ...";
                print_clusters(g_clusters);
#endif
                if (g_kmspp_cluster_idx+1 == K) { break; } // skip the distance comp since we picked clusters
                g_kmspp_stage = DIST;
                mat->start_all(); // Only need a vanilla vertex_program
                mat->wait4complete();
#if KM_TEST
                BOOST_LOG_TRIVIAL(info) << "Printing updated distances";
                print_vector<double>(g_kmspp_distance);
#endif

                g_kmspp_next_cluster = kmeanspp_get_next_cluster_id(mat);
            }
        }

#if KM_TEST
        BOOST_LOG_TRIVIAL(info) << "Printing cluster means:";
        print_clusters(g_clusters);
        BOOST_LOG_TRIVIAL(info) << "Printing cluster assignments:";
        get_membership(mat)->print(NUM_COLS);
#endif

        g_stage = ESTEP;
        BOOST_LOG_TRIVIAL(info) << "SEM-K||means starting ...";

        bool converged = false;

        std::string str_iters = max_iters == std::numeric_limits<unsigned>::max() ?
            "until convergence ...":
            std::to_string(max_iters) + " iterations ...";
        BOOST_LOG_TRIVIAL(info) << "Computing " << str_iters;
        g_iter = 1;

        while (g_iter < max_iters) {
            g_even_iter = true ? (g_iter % 2 == 0) : false;
            BOOST_LOG_TRIVIAL(info) << "E-step Iteration " << g_iter <<
                " . Computing cluster assignments ...";

            mat->start_all(vertex_initializer::ptr(),
                    vertex_program_creater::ptr(new kmeans_vertex_program_creater()));
            mat->wait4complete();
            BOOST_LOG_TRIVIAL(info) << "M-step Updating cluster means ...";
            update_clusters(mat, num_members_v);

#if KM_TEST
            BOOST_LOG_TRIVIAL(info) << "Printing cluster means:";
            print_clusters(g_clusters);
            BOOST_LOG_TRIVIAL(info) << "Getting cluster membership ...";
            get_membership(mat)->print(NUM_COLS);
            BOOST_LOG_TRIVIAL(info) << "** Samples changes cluster: " << g_num_changed << " **\n";
#endif

            if (g_num_changed == 0 || ((g_num_changed/(double)NUM_ROWS)) <= tolerance) {
                converged = true;
                break;
            } else {
                g_num_changed = 0;
            }
            g_iter++;
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
                "K-means converged in " << g_iter << " iterations";
        } else {
            BOOST_LOG_TRIVIAL(warning) << "[Warning]: K-means failed to converge in "
                << g_iter << " iterations";
        }
        BOOST_LOG_TRIVIAL(info) << "\n******************************************\n";

        print_vector<unsigned>(num_members_v);

        std::vector<std::vector<double>> means;
        get_means(means);
        cluster_assignments = get_membership(mat);

        return sem_kmeans_ret::create(cluster_assignments, means, num_members_v, g_iter);
    }
}
