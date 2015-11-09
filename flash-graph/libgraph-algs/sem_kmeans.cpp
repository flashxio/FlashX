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

#define KM_TEST 1
#define VERBOSE 0
#define INVALID_CLUST_ID -1

// TODO: Opt Assign cluster ID to kms++ selected vertices in ADDMEAN phase

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

    enum dist_type_t { EUCL, COS }; // Euclidean, Cosine distance
    enum init_type_t { RANDOM, FORGY, PLUSPLUS } g_init; // May have to use
    enum kmspp_stage_t { ADDMEAN, DIST } g_kmspp_stage; // Either adding a mean / computing dist
    enum kms_stage_t { INIT, ESTEP } g_stage; // What phase of the algo we're in

    static std::vector<cluster::ptr> g_clusters; // cluster means/centers

    // Begin Helpers //
    template <typename T>
        static void print_vector(typename std::vector<T> v) {
            std::cout << "[";
            typename std::vector<T>::iterator itr = v.begin();
            for (; itr != v.end(); itr++) {
                std::cout << " "<< *itr;
            }
            std::cout <<  " ]\n";
        }

    static void print_clusters(std::vector<cluster::ptr>& clusters) {
        for (std::vector<cluster::ptr>::iterator it = clusters.begin();
                it != clusters.end(); ++it) {
            print_vector<double>((*it)->get_mean());
        }
        std::cout << "\n";
    }

    static void print_sample(vertex_id_t my_id, data_seq_iter& count_it, edge_seq_iterator& id_it) {
        std::vector<std::string> v;
        vertex_id_t nid = 0;
        while (count_it.has_next()) {
            char buffer [50];
#ifdef PAGE_ROW
            double e = count_it.next();
            assert(sprintf(buffer, "%u:%.3f", nid++, e));
#else
            edge_count e = count_it.next();
            vertex_id_t nid = id_it.next();
            assert(sprintf(buffer, "%u:%i",nid, e.get_count()));
#endif
            v.push_back(std::string(buffer));
        }
        printf("V%u's vector: \n", my_id); print_vector<std::string>(v);
    }
    // End helpers //

#if 0
    static double const eucl_dist(const cluster::ptr l_clust, const cluster::ptr r_clust) {
        double dist = 0;
        BOOST_VERIFY(l_clust->size() == r_clust->size());

        for (unsigned col = 0; col < NUM_COLS; col++) {
            double diff = (*l_clust)[col] - (*r_clust)[col];
            dist += diff * diff;
        }
        return  dist;
    }
#else
    template <typename T>
        static double const eucl_dist(const T* lhs, const T* rhs) {
            double dist = 0;

            BOOST_VERIFY(lhs->size() == rhs->size());

            for (unsigned col = 0; col < lhs->size(); col++) {
                double diff = (*lhs)[col] - (*rhs)[col];
                dist += diff * diff;
            }
            return  dist;
        }
#endif

#if PRUNE
    // NOTE: Creates a matrix like this e.g for K = 5
    /* - Don't store full matrix, don't store dist to myself -> space: (k/2)-k
       0 ==> 1 2 3 4
       1 ==> 2 3 4
       2 ==> 3 4
       3 ==> 4
       */
    class dist_matrix
    {
        private:
            std::vector<std::vector<double>> mat;
            unsigned rows;

            dist_matrix(const unsigned rows) {
                BOOST_VERIFY(rows > 1);

                this->rows = rows-1;
                // Distance to everyone other than yourself
                for (unsigned i = this->rows; i > 0; i--) {
                    std::vector<double> dist_row;
                    dist_row.assign(i, std::numeric_limits<double>::max());
                    mat.push_back(dist_row);
                }
            }

        public:
            typedef typename std::shared_ptr<dist_matrix> ptr;

            static ptr create(const unsigned rows) {
                return ptr(new dist_matrix(rows));
            }

            void translate(unsigned& row, unsigned& col) {
                // First make sure the smaller is the row
                if (row > col) {
                    std::swap(row, col);
                }

                BOOST_VERIFY(row < rows);
                col = col - row - 1; // Translation
                BOOST_VERIFY(col < (rows - row));
            }

            /* Do a translation from raw id's to indexes in the distance matrix */
            double get(unsigned row, unsigned col) {
                if (row == col) { return 0; }
                translate(row, col);
#if VERBOSE
                BOOST_LOG_TRIVIAL(info) << "Getting row=" << row << "col=" << col;
#endif
                return mat[row][col];
            }

            void set(unsigned row, unsigned col, double val) {
                BOOST_VERIFY(row != col);
                translate(row, col);
                mat[row][col] = val;
            }

            void print() {
                for (unsigned row = 0; row < rows; row++) {
                    std::cout << row << " ==> ";
                    print_vector<double>(mat[row]);
                }
            }
    };

    static bool g_prune_init = false;
    static dist_matrix::ptr g_cluster_dist;
#endif

    class kmeans_vertex: public compute_vertex
    {
        unsigned cluster_id;
        double dist;
        std::vector<double> lwr_bnd;
        double uppr_bnd;
#if PRUNE
        unsigned nxt_clstr; // If prune fails .. which cluster do we start from?
        bool r;
        bool in_x_prev_iter;
#endif
        public:
        kmeans_vertex(vertex_id_t id):
            compute_vertex(id) {
                dist = std::numeric_limits<double>::max(); // Start @ max
#if PRUNE
                uppr_bnd = dist;
                lwr_bnd.resize(K); // Set K items to 0
                nxt_clstr = 0;
                r = false;
                in_x_prev_iter = false;
#endif
            }

#if PRUNE
        void set_in_x_prev_iter(bool in=true) {
            in_x_prev_iter = false;
        }

        bool const get_in_x_prev_iter() {
            return in_x_prev_iter;
        }

#endif
        unsigned get_result() const {
            return cluster_id;
        }

        const vsize_t get_cluster_id() const {
            return cluster_id;
        }

        void run(vertex_program &prog);

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
#ifndef PAGE_ROW
            edge_seq_iterator id_it = vertex.get_neigh_seq_it(OUT_EDGE);
#endif

            vertex_id_t nid = 0;
            data_seq_iter count_it = ((const page_row&)vertex).
                get_data_seq_it<double>();

            // Build the setter vector that we assign to a cluster center
            std::vector<double> setter;
            setter.assign(NUM_COLS, 0);
            while (count_it.has_next()) {
#ifdef PAGE_ROW
                double e = count_it.next();
                setter[nid++] = e;
#else
                edge_count e = count_it.next();
                nid = id_it.next();
                setter[nid] = (double) e.get_count();
#endif
            }
            g_clusters[to_cluster_id]->set_mean(setter);
        }

        void run_on_message(vertex_program& prog, const vertex_message& msg) { }
        void run_init(vertex_program& prog, const page_vertex &vertex, init_type_t init);
        void run_distance(vertex_program& prog, const page_vertex &vertex);
        double get_distance(unsigned cl, edge_seq_iterator& id_it, data_seq_iter& count_it);
    };

    /* Used in per thread cluster formation */
    class kmeans_vertex_program : public vertex_program_impl<kmeans_vertex>
    {
        std::vector<cluster::ptr> pt_clusters;
        unsigned pt_changed;

        public:
        typedef std::shared_ptr<kmeans_vertex_program> ptr;

        //TODO: Optimize only add cluster when a vertex joins it
        kmeans_vertex_program() {
            for (unsigned thd = 0; thd < K; thd++) {
                pt_clusters.push_back(cluster::create(NUM_COLS));
                pt_changed = 0;
            }
        }

        static ptr cast2(vertex_program::ptr prog) {
            return std::static_pointer_cast<kmeans_vertex_program, vertex_program>(prog);
        }

        const std::vector<cluster::ptr>& get_pt_clusters() {
            return pt_clusters;
        }

        void add_member(const unsigned id, edge_seq_iterator& id_it, data_seq_iter& count_it) {
            pt_clusters[id]->add_member(id_it, count_it);
        }

#if PRUNE
        void remove_member(const unsigned id, edge_seq_iterator& id_it, data_seq_iter& count_it) {
            pt_clusters[id]->remove_member(id_it, count_it);
        }
#endif

        const unsigned get_pt_changed() {
            return pt_changed;
        }

        void pt_changed_pp() {
            pt_changed++;
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


    void kmeans_vertex::run(vertex_program &prog) {
#if KM_TEST
        vertex_id_t my_id = prog.get_vertex_id(*this);
#endif
#if PRUNE
        this->nxt_clstr = 0;

        if (g_prune_init) {
            for (nxt_clstr = 0; nxt_clstr < K; nxt_clstr++) {
                // Prunes `some cluster subset' but still needs to
                //  do I/O because `some other subset' of clusters needs to be computed :(
                // Use Lemma 1 to prune
                if (g_cluster_dist->get(this->cluster_id, nxt_clstr) < (2*this->dist)) {
#if KM_TEST
                    printf("Lemma 1 failed for v:%u at nxt_clstr = %u in init\n", my_id, nxt_clstr);
#endif
                    break; // Because d(this,cl) >= d(this,this->cluster_id)
                }
            }
        } else {
            if (get_in_x_prev_iter()) {
                // #5
                for (unsigned cl = 0; cl < K; cl++) {
                    std::max((lwr_bnd[cl] - g_clusters[cl]->get_prev_dist()), 0.0);
                }
#if 1
                BOOST_LOG_TRIVIAL(info) << "After #5 v:" << my_id << " lwr_bnd = ";
                print_vector(lwr_bnd);
#endif
                // #6
                uppr_bnd = uppr_bnd + g_clusters[cluster_id]->get_prev_dist();
#if 1
                BOOST_LOG_TRIVIAL(info) << "After #5 v:" << my_id << " uppr_bnd = " << uppr_bnd;
#endif
                r = true;
            } else if (uppr_bnd <= g_clusters[cluster_id]->get_s_val()) {
                set_in_x_prev_iter();
                // These may still be prunable
                for (nxt_clstr = 0; nxt_clstr < K; nxt_clstr++) {
                    if (nxt_clstr != cluster_id) { // #3
                        if (uppr_bnd > lwr_bnd[nxt_clstr] && uppr_bnd >
                                (g_cluster_dist->get(cluster_id, nxt_clstr)*.5)) {
                            // #3a
                            if (r) {
                                r = false;
#if KM_TEST
                                BOOST_LOG_TRIVIAL(info) << "3a failed for v:" << my_id
                                    << " at nxt_clstr = " << nxt_clstr;
#endif
                                break; // do I/O request and compute d(x,c(x))
                            } else {
                                BOOST_VERIFY(dist >= uppr_bnd);
                                this->dist = uppr_bnd; // d(x,c(x)) = u(x)
                            }

                            // #3b
                            if (this->dist > lwr_bnd[nxt_clstr] ||
                                    this->dist > (.5*g_cluster_dist->get(this->cluster_id,nxt_clstr))) {
#if KM_TEST
                                BOOST_LOG_TRIVIAL(info) << "3b failed for v:" << my_id
                                    << " at nxt_clstr = " << nxt_clstr;
#endif
                                break; // do I/O request and compute this d(x,c)
#if 0
                                double _dist = get_distance(nxt_clstr, id_it, count_it);
                                if (_dist < this->dist or best) {
                                    this->cluster_id = nxt_clstr;
                                }
#endif
                            }
                        }
                    }
                }
            } else {
                set_in_x_prev_iter(false);
            }
        }

        // No I/O
        if (nxt_clstr == K) {
#if KM_TEST
            printf("No I/O for v:%u\n", my_id);
#endif
            return;
        }
#endif

#if KM_TEST
        printf("I/O for v:%u\n", my_id);
#endif
        vertex_id_t id = prog.get_vertex_id(*this);
        request_vertices(&id, 1);
    }


    void kmeans_vertex::run_init(vertex_program& prog, const page_vertex &vertex, init_type_t init) {
        switch (g_init) {
            case RANDOM:
                {
                    unsigned new_cluster_id = random() % K;
                    kmeans_vertex_program& vprog = (kmeans_vertex_program&) prog;
#if KM_TEST
                    printf("Random init: v%u assigned to cluster: c%x\n",
                            prog.get_vertex_id(*this), new_cluster_id);
#endif
                    this->cluster_id = new_cluster_id;
                    edge_seq_iterator id_it = vertex.get_neigh_seq_it(OUT_EDGE);
                    data_seq_iter count_it = ((const page_row&)vertex).
                        get_data_seq_it<double>();
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
                    data_seq_iter count_it = ((const page_row&)vertex).
                        get_data_seq_it<double>();

                    if (g_kmspp_stage == ADDMEAN) {
#ifndef PAGE_ROW
#if KM_TEST
                        vertex_id_t my_id = prog.get_vertex_id(*this);
                        printf("kms++ v%u making itself c%u\n", my_id, g_kmspp_cluster_idx);
#endif
#endif
                        g_clusters[g_kmspp_cluster_idx]->add_member(id_it, count_it);
                    } else {
                        // FIXME: Opt Test putting if (my_id != g_kmspp_next_cluster) test
                        vertex_id_t my_id = prog.get_vertex_id(*this);
                        double dist = get_distance(g_kmspp_cluster_idx, id_it, count_it);
                        if (dist < g_kmspp_distance[my_id]) {
#if VERBOSE
                            printf("kms++ v%u updating dist from: %.3f to %.3f\n",
                                    my_id, g_kmspp_distance[my_id], dist);
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
            data_seq_iter& count_it) {
        double dist = 0;
        double diff;
        vertex_id_t nid = 0;

        while(count_it.has_next()) {
#ifdef PAGE_ROW
            double e = count_it.next();
            diff = e - (*g_clusters[cl])[nid++];
#else
            nid = id_it.next();
            edge_count e = count_it.next();
            diff = e.get_count() - (*g_clusters[cl])[nid];
#endif
            dist += diff*diff;
        }
        return dist;
    }

    void kmeans_vertex::run_distance(vertex_program& prog, const page_vertex &vertex) {
        kmeans_vertex_program& vprog = (kmeans_vertex_program&) prog;
        double best = std::numeric_limits<double>::max();
        unsigned new_cluster_id = cluster_id;

#if PRUNE
#if KM_TEST
        vertex_id_t my_id = prog.get_vertex_id(*this);
        BOOST_LOG_TRIVIAL(info) << "Restarting v:" << my_id << " prune from cl:" << nxt_clstr;
#endif
        for (unsigned cl = nxt_clstr; cl < K; cl++) {
            // Prunes `some cluster subset' but still needs to
            //  do I/O because `some other subset' of clusters needs to be computed :(
            if (g_prune_init) {
                // Use Lemma 1 to prune
                if (g_cluster_dist->get(this->cluster_id, cl) >= (2*this->dist)) {
                    continue; // Because d(this,cl) >= d(this,this->cluster_id)
                }
            }
#else
            for (unsigned cl = 0; cl < K; cl++) {
#endif
                // TODO: Better access pattern than getting a new iterator every time
                edge_seq_iterator id_it = vertex.get_neigh_seq_it(OUT_EDGE);
                data_seq_iter count_it =
                    ((const page_row&)vertex).get_data_seq_it<double>();

                double dist = get_distance(cl, id_it, count_it);
                if (dist < best) { // Get the distance to cluster `cl'
                    new_cluster_id = cl;
                    best = dist;
                }
#if PRUNE
                if (g_prune_init) {
                    lwr_bnd[cl] = dist;
                    uppr_bnd = best;
                }
#endif
            }

            BOOST_VERIFY(new_cluster_id >= 0 && new_cluster_id < K);
            // TODO: Optmize by making local
            edge_seq_iterator id_it = vertex.get_neigh_seq_it(OUT_EDGE);
            data_seq_iter count_it = ((const page_row&)vertex).get_data_seq_it<double>();

#if PRUNE
            if (g_prune_init) {
#if KM_TEST
                printf("v%u lwr_bnd = ", my_id); print_vector<double>(lwr_bnd);
#endif
                vprog.pt_changed_pp(); // Add a vertex to the count of changed ones
                this->cluster_id = new_cluster_id;
                vprog.add_member(cluster_id, id_it, count_it);
            } else if (new_cluster_id != this->cluster_id) {
#if KM_TEST
                printf("Vertex%u changed membership from c%u to c%u with best-dist: %.4f\n",
                        my_id, cluster_id, new_cluster_id, best);
#endif
                vprog.pt_changed_pp(); // Add a vertex to the count of changed ones
                vprog.remove_member(cluster_id, id_it, count_it);

                this->cluster_id = new_cluster_id;
                id_it = vertex.get_neigh_seq_it(OUT_EDGE);
                count_it = ((const page_row&)vertex).get_data_seq_it<double>();
                vprog.add_member(cluster_id, id_it, count_it);
            } else {
                // Not g_prune init & my_cluster remained the same .. Do nothing
#if KM_TEST
                printf("\tVertex%u stayed in c%u\n", my_id, cluster_id);
#endif
            }
#else
            if (this->cluster_id != new_cluster_id) {
#if KM_TEST
                printf("Vertex%u changed membership from c%u to c%u with best-dist: %.4f\n",
                        my_id, cluster_id, new_cluster_id, best);
#endif
                vprog.pt_changed_pp(); // Add a vertex to the count of changed ones
            }

            this->cluster_id = new_cluster_id;
            vprog.add_member(cluster_id, id_it, count_it);
#endif
            // Done by all : TODO: Verify I really need this
            this->dist = best;
        }

        static FG_vector<unsigned>::ptr get_membership(graph_engine::ptr mat) {
            FG_vector<unsigned>::ptr vec = FG_vector<unsigned>::create(mat);
            mat->query_on_all(vertex_query::ptr(new save_query<unsigned, kmeans_vertex>(vec)));
            return vec;
        }

        static void clear_clusters() {
            for (unsigned cl = 0; cl < g_clusters.size(); cl++) {
#if PRUNE
                g_clusters[cl]->set_prev_mean();

                if (g_prune_init) {
                    g_clusters[cl]->clear();
                } else {
                    g_clusters[cl]->unfinalize();
#if KM_TEST
                    std::cout << "Unfinalized g_clusters[thd] ==> ";
                    print_vector<double>(g_clusters[cl]->get_mean());
#endif
                }
#else
                g_clusters[cl]->clear();
#endif
            }
        }

        static void update_clusters(graph_engine::ptr mat, std::vector<unsigned>& num_members_v) {
            clear_clusters();
            std::vector<vertex_program::ptr> kms_clust_progs;
            mat->get_vertex_programs(kms_clust_progs);

            for (unsigned thd = 0; thd < kms_clust_progs.size(); thd++) {
                kmeans_vertex_program::ptr kms_prog = kmeans_vertex_program::cast2(kms_clust_progs[thd]);
                std::vector<cluster::ptr> pt_clusters = kms_prog->get_pt_clusters();
                g_num_changed += kms_prog->get_pt_changed();

                BOOST_VERIFY(g_num_changed <= NUM_ROWS);
                /* Merge the per-thread clusters */
                for (unsigned cl = 0; cl < K; cl++) {
#if KM_TEST
                    std::cout << "pt_clusters[" << cl << "] #" << pt_clusters[cl]->get_num_members() << " ==> ";
                    print_vector<double>(pt_clusters[cl]->get_mean());
#endif
                    *(g_clusters[cl]) += *(pt_clusters[cl]);
                    if (thd == kms_clust_progs.size()-1) {
                        g_clusters[cl]->finalize();
                        num_members_v[cl] = g_clusters[cl]->get_num_members();
#if PRUNE
                        double dist = eucl_dist<std::vector<double>>(&(g_clusters[cl]->get_mean()),
                                &(g_clusters[cl]->get_prev_mean()));
#if KM_TEST
                        BOOST_LOG_TRIVIAL(info) << "Distance to prev mean for c:"
                            << cl << " is " << dist;
#endif
                        g_clusters[thd]->set_prev_dist(dist);
#endif
                    }
                }
            }
        }

        // This ignore the cluster_id == cluster_id & any repetitive computations
#if PRUNE
        static void compute_cluster_dist() {
            //#pragma omp parallel for collapse(2) // FIXME: Coalese perhaps
            for (unsigned i = 0; i < K; i++) {
                for (unsigned j = i+1; j < K-1; j++) {
                    double dist = eucl_dist<std::vector<double>>(&(g_clusters[i]->get_mean()),
                            &(g_clusters[j]->get_mean()));
                    g_cluster_dist->set(i,j, dist);

                    // Set s(x) for each cluster
                    if (dist < g_clusters[i]->get_s_val()) {
                        g_clusters[i]->set_s_val(dist);
                    }

                    if (dist < g_clusters[j]->get_s_val()) {
                        g_clusters[j]->set_s_val(dist);
                    }
                }
            }

            // Serial
            for (unsigned cl = 0; cl < K; cl++) {
                g_clusters[cl]->finalize_s_val();
            }
        }
#endif

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
#if VERBOSE
                BOOST_LOG_TRIVIAL(info) << "cuml_sum = " << cuml_sum;
#endif
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
                const unsigned max_iters, const double tolerance, const double comp_thresh,
                const unsigned num_rows, const unsigned num_cols) {
#ifdef PROFILER
            ProfilerStart("/home/disa/FlashGraph/flash-graph/libgraph-algs/sem_kmeans.perf");
#endif
            K = k;

            // Check Initialization
            if (init.compare("random") && init.compare("kmeanspp") &&
                    init.compare("forgy")) {
                BOOST_LOG_TRIVIAL(fatal)
                    << "[ERROR]: param init must be one of: 'random', 'forgy', 'kmeanspp'.It is '"
                    << init << "'";
                exit(EXIT_FAILURE);
            }

            graph_index::ptr index = NUMA_graph_index<kmeans_vertex>::create(
                    fg->get_graph_header());
            graph_engine::ptr mat = fg->create_engine(index);

            NUM_ROWS = mat->get_max_vertex_id() + 1;
            NUM_COLS = num_cols;

            // Check k
            if (K > NUM_ROWS || K < 2 || K == (unsigned)-1) {
                BOOST_LOG_TRIVIAL(fatal)
                    << "'k' must be between 2 and the number of rows in the matrix " <<
                    "k = " << K;
                exit(EXIT_FAILURE);
            }

            BOOST_VERIFY(num_cols > 0);

#if KM_TEST
            BOOST_LOG_TRIVIAL(info) << "We have rows = " << NUM_ROWS << ", cols = " <<
                NUM_COLS;
#endif

            gettimeofday(&start , NULL);
            /*** Begin VarInit of data structures ***/
            FG_vector<unsigned>::ptr cluster_assignments; // Which cluster a sample is in
            for (size_t cl = 0; cl < k; cl++)
                g_clusters.push_back(cluster::create(NUM_COLS));

            std::vector<unsigned> num_members_v;
            num_members_v.resize(K);

#if PRUNE
            BOOST_LOG_TRIVIAL(info) << "Init of g_cluster_dist";
            // Distance to everyone other than yourself
            g_cluster_dist = dist_matrix::create(K);
#endif
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
#if KM_TEST
                BOOST_LOG_TRIVIAL(info) << "Assigning v:" << g_kmspp_next_cluster << " as first cluster";
#endif
                g_kmspp_distance[g_kmspp_next_cluster] = 0;

                // Fire up K engines with 2 iters/engine
                while (true) {
#if KM_TEST
                    BOOST_LOG_TRIVIAL(info) << "Printing updated distances";
                    print_vector<double>(g_kmspp_distance);
#endif
                    // TODO: Start 1 vertex which will activate all
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

                    g_kmspp_next_cluster = kmeanspp_get_next_cluster_id(mat);
                }
            }

#if PRUNE
            if (init == "forgy" || init == "kmeanspp") {
                g_prune_init = true; // set
                g_stage = ESTEP;
                BOOST_LOG_TRIVIAL(info) << "Init: Computing cluster distance matrix ...";
                compute_cluster_dist();
#if KM_TEST
                BOOST_LOG_TRIVIAL(info) << "Printing inited cluster distance matrix ...";
                g_cluster_dist->print();
#endif

                BOOST_LOG_TRIVIAL(info) << "Init: Running an engine for PRUNE since init is " << init;

                mat->start_all(vertex_initializer::ptr(),
                        vertex_program_creater::ptr(new kmeans_vertex_program_creater()));
                mat->wait4complete();
                BOOST_LOG_TRIVIAL(info) << "Init: M-step Updating cluster means ...";
#if KM_TEST
                BOOST_LOG_TRIVIAL(info) << "Init: Printing cluster means:";
                print_clusters(g_clusters);
#endif
                update_clusters(mat, num_members_v);
                g_prune_init = false; // reset
                g_num_changed = 0; // reset
            }
#endif

#if KM_TEST
            BOOST_LOG_TRIVIAL(info) << "Printing cluster means:";
            print_clusters(g_clusters);
            BOOST_LOG_TRIVIAL(info) << "Printing cluster assignments:";
            get_membership(mat)->print(NUM_ROWS);
#endif

            g_stage = ESTEP;
            BOOST_LOG_TRIVIAL(info) << "SEM-K||means starting ...";

            bool converged = false;

            std::string str_iters = max_iters == std::numeric_limits<unsigned>::max() ?
                "until convergence ...":
                std::to_string(max_iters) + " iterations ...";
            BOOST_LOG_TRIVIAL(info) << "Computing " << str_iters;
            unsigned iter = 1;

            while (iter < max_iters) {
                BOOST_LOG_TRIVIAL(info) << "E-step Iteration " << iter <<
                    " . Computing cluster assignments ...";
#if PRUNE
                BOOST_LOG_TRIVIAL(info) << "Main: Computing cluster distance matrix ...";
                compute_cluster_dist();
#endif
                mat->start_all(vertex_initializer::ptr(),
                        vertex_program_creater::ptr(new kmeans_vertex_program_creater()));
                mat->wait4complete();
                BOOST_LOG_TRIVIAL(info) << "Main: M-step Updating cluster means ...";
                update_clusters(mat, num_members_v);

#if KM_TEST
                BOOST_LOG_TRIVIAL(info) << "Printing cluster means:";
                print_clusters(g_clusters);
                BOOST_LOG_TRIVIAL(info) << "Getting cluster membership ...";
                get_membership(mat)->print(NUM_ROWS);
#endif
                BOOST_LOG_TRIVIAL(info) << "Printing cluster counts ...";
                print_vector<unsigned>(num_members_v);

                BOOST_LOG_TRIVIAL(info) << "** Samples changes cluster: " << g_num_changed << " **\n";

                if (g_num_changed == 0 || ((g_num_changed/(double)NUM_ROWS)) <= tolerance) {
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

            print_vector<unsigned>(num_members_v);

            std::vector<std::vector<double>> means;
            get_means(means);
            cluster_assignments = get_membership(mat);
#if KM_TEST
            BOOST_LOG_TRIVIAL(info) << "Printing updated distances";
            print_vector<double>(g_kmspp_distance);
#endif

            return sem_kmeans_ret::create(cluster_assignments, means, num_members_v, iter);
        }
    }
