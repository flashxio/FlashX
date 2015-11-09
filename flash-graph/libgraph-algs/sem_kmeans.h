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
#ifndef __SEM_KMEANS_H__
#define __SEM_KMEANS_H__

#include <math.h>

#include <vector>
#include <algorithm>
#include <map>

#include "graph_engine.h"
#include "graph_config.h"
#include "FGlib.h"
#include "save_result.h"
#define PAGE_ROW
#define INVALID_CLUST_ID -1

#define PRUNE 1

using namespace fg;

namespace {
    typedef safs::page_byte_array::seq_const_iterator<double> data_seq_iter;
    class cluster
    {
        private:
            std::vector<double> mean; // Cluster mean
            int num_members; // Cluster assignment count
            bool complete; // Have we already divided by num_members
#if PRUNE
            double s_val;
            std::vector<double> prev_mean; // For lwr_bnd
            double prev_dist; // Distance to prev mean
#endif

            void div(const unsigned val) {
                if (num_members > 0) {
                    for (unsigned i = 0; i < mean.size(); i++) {
                        mean[i] /= double(val);
                    }
                }
                complete = true;
            }

            cluster(const unsigned len) {
                mean.assign(len, 0);
                num_members = 0;
                complete = false;
#if PRUNE
                prev_mean.assign(len, 0);
#endif
            }

            cluster(const std::vector<double>& mean) {
                set_mean(mean);
                num_members = 1;
                complete = true;
#if PRUNE
                prev_mean.assign(mean.size(), 0);
#endif
            }

        public:
#if PRUNE
            void set_s_val(double val) {
                s_val = val;
            }

            double const get_s_val() { return s_val; }

            void finalize_s_val() {
                s_val /= 2.0;
            }

            const std::vector<double>& get_prev_mean() const {
                return prev_mean;
            }

            void set_prev_mean() {
                if (!is_complete()) {
                    BOOST_LOG_TRIVIAL(warning) << "WARNING: Doing nothing for "
                        "unfinalized mean. Permissible once";
                    return;
                }
                prev_mean = mean;
            }

            void set_prev_dist(double dist) {
                this->prev_dist = dist;
            }

            double get_prev_dist() {
                return this->prev_dist;
            }
#endif
            typedef typename std::shared_ptr<cluster> ptr;

            static ptr create(const unsigned len) {
                return ptr(new cluster(len));
            }

            static ptr create(const std::vector<double>& mean) {
                return ptr(new cluster(mean));
            }

            void init(const unsigned len) {
                mean.assign(len, 0);
                num_members = 0;
            }

            void clear() {
                std::fill(this->mean.begin(), this->mean.end(), 0);
                this->num_members = 0;
                complete = false;
            }

            const std::vector<double>& get_mean() const {
                return mean;
            }

            void set_mean(const std::vector<double>& mean) {
                this->mean = mean;
            }

            const int get_num_members() const {
                return num_members;
            }

            const bool is_complete() const {
                return complete;
            }

            const unsigned size() const {
                return mean.size();
            }

            void finalize() {
                if (is_complete()) {
                    BOOST_LOG_TRIVIAL(warning) << "WARNING: Calling finalize() on a"
                        " finalized object";
                    return;
                }
                this->div(this->num_members);
            }

            void unfinalize() {
                if (!is_complete()) {
                    BOOST_LOG_TRIVIAL(warning) << "WARNING: Calling unfinalize() on an"
                        " unfinalized object";
                    return;
                }
                complete = false;

                for (unsigned i = 0; i < size(); i++) {
                    this->mean[i] *= (double)num_members;
                }
            }

            void add_member(edge_seq_iterator& id_it, data_seq_iter& count_it) {
                vertex_id_t nid = 0;
                while(count_it.has_next()) {
#ifdef PAGE_ROW
                    double e = count_it.next();
                    mean[nid++] += e;
#else
                    nid = id_it.next();
                    edge_count e = count_it.next();
                    mean[nid] += e.get_count();
#endif
                }
                num_members++;
            }

#if PRUNE
            void remove_member(edge_seq_iterator& id_it, data_seq_iter& count_it) {
                vertex_id_t nid = 0;
                while(count_it.has_next()) {
#ifdef PAGE_ROW
                    double e = count_it.next();
                    mean[nid++] -= e;
#else
                    nid = id_it.next();
                    edge_count e = count_it.next();
                    mean[nid] -= e.get_count();
#endif
                }
                num_members--;
            }
#endif

            cluster& operator=(const cluster& other) {
                this->mean = other.get_mean();
                this->num_members = other.get_num_members();
                return *this;
            }

            double& operator[](const unsigned index) {
                assert(index < mean.size());
                return mean[index];
            }

            cluster& operator+=(cluster& rhs) {
                assert(rhs.size() == size());
                // TODO vectorize perhaps
                for (unsigned i = 0; i < mean.size(); i++) {
                    this->mean[i] += rhs[i];
                }
                this->num_members += rhs.get_num_members();
                return *this;
            }
    };
}

namespace fg
{
    // A class used a return object for R bindings
    class sem_kmeans_ret
    {
        private:
            FG_vector<unsigned>::ptr cluster_assignments;
            std::vector<std::vector<double>> centers;
            std::vector<unsigned> size;
            unsigned iters;

            sem_kmeans_ret(const FG_vector<unsigned>::ptr cluster_assignments,
                    const std::vector<std::vector<double>> centers,
                    const std::vector<unsigned>& size, const unsigned iters) {
                this->cluster_assignments = cluster_assignments;
                this->centers = centers;
                this->size = size;
                this->iters = iters;
            }

        public:
            typedef typename std::shared_ptr<sem_kmeans_ret> ptr;

            static ptr create(const FG_vector<unsigned>::ptr cluster_assignments,
                    const std::vector<std::vector<double>> centers,
                    const std::vector<unsigned>& size, const unsigned iters) {
                return ptr(new sem_kmeans_ret(cluster_assignments, centers, size, iters));
            }

            const FG_vector<unsigned>::ptr get_cluster_assignments() const {
                return this->cluster_assignments;
            }

            const unsigned get_iters() const {
                return this->iters;
            }

            const std::vector<unsigned>& get_size() const {
                return this->size;
            }

            const std::vector<std::vector<double>>& get_centers() const {
                return this->centers;
            }
    };

    /**
     * \brief Compute Semi-External Memory kmeans
     * \param fg The FlashGraph graph object for which you want to compute.
     * \param k The number of clusters.
     * \param init Initialization type [random, forgy, kmeanspp].
     * \param max_iters The max number of iterations to compute for.
     * \param tolerance The min fraction of changes from 1 iter to next required to converge.
     * \param comp_thresh Used to prune computation if specified. TODO: Explain.
     */
    sem_kmeans_ret::ptr compute_sem_kmeans(FG_graph::ptr fg, const size_t k, const std::string init,
            const unsigned max_iters, const double tolerance, const double comp_thresh=0,
            const unsigned num_rows=0, const unsigned num_cols=0);
}
#endif
