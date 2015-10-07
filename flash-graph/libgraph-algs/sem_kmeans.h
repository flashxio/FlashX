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

using namespace fg;

namespace {
    typedef safs::page_byte_array::seq_const_iterator<edge_count> data_seq_iterator;
    class cluster
    {
        private:
            std::vector<double> mean; // cluster mean
            unsigned num_members; // cluster assignment count
            bool complete; // Have we already divided by num_members

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
            }

            cluster(const std::vector<double>& mean) {
                set_mean(mean);
                num_members = 1;
                complete = true;
            }

        public:
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

            void assign(const double val) {
                mean.assign(mean.size(), val);
            }

            void clear() {
                this->assign(0);
                this->num_members = 0;
                complete = false;
            }

            const std::vector<double>& get_mean() const {
                return mean;
            }

            void set_mean(const std::vector<double>& mean) {
                this->mean = mean;
            }

            const unsigned get_num_members() const {
                return num_members;
            }

            const bool is_complete() const {
                return complete;
            }

            const unsigned size() const {
                return mean.size();
            }

            void finalize() {
                assert(!complete);
                this->div(this->num_members);
            }

            void add_member(edge_seq_iterator& id_it, data_seq_iterator& count_it) {
                while(id_it.has_next()) {
                    vertex_id_t nid = id_it.next();
                    edge_count e = count_it.next();
                    mean[nid] += e.get_count();
                }
                num_members++;
            }

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
            const unsigned max_iters, const double tolerance, const double comp_thresh=0);
}
#endif
