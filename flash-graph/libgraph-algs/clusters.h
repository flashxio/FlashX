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

#ifndef __CLUSTERS_H__
#define __CLUSTERS_H__

#include <limits>
#include <vector>
#include <memory>
#include <algorithm>

#include "log.h"

namespace km {
    typedef std::vector<double> kmsvector;
    typedef std::vector<double>::iterator kmsiterator;

    class clusters
    {
        friend class prune_clusters;
        private:
            // Together are nXd matrix
            unsigned ncol;
            unsigned nclust;
            std::vector<int> num_members_v; // Cluster assignment counts
            std::vector<bool> complete_v; // Have we already divided by num_members

            kmsvector means; // Cluster means

            double& operator[](const unsigned index) {
                return means[index];
            }
        public:
            typedef typename std::shared_ptr<clusters> ptr;

            clusters(const unsigned nclust, const unsigned ncol);
            clusters(const unsigned nclust, const unsigned ncol,
                    const kmsvector& means);

            static ptr create(const unsigned nclust, const unsigned ncol) {
                return ptr(new clusters(nclust, ncol));
            }

            static ptr create(const unsigned nclust, const unsigned ncol,
                    const kmsvector& means) {
                return ptr(new clusters(nclust, ncol, means));
            }

            const kmsvector& get_means() const {
                return means;
            }

            const int get_num_members(const unsigned idx) const {
                return num_members_v[idx];
            }

            const std::vector<int>& get_num_members_v() const {
                return num_members_v;
            }

            const bool is_complete(const unsigned idx) const {
                return complete_v[idx];
            }

            const unsigned size() const {
                return means.size();
            }

            void num_members_peq(const int val, const unsigned idx) {
                num_members_v[idx] += val;
            }

            double get(const unsigned index) {
                return means[index];
            }

            // Get an index (based on the entire chunck)
            const unsigned get_ncol() const {
                return ncol;
            }

            const unsigned get_nclust() const {
                return nclust;
            }

            const std::vector<bool>& get_complete_v() const {
                return complete_v;
            }

            clusters& operator+=(clusters& rhs);
            clusters& operator=(const clusters& other);
            bool operator==(const clusters& other);
            void peq(ptr rhs);
            const void print_means() const;
            void clear();
            /** \param idx the cluster index.
              */
            void set_mean(const kmsvector& mean, const int idx=-1);
            void set_mean(const double* mean, const int idx=-1);
            template<typename T>
                void set_mean(T& it, const int idx);
            void finalize(const unsigned idx);
            void unfinalize(const unsigned idx);
            template <typename T>
                void add_member(T& count_it, const unsigned idx);
            template <typename T>
                void add_member(const T* arr, const unsigned idx);
            template <typename T>
                void remove_member(T& count_it, const unsigned idx);
            template <typename T>
                void remove_member(const T* arr, const unsigned idx);
            template <typename T>
                void swap_membership(const T* arr, const unsigned from_idx,
                        const unsigned to_idx);
            template <typename T>
                void swap_membership(T& count_it,
                        const unsigned from_id, const unsigned to_id);
    };

    class prune_clusters : public clusters
    {
        private:
            kmsvector s_val_v;
            kmsvector prev_means;
            kmsvector prev_dist_v; // Distance to prev mean

            void init() {
                prev_means.resize(ncol*nclust);
                prev_dist_v.resize(nclust);
                s_val_v.assign(nclust, std::numeric_limits<double>::max());
            }

            prune_clusters(const unsigned nclust, const unsigned ncol):
                clusters(nclust, ncol) {
                    init();
            }

            prune_clusters(const unsigned nclust, const unsigned ncol,
                    const kmsvector means): clusters(nclust, ncol, means) {
                init();
            }

        public:
            typedef typename std::shared_ptr<prune_clusters> ptr;

            static ptr create(const unsigned nclust, const unsigned ncol) {
                return ptr(new prune_clusters(nclust, ncol));
            }

            static ptr create(const unsigned nclust, const unsigned ncol,
                    const kmsvector& mean) {
                return ptr(new prune_clusters(nclust, ncol, mean));
            }

            void set_s_val(const double val, const unsigned idx) {
                s_val_v[idx] = val;
            }

            double const get_s_val(const unsigned idx) { return s_val_v[idx]; }

            const kmsvector& get_prev_means() const {
                return prev_means;
            }

            void set_prev_means() {
                this->prev_means = means;
            }

            void set_prev_dist(const double dist, const unsigned idx) {
                prev_dist_v[idx] = dist;
            }

            double get_prev_dist(const unsigned idx) {
                return prev_dist_v[idx];
            }

            const void print_prev_means_v() const;
            void reset_s_val_v();
    };
}
#endif
