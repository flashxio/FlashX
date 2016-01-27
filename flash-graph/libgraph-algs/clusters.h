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

namespace {

    typedef std::vector<double> kmsvector;
    typedef std::vector<double>::iterator kmsiterator;

    /**
     * \brief Print an arry of some length `len`.
     *   \param len The length of the array.
     */
    template <typename T>
        static void print_arr(const T* arr, const unsigned len) {
            printf("[ ");
            for (unsigned i = 0; i < len; i++) {
                std::cout << arr[i] << " ";
            }
            printf("]\n");
        }

    // Vector equal function
    template <typename T>
        static bool v_eq(const T& lhs, const T& rhs) {
            return std::equal(lhs.begin(), lhs.end(), rhs.begin());
        }


    template <typename T>
        static const bool v_eq_const(const std::vector<T>& v, const T var) {
            for (unsigned i=0; i < v.size(); i++) {
                if (v[i] != var) return false;
            }
            return true;
        }

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

            clusters(const unsigned nclust, const unsigned ncol) {
                this->nclust = nclust;
                this->ncol = ncol;

                means.resize(ncol*nclust);
                num_members_v.resize(nclust);
                complete_v.assign(nclust, false);
            }

            clusters(const unsigned nclust, const unsigned ncol,
                    const kmsvector& means) {
                this->nclust = nclust;
                this->ncol = ncol;

                set_mean(means);
                num_members_v.resize(nclust);
                complete_v.assign(nclust, true);
            }

        public:
            typedef typename std::shared_ptr<clusters> ptr;

            static ptr create(const unsigned nclust, const unsigned ncol) {
                return ptr(new clusters(nclust, ncol));
            }

            static ptr create(const unsigned nclust, const unsigned ncol,
                    const kmsvector& means) {
                return ptr(new clusters(nclust, ncol, means));
            }

            void clear() {
                std::fill(means.begin(), means.end(), 0);
                std::fill(num_members_v.begin(), num_members_v.end(), 0);
                std::fill(complete_v.begin(), complete_v.end(), false);
            }

            const kmsvector& get_means() const {
                return means;
            }

            void set_mean(const kmsvector& mean, const int idx=-1) {
                if (idx == -1) { // Set all means
                    means = mean;
                } else {
                    std::copy(mean.begin(), mean.end(), this->means.begin()+(idx*ncol));
                }
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

            void finalize(const unsigned idx) {
                if (is_complete(idx)) {
                    BOOST_LOG_TRIVIAL(info) << "WARNING: Calling finalize() on a"
                        " finalized object";
                    return;
                }

                if (num_members_v[idx] > 1) { // Less than 2 is the same result
                    for (unsigned i = 0; i < ncol; i++) {
                        means[(idx*ncol)+i] /= double(num_members_v[idx]);
                    }
                }
                complete_v[idx] = true;
            }

            void unfinalize(const unsigned idx) {
                if (!is_complete(idx)) {
                    BOOST_LOG_TRIVIAL(info) << "WARNING: Calling unfinalize() on an"
                        " UNfinalized object";
                    return;
                }
                complete_v[idx] = false;

                for (unsigned col = 0; col < ncol; col++) {
                    this->means[(ncol*idx) + col] *= (double)num_members_v[idx];
                }
            }

            template <typename T>
                void add_member(T& count_it, const unsigned idx) {
                    unsigned nid = 0;
                    while(count_it.has_next()) {
                        double e = count_it.next();
                        means[(idx*ncol)+(nid++)] += e;
                    }
                    num_members_v[idx]++;
                }

            template <typename T>
                void add_member(const T* arr, const unsigned idx) {
                    for (unsigned i=0; i < ncol; i++) {
                        means[(idx*ncol)+i] += arr[i];
                    }
                    num_members_v[idx]++;
                }

            template <typename T>
                void remove_member(T& count_it, const unsigned idx) {
                    unsigned nid = 0;
                    while(count_it.has_next()) {
                        double e = count_it.next();
                        means[(idx*ncol)+nid++] -= e;
                    }
                    num_members_v[idx]--;
                }

            template <typename T>
                void remove_member(const T* arr, const unsigned idx) {
                    for (unsigned i=0; i < ncol; i++) {
                        means[(idx*ncol)+i] -= arr[i];
                    }
                    num_members_v[idx]--;
                }

            clusters& operator=(const clusters& other) {
                this->means = other.get_means();
                this->num_members_v = other.get_num_members_v();
                this->ncol = other.get_ncol();
                this->nclust = other.get_nclust();
                return *this;
            }

            bool operator==(const clusters& other) {
                return (get_ncol() == other.get_ncol() &&
                        get_nclust() == other.get_nclust() &&
                        v_eq(get_num_members_v(), other.get_num_members_v()) &&
                        v_eq(get_means(), other.get_means())
                );
            }

            double& operator[](const unsigned index) {
                BOOST_VERIFY(index < size()); // TODO: rm -- called to often!
                return means[index];
            }

            clusters& operator+=(clusters& rhs) {
                // TODO: rm -- called to often!
                BOOST_VERIFY(rhs.size() == size() &&
                        (v_eq(get_complete_v(), rhs.get_complete_v())) &&
                        v_eq_const(get_complete_v(), false));

                // TODO vectorize perhaps
                for (unsigned i = 0; i < size(); i++)
                    this->means[i] += rhs[i];

                for (unsigned idx = 0; idx < nclust; idx++)
                    num_members_peq(rhs.get_num_members(idx), idx);

                return *this;
            }

            const unsigned get_ncol() const {
                return ncol;
            }

            const unsigned get_nclust() const {
                return nclust;
            }
            
            const std::vector<bool>& get_complete_v() const {
                return complete_v;
            }

            // Begin Helpers //
            const void print(const kmsvector& v) const {
                for (unsigned cl_idx = 0; cl_idx < get_nclust(); cl_idx++) {
                    std::cout << "#memb = " << get_num_members(cl_idx) << " ";
                        print_arr<double>(&(v[cl_idx*ncol]), ncol);
                }
                std::cout << "\n";
            }
    };

    class prune_clusters : public clusters
    {
        private:
            kmsvector s_val_v;
            kmsvector prev_means;
            kmsvector prev_dist_v; // Distance to prev mean

            prune_clusters(const unsigned nclust, const unsigned ncol):
                clusters(nclust, ncol) {
                prev_means.resize(ncol*nclust);
                prev_dist_v.resize(nclust);
            }

            prune_clusters(const unsigned nclust, const unsigned ncol,
                    const kmsvector means): clusters(nclust, ncol, means) {
                prev_means.resize(ncol*nclust);
                s_val_v.assign(nclust, std::numeric_limits<double>::max());
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

            static ptr cast2(clusters::ptr cl) {
                return std::static_pointer_cast<prune_clusters, clusters>(cl);
            }

            void reset_s_val_v() {
                std::fill(s_val_v.begin(), s_val_v.end(),
                        std::numeric_limits<double>::max());
            }

            void set_s_val(const double val, const unsigned idx) {
                s_val_v[idx] = val;
            }

            double const get_s_val(const unsigned idx) { return s_val_v[idx]; }

            const kmsvector& get_prev_means() const {
                return prev_means;
            }

            // TODO: Test vs double* &(prev_means[idx*ncol])
            // TODO: Return reference
            const kmsiterator get_prev_mean(const unsigned idx) {
                return prev_means.begin()+(ncol*idx);
            }

            void set_prev_mean(const unsigned idx) {
                if (!is_complete(idx)) {
                    BOOST_LOG_TRIVIAL(warning) << "WARNING: Doing nothing for "
                        "unfinalized mean. Permissible once";
                    return;
                }

                // TODO: Check me
                std::copy(means.begin()+(idx*ncol), means.end()+((idx+1)*ncol),
                        this->prev_means.begin()+(idx*ncol));
            }

            void set_prev_dist_v(const double dist, const unsigned idx) {
                prev_dist_v[idx] = dist;
            }

            double get_prev_dist(const unsigned idx) {
                return prev_dist_v[idx];
            }
    };
}
#endif
