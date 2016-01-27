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

#ifndef __CLUSTER_H__
#define __CLUSTER_H__

#include <vector>
#include <memory>
#include <boost/assert.hpp>

namespace {
    class cluster
    {
        friend class prune_cluster;
        private:
            std::vector<double> mean; // Cluster mean
            int num_members; // Cluster assignment count
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
                mean.resize(len);
                num_members = 0;
                complete = false;
            }

            cluster(const std::vector<double> mean) {
                set_mean(mean);
                num_members = 0;
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

            /*
            void init(const unsigned len) {
                mean.resize(len);
                num_members = 0;
            }*/

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

            void num_members_peq(const int val) {
                num_members += val;
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
                        " UNfinalized object";
                    return;
                }
                complete = false;

                for (unsigned i = 0; i < size(); i++) {
                    this->mean[i] *= (double)num_members;
                }
            }

            template <typename T>
                void add_member(T& count_it) {
                    unsigned nid = 0;
                    while(count_it.has_next()) {
                        double e = count_it.next();
                        mean[nid++] += e;
                    }
                    num_members++;
                }

            cluster& operator=(const cluster& other) {
                this->mean = other.get_mean();
                this->num_members = other.get_num_members();
                return *this;
            }

            double& operator[](const unsigned index) {
                BOOST_VERIFY(index < mean.size());
                return mean[index];
            }

            cluster& operator+=(cluster& rhs) {
                BOOST_VERIFY(rhs.size() == size());
                // TODO vectorize perhaps
                for (unsigned i = 0; i < mean.size(); i++) {
                    this->mean[i] += rhs[i];
                }
                this->num_members += rhs.get_num_members();
                return *this;
            }
    };

    class prune_cluster : public cluster
    {
        private:
            double s_val;
            std::vector<double> prev_mean; // For lwr_bnd
            double prev_dist; // Distance to prev mean

            prune_cluster(const unsigned len): cluster(len) {
                prev_mean.resize(len);
            }

            prune_cluster(const std::vector<double> mean): cluster(mean) {
                prev_mean.resize(mean.size());
                reset_s_val();
            }

        public:
            void reset_s_val() {
                s_val = std::numeric_limits<double>::max();
            }

            void set_s_val(double val) {
                s_val = val;
            }

            double const get_s_val() { return s_val; }

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

            typedef typename std::shared_ptr<prune_cluster> ptr;

            static ptr create(const unsigned len) {
                return ptr(new prune_cluster(len));
            }

            static ptr create(const std::vector<double>& mean) {
                return ptr(new prune_cluster(mean));
            }

            template <typename T>
                void remove_member(T& count_it) {
                    unsigned nid = 0;
                    while(count_it.has_next()) {
                        double e = count_it.next();
                        mean[nid++] -= e;
                    }
                    num_members--;
                }
    };
}
#endif
