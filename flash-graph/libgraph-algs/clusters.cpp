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

#include <limits>
#include <vector>
#include <memory>
#include <algorithm>

#include "log.h"
#include "sem_kmeans_util.h"
#include "clusters.h"

using namespace km;

void clusters::clear() {
    std::fill(means.begin(), means.end(), 0);
    std::fill(num_members_v.begin(), num_members_v.end(), 0);
    std::fill(complete_v.begin(), complete_v.end(), false);
}

/** \param idx the cluster index.
*/
void clusters::set_mean(const kmsvector& mean, const int idx) {
    if (idx == -1) { // Set all means
        means = mean;
    } else {
        std::copy(mean.begin(), mean.end(),
                this->means.begin()+(idx*ncol));
    }
}

void clusters::set_mean(const double* mean, const int idx) {
    if (idx == -1) { // Set all means
        if (means.size() != (ncol*nclust))
            means.resize(ncol*nclust);
        std::copy(&(mean[0]), &(mean[ncol*nclust]), this->means.begin());
    } else {
        std::copy(&(mean[0]), &(mean[ncol]), this->means.begin()+(idx*ncol));
    }
}

template<typename T>
void clusters::set_mean(T& it, const int idx) {
    unsigned offset = idx*ncol;
    unsigned nid = 0;
    while(it.has_next()) {
        double e = it.next();
        means[offset+nid] = e;
    }
}

void clusters::finalize(const unsigned idx) {
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

void clusters::unfinalize(const unsigned idx) {
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
void clusters::add_member(T& count_it, const unsigned idx) {
    unsigned nid = 0;
    while(count_it.has_next()) {
        double e = count_it.next();
        means[(idx*ncol)+(nid++)] += e;
    }
    num_members_v[idx]++;
}

template <typename T>
void clusters::add_member(const T* arr, const unsigned idx) {
    unsigned offset = idx * ncol;
    for (unsigned i=0; i < ncol; i++) {
        means[offset+i] += arr[i];
    }
    num_members_v[idx]++;
}

template <typename T>
void clusters::remove_member(T& count_it, const unsigned idx) {
    unsigned nid = 0;
    while(count_it.has_next()) {
        double e = count_it.next();
        means[(idx*ncol)+nid++] -= e;
    }
    num_members_v[idx]--;
}

template <typename T>
void clusters::remove_member(const T* arr, const unsigned idx) {
    unsigned offset = idx * ncol;
    for (unsigned i=0; i < ncol; i++) {
        means[offset+i] -= arr[i];
    }
    num_members_v[idx]--;
}

template <typename T>
void clusters::swap_membership(const T* arr,
        const unsigned from_idx, const unsigned to_idx) {
    remove_member(arr, from_idx);
    add_member(arr, to_idx);
}

template <typename T>
void clusters::swap_membership(T& count_it,
        const unsigned from_id, const unsigned to_id) {
    unsigned nid = 0;
    unsigned from_offset = from_id * ncol;
    unsigned to_offset = to_id * ncol;
    while(count_it.has_next()) {
        double e = count_it.next();
        means[from_offset+nid] -= e;
        means[to_offset+nid++] += e;
    }
    num_members_v[from_id]--;
    num_members_v[to_id]++;
}

clusters& clusters::operator=(const clusters& other) {
    this->means = other.get_means();
    this->num_members_v = other.get_num_members_v();
    this->ncol = other.get_ncol();
    this->nclust = other.get_nclust();
    return *this;
}

bool clusters::operator==(const clusters& other) {
    return (get_ncol() == other.get_ncol() &&
            get_nclust() == other.get_nclust() &&
            v_eq(get_num_members_v(), other.get_num_members_v()) &&
            v_eq(get_means(), other.get_means())
           );
}

clusters& clusters::operator+=(clusters& rhs) {
    // TODO vectorize perhaps OR omp parallel
    for (unsigned i = 0; i < size(); i++)
        this->means[i] += rhs[i];

    for (unsigned idx = 0; idx < nclust; idx++)
        num_members_peq(rhs.get_num_members(idx), idx);
    return *this;
}

void clusters::peq(ptr rhs) {
    BOOST_VERIFY(rhs->size() == size());
    for (unsigned i = 0; i < size(); i++)
        this->means[i] += rhs->get(i);

    for (unsigned idx = 0; idx < nclust; idx++)
        num_members_peq(rhs->get_num_members(idx), idx);
}

// Begin Helpers //
const void clusters::print_means() const {
    for (unsigned cl_idx = 0; cl_idx < get_nclust(); cl_idx++) {
        std::cout << "#memb = " << get_num_members(cl_idx) << " ";
        print_arr<double>(&(means[cl_idx*ncol]), ncol);
    }
    std::cout << "\n";
}

clusters::clusters(const unsigned nclust, const unsigned ncol) {
    this->nclust = nclust;
    this->ncol = ncol;

    means.resize(ncol*nclust);
    num_members_v.resize(nclust);
    complete_v.assign(nclust, false);
}

clusters::clusters(const unsigned nclust, const unsigned ncol,
        const kmsvector& means) {
    this->nclust = nclust;
    this->ncol = ncol;

    set_mean(means);
    num_members_v.resize(nclust);
    complete_v.assign(nclust, true);
}

// Pruning clusters //
void prune_clusters::reset_s_val_v() {
    std::fill(s_val_v.begin(), s_val_v.end(),
            std::numeric_limits<double>::max());
}

const void prune_clusters::print_prev_means_v() const {
    for (unsigned cl_idx = 0; cl_idx < get_nclust(); cl_idx++) {
        print_arr<double>(&(prev_means[cl_idx*ncol]), ncol);
    }
    std::cout << "\n";
}
