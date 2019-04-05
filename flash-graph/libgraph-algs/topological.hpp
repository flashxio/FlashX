/*
 * Copyright 2017 Neurodata (https://neurodata.io)
 * Written by Disa Mhembere (disa@cs.jhu.edu)
 *
 * This file is part of Monya.
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

#ifndef FG_TOPOLOGICAL_HPP__
#define FG_TOPOLOGICAL_HPP__


#include <vector>
#include <cassert>
#include <algorithm>
#include <mutex>

#ifdef _OPENMP
#include <parallel/algorithm>
#endif

#include <utility>

namespace monya {
    typedef unsigned IndexType; // Tree node identifier type

    template <typename IndexType, typename ValType>
    class IndexVal {
        private:
            IndexType index;
            ValType val;
        public:
            IndexVal() {
            }

            IndexVal(const IndexType index, const ValType val) :
                index(index), val(val) {
            }

            void set(const IndexType index, const ValType val) {
                this->index = index;
                this->val = val;
            }

            void set_index(const IndexType index) {
                this->index = index;
            }

            const IndexType& get_index() const {
                return this->index;
            }

            void set_val(const ValType val) {
                this->val = val;
            }

            const ValType& get_val() const {
                return val;
            }

            bool operator< (const IndexVal& other) const {
                if (val == other.get_val())
                    return index < other.get_index();
                return val < other.get_val();
            }

            // Equivalence is only based on the value
            bool operator== (const IndexVal& other) const {
                return val == other.get_val();
            }

            bool operator!= (const IndexVal& other) const {
                return val != other.get_val();
            }

            std::string to_string() {
                return std::string("(") + std::to_string(get_index()) +
                    std::string(", ") + std::to_string(get_val()) +
                    std::string(")");
            }

            void print() {
               printf("%s ", to_string().c_str());
            }
    };

    template <typename T, typename IndexType, typename ValType>
    class IndexVector {
        private:
            std::vector<T> _;
            std::vector<std::mutex*> locks;
            std::vector<size_t> ranges;
            IndexType max_part_id;
            bool sorted;

            void add_lock() {
                locks.push_back(new std::mutex());
            }

        public:
            typedef typename std::vector<T>::iterator iterator;

            IndexVector() : sorted(false) { } // Default ctor

            std::vector<size_t>& get_ranges() {
                return ranges;
            }

            void lock(const size_t part_num) {
                locks[part_num]->lock();
            }

            void unlock(const size_t part_num) {
                locks[part_num]->unlock();
            }

            IndexVector(const std::vector<IndexType>& v) : IndexVector() {
                for (auto idx : v)
                    _.push_back(T(idx, 0)); // 0 is a place holder
            }

            void swap(IndexType arg0, IndexType arg1, size_t start_idx,
                    size_t end_offset) {

                // FIXME This is the simple case when we do linear search
                bool found0, found1;
                size_t pos0, pos1;
                found0 = found1 = false;
                for (size_t idx = start_idx; idx < start_idx+end_offset;
                        idx++) {
                    if (!found0 && _[idx].get_index() == arg0) {
                        pos0 = idx; found0 = true;
                    } else if (!found1 && _[idx].get_index() == arg1) {
                        pos1 = idx; found1 = true;
                    }
                    if (found0 && found1) break;
                }
                assert(found0 && found1);
                std::swap(_[pos0], _[pos1]);
            }

            void resize(const size_t nelem) {
                _.resize(nelem);
            }

            T& operator[](const int index) {
                return this->_[index];
            }

            void print(const size_t max_psize=100) {
                auto psize = _.size() < max_psize ? _.size() : max_psize;

                printf("[ ");
                for (size_t i = 0; i < psize; i++)
                    _[i].print();
                if (max_psize == psize && _.size() != max_psize)
                    printf("... ");
                printf("]\n");
            }

            std::string to_string() {
                std::string __repr__ = "";
                for (auto item : _)
                    __repr__ + item.to_string() + std::string("\n");
                return __repr__;
            }

            // Insert values with contiguous indexes
            void set_indexes(const ValType* vals, const size_t nelem) {
                if (_.size())
                    _.clear();

                for (size_t i = 0; i < nelem; i++)
                    _.push_back(T(i, vals[i]));
            }

            void get_indexes(std::vector<IndexType>& v) {
                assert(v.empty());

                for (auto i = begin(); i != end(); ++i) {
                    v.push_back(i->get_index());
                }
            }

            void append(const IndexType id, const ValType val) {
                _.push_back(T(id, val));
            }

            void append(T& iv) {
                _.push_back(iv);
            }

            const size_t size() const {
                return _.size();
            }

            void insert(T item, const size_t offset) {
                _.insert(begin()+offset, item);
            }

            void trim(const size_t nelem) {
                _.resize(nelem);
            }

            bool empty() const { return _.empty(); }
            iterator begin() { return _.begin(); }
            iterator end() { return _.end(); }
            iterator rbegin() { return _.rbegin(); }
            iterator rend() { return _.rend(); }

            T* find(T& iv) {
                if (!is_sorted())
                    sort();
                return std::binary_search(begin(), end(), iv);
            }

            T* find(T& iv, T* _start, T* _end) {
                if (!is_sorted())
                    sort();
                return std::binary_search(_start, _end, iv);
            }

            struct IndexComp {
                inline bool operator() (T& arg1, T& arg2) {
                    return (arg1.get_index() < arg2.get_index());
                }
            };

            void sort() {
#ifdef _OPENMP
                omp_set_num_threads(omp_get_max_threads());
                __gnu_parallel::sort(begin(), end());
#else
                std::sort(begin(), end());
#endif
                sorted = true;
            }

            void set_sorted(const bool sorted) {
                this->sorted = sorted;
            }

            const bool is_sorted() const {
                return sorted;
            }

            void clear() {
                sorted = false;
                _.clear();
            }

            const size_t locks_size() {
                return locks.size();
            }

            // Obtain the partition mapping
            void partition(std::vector<IndexType>& part_idx) {
                if (_.size()) {
                    ranges.push_back(0);
                    add_lock();
                } else {
                    return;
                }

                IndexType part = 0;

                auto prev_val = _[0].get_val();
                part_idx[_[0].get_index()] = part;

                for (size_t i = 1; i < _.size(); i++) {
                    if (_[i].get_val() != prev_val) {
                        part++;

                        ranges.push_back(i);
                        add_lock();
                        prev_val = _[i].get_val();
                    }
                    part_idx[_[i].get_index()] = part;
                }
                max_part_id = part;
            }

            const IndexType& get_max_part_id() {
                return max_part_id;
            }

            ~IndexVector() {
                for (auto& lock : locks)
                    delete lock;
            }
    };
} // End monya

#endif
