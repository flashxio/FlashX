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
#ifndef __SEM_KMEANS_UTIL_H__
#define __SEM_KMEANS_UTIL_H__

#include <stdlib.h>
#include <stdio.h>

#include <vector>
#include <iostream>
#include <boost/assert.hpp>

using namespace fg;

namespace {

    template <typename T>
        static void print_vector(typename std::vector<T> v, unsigned max_print=100) {
            unsigned print_len = v.size() > max_print ? max_print : v.size();

            std::cout << "[";
            typename std::vector<T>::iterator itr = v.begin();
            for (; itr != v.begin()+print_len; itr++) {
                std::cout << " "<< *itr;
            }

            if (v.size() > print_len) std::cout << " ...";
            std::cout <<  " ]\n";
        }

    // A very C-style binary data reader
    template <typename T>
        class bin_reader {
            private:
                FILE* f;
                unsigned nrow, ncol;

                void cat(const T* arr) {
                    std::cout << "[ ";
                    for (unsigned i = 0; i < ncol; i++) {
                        std::cout << arr[i] << " ";
                    }
                    std::cout << "]\n";
                }

            public:
                bin_reader(std::string fn, unsigned nrow, unsigned ncol) {
                    f = fopen(fn.c_str(), "rb");
                    BOOST_VERIFY(NULL != f);
                    this->nrow = nrow;
                    this->ncol = ncol;
                }

                // Read data and cat in a viewer friendly fashion
                void read_cat() {
                    T arr [ncol];
                    for (unsigned i = 0; i < nrow; i++) {
                        size_t num_read = fread(&arr[0], sizeof(T)*ncol, 1, f);
                        BOOST_ASSERT_MSG(num_read == 1, "Error reading file!\n");
                        cat(arr);
                    }
                }

                std::vector<double> readline() {
                    std::vector<double> v;
                    v.resize(ncol);
                    size_t num_read = fread(&v[0], sizeof(v[0])*ncol, 1, f);
                    BOOST_ASSERT_MSG(num_read == 1, "Error reading file!\n");
                    return v;
                }

                void readline(T* v) {
                    size_t num_read = fread(&v[0], sizeof(v[0])*ncol, 1, f);
                    BOOST_ASSERT_MSG(num_read == 1, "Error reading file!\n");
                }

                ~bin_reader() {
                    fclose(f);
                }
        };

    // Class to hold stats on the effectiveness of pruning
    class prune_stats
    {
        private:
            // Counts per iteration
            unsigned lemma1, _3a, _3b, _3c, _4;

            // Total counts
            unsigned tot_lemma1, tot_3a, tot_3b, tot_3c, tot_4, iter;
            unsigned nrow;
            unsigned K;

            prune_stats(const unsigned nrows, const unsigned nclust) {
                _3a = _3b = lemma1 = _3c = _4 = 0;
                tot_lemma1 = tot_3a = tot_3b = tot_3c = tot_4 = iter = 0;

                this->nrow = nrows;
                this->K = nclust;
            }

        public:
            typedef std::shared_ptr<prune_stats> ptr;

            static ptr create(const unsigned nrows, const unsigned nclust) {
                return ptr(new prune_stats(nrows, nclust));
            }
            void pp_lemma1(unsigned var=1) {
                lemma1 = lemma1 + var;
            }
            void pp_3a() {
                _3a = _3a + 1;
            }
            void pp_3b() {
                _3b = _3b + 1;
            }
            void pp_3c() {
                _3c = _3c + 1;
            }
            void pp_4() {
                _4 = _4 + 1;
            }

            const unsigned get_lemma1() const { return lemma1; }
            const unsigned get_3a() const { return _3a; }
            const unsigned get_3b() const { return _3b; }
            const unsigned get_3c() const { return _3c; }
            const unsigned get_4() const { return _4; }

            prune_stats& operator+=(prune_stats& other) {
                lemma1 += other.get_lemma1();
                _3a += other.get_3a();
                _3b += other.get_3b();
                _3c += other.get_3c();
                _4 += other.get_4();
                return *this;
            }

            void finalize() {
                iter++;
                BOOST_VERIFY((lemma1 + _3a + _3b + _3c + _4) <=  nrow*K);
                BOOST_LOG_TRIVIAL(info) << "\n\nPrune stats count:\n"
                    "lemma1 = " << lemma1 << ", 3a = " << _3a
                    << ", 3b = " << _3b << ", 3c = " << _3c << ", 4 = " << _4;

                BOOST_LOG_TRIVIAL(info) << "\n\nPrune stats \%s:\n"
                    "lemma1 = " << (lemma1 == 0 ? 0 : ((double)lemma1/(nrow*K))*100) <<
                    "\%, 3a = " << (_3a == 0 ? 0 : ((double)_3a/(nrow*K))*100) <<
                    "\%, 3b = " << (_3b == 0 ? 0 : ((double) _3b/(nrow*K))*100) <<
                    "\%, 3c = " << (_3c == 0 ? 0 : ((double) _3c/(nrow*K))*100) <<
                    "\%, 4 = " << (_4 == 0 ? 0 : ((double) _4/(nrow*K))*100) << "\%";

                tot_lemma1 += lemma1;
                tot_3a += _3a;
                tot_3b += _3b;
                tot_3c += _3c;
                tot_4 += _4;

                lemma1 = _3a = _3b = _3c = _4 = 0; // reset
            }

            std::vector<double> get_stats() {
                double perc_lemma1 = tot_lemma1 / ((double)(nrow*this->iter*K))*100;
                double perc_3a = tot_3a / ((double)(nrow*this->iter*K))*100;
                double perc_3b = tot_3b / ((double)(nrow*this->iter*K))*100;
                double perc_3c = tot_3c / ((double)(nrow*this->iter*K))*100;
                double perc_4 = tot_4 / ((double)(nrow*this->iter*K))*100;
                // Total percentage
                double perc = ((tot_3b + tot_3a + tot_3c + tot_4 + tot_lemma1) /
                        ((double)((nrow*this->iter*K)) + ((K*(K-1))/2.0)))*100;
                BOOST_LOG_TRIVIAL(info) << "\n\nPrune stats total:\n"
                    "Tot = " << perc << "\%, 3a = " << perc_3a <<
                    "\%, 3b = " << perc_3b << "\%, 3c = " << perc_3c
                    << "\%, 4 = " << perc_4 << "\%, lemma1 = " << perc_lemma1 << "\%";

                std::vector<double> ret {perc_lemma1, perc_3a, perc_3b, perc_3c, perc_4, perc};
                return ret;
            }
    };
}
#endif
