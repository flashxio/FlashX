/*
 * Copyright 2016 Open Connectome Project (http://openconnecto.me)
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

#ifndef __PRUNE_STATS_H__
#define __PRUNE_STATS_H__

#include <vector>
#include <map>

namespace {
    // Class to hold stats on the effectiveness of pruning
    class prune_stats
    {
        private:
            // Counts per iteration
            unsigned lemma1, _3a, _3b, _3c, _4;

            // Total counts
            unsigned long tot_lemma1, tot_3a, tot_3b, tot_3c, tot_4, iter;
            unsigned nrow;
            unsigned nclust;

            prune_stats(const unsigned nrows, const unsigned nclust) {
                _3a = 0; _3b = 0; lemma1 = 0; _3c = 0; _4 = 0;
                tot_lemma1 = 0; tot_3a = 0; tot_3b = 0;
                tot_3c = 0; tot_4 = 0; iter = 0;

                this->nrow = nrows;
                this->nclust = nclust;
            }

        public:
            typedef std::shared_ptr<prune_stats> ptr;

            static ptr create(const unsigned nrows, const unsigned nclust) {
                return ptr(new prune_stats(nrows, nclust));
            }
            void pp_lemma1(const unsigned var=1) { lemma1 += var; }
            void pp_3a() { _3a++; }
            void pp_3b() { _3b++; }
            void pp_3c() { _3c++; }
            void pp_4() { _4++; }

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
                BOOST_VERIFY((lemma1 + _3a + _3b + _3c + _4) <=  nrow*nclust);
                BOOST_LOG_TRIVIAL(info) << "\n\nPrune stats count:\n"
                    "lemma1 = " << lemma1 << ", 3a = " << _3a
                    << ", 3b = " << _3b << ", 3c = " << _3c << ", 4 = " << _4;

                BOOST_LOG_TRIVIAL(info) << "\n\nPrune stats \%s:\n"
                    "lemma1 = " << (lemma1 == 0 ? 0 : ((double)lemma1/(nrow*nclust))*100) <<
                    "\%, 3a = " << (_3a == 0 ? 0 : ((double)_3a/(nrow*nclust))*100) <<
                    "\%, 3b = " << (_3b == 0 ? 0 : ((double) _3b/(nrow*nclust))*100) <<
                    "\%, 3c = " << (_3c == 0 ? 0 : ((double) _3c/(nrow*nclust))*100) <<
                    "\%, 4 = " << (_4 == 0 ? 0 : ((double) _4/(nrow*nclust))*100) << "\%";

                tot_lemma1 += lemma1;
                tot_3a += _3a;
                tot_3b += _3b;
                tot_3c += _3c;
                tot_4 += _4;

                lemma1 = 0; _3a = 0; _3b = 0; _3c = 0; _4 = 0; // reset
            }

            std::vector<double> get_stats() {
                double perc_lemma1 = (tot_lemma1 / ((double)(nrow*iter*nclust)))*100;
                double perc_3a = (tot_3a / ((double)(nrow*iter*nclust)))*100;
                double perc_3b = (tot_3b / ((double)(nrow*iter*nclust)))*100;
                double perc_3c = (tot_3c / ((double)(nrow*iter*nclust)))*100;
                double perc_4 = (tot_4 / ((double)(nrow*iter*nclust)))*100;
                // Total percentage
                double perc = ((tot_3b + tot_3a + tot_3c + tot_4 + tot_lemma1) /
                        ((double)(nrow*iter*nclust)))*100;

                printf("tot_lemma1 = %lu, tot_3a = %lu, tot_3b = %lu, tot_3c = %lu, tot_4 = %lu\n",
                        tot_lemma1, tot_3a, tot_3b, tot_3c, tot_4);

                BOOST_LOG_TRIVIAL(info) << "\n\nPrune stats total:\n"
                    "Tot = " << perc << "\%, 3a = " << perc_3a <<
                    "\%, 3b = " << perc_3b << "\%, 3c = " << perc_3c
                    << "\%, 4 = " << perc_4 << "\%, lemma1 = " << perc_lemma1 << "\%";

                std::vector<double> ret {perc_lemma1, perc_3a, perc_3b, perc_3c, perc_4, perc};
                return ret;
            }
    };

    class active_counter {
        private:
        std::vector<bool> prev_active; // Was a vertex active last iter
        // Active in this iter & active in last
        std::vector<std::vector<bool>> active;
        size_t nrow;

        active_counter(const size_t nrow) {
            this->nrow = nrow;
            prev_active.assign(nrow, false);
            init_iter();
        }

        // Was active in the prev iteration
        const bool was_active(const size_t row) const {
            return prev_active[row];
        }

        void consolidate_samples(std::map<std::vector<bool>, size_t>& data_hash,
                const size_t rows) {
            for (size_t row = 0; row < rows; row++) {
                std::vector<bool> sample(active.size());
                for (unsigned iter = 0; iter < active.size(); iter++) {
                    sample[iter] = active[iter][row];
                }

                std::map<std::vector<bool>, size_t>::iterator it = data_hash.find(sample);
                if (it != data_hash.end())
                    it->second++;
                else
                    data_hash[sample] = 1;
            }
        }

        public:
        typedef std::shared_ptr<active_counter> ptr;
        static ptr create(const size_t nrow) {
            return ptr(new active_counter(nrow));
        }

        void init_iter() {
            std::vector<bool> v;
            v.assign(nrow, false);
            active.push_back(v); // iteration i all are initially false
        }

        void is_active(const size_t row, const bool val) {
            if (active.size() == 1 && was_active(row))
                BOOST_ASSERT_MSG(false, "In first iter the row cannot be active previously");

            if (val && was_active(row)) {
                // 1. Grow the rows active vec, to add a true
                active.back()[row] = true;
            } else {
                active.back()[row] = false;
            }

            prev_active[row] = val; // Seen in next iteration
        }

        void write_raw(std::string fn, size_t print_row_cnt) {

            if (print_row_cnt > nrow)
                print_row_cnt = nrow;

            std::string out = "";
            for (size_t row = 0; row < print_row_cnt; row++) {
                for (unsigned iter = 0; iter < active.size(); iter++) {
                    if (iter == 0)
                        out += std::to_string(row) + ", ";

                    if (iter + 1 == active.size())
                        out += std::to_string(active[iter][row]) + "\n";
                    else
                        out += std::to_string(active[iter][row]) + ", ";
                }
            }

            FILE* f = fopen(fn.c_str(), "wb");
            fwrite((char*)out.c_str(), out.size(), 1, f);
            fclose(f);

            //printf("Printing the activation:\n%s\n", out.c_str());
        }

        void write_consolidated(std::string fn, size_t print_row_cnt) {
            if (print_row_cnt > nrow)
                print_row_cnt = nrow;

            std::string out = "";
            std::map<std::vector<bool>, size_t> dhash;
            consolidate_samples(dhash, print_row_cnt);

            // Iterate and print
            for (std::map<std::vector<bool>, size_t>::iterator it = dhash.begin();
                    it != dhash.end(); ++it) {

                out += std::to_string(it->second) + ", ";
                std::vector<bool> v = it->first;
                for (std::vector<bool>::iterator vit = v.begin();
                        vit != v.end(); ++vit) {
                    if (vit + 1 == v.end())
                        out += std::to_string(*vit) + "\n";
                    else
                        out += std::to_string(*vit) + ", ";
                }
            }

            FILE* f = fopen(fn.c_str(), "wb");
            fwrite((char*)out.c_str(), out.size(), 1, f);
            fclose(f);

            //printf("Consolidated activation:\n%s\n", out.c_str());
        }
    };
}
#endif
