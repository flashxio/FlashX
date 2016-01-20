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

#ifndef __DIST_MATRIX_H__
#define __DIST_MATRIX_H__

#include "cluster.h"
#include "sem_kmeans_util.h"

namespace
{
    // NOTE: Creates a matrix like this e.g for K = 5
    /* - Don't store full matrix, don't store dist to myself -> space: (k*k-1)/2
       0 ==> 1 2 3 4
       1 ==> 2 3 4
       2 ==> 3 4
       3 ==> 4
       (4 ==> not needed)
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

            void translate(unsigned& row, unsigned& col) {
                // First make sure the smaller is the row
                if (row > col) {
                    std::swap(row, col);
                }

                BOOST_VERIFY(row < rows);
                col = col - row - 1; // Translation
                BOOST_VERIFY(col < (rows - row));
            }

        public:
            typedef typename std::shared_ptr<dist_matrix> ptr;

            static ptr create(const unsigned rows) {
                return ptr(new dist_matrix(rows));
            }

            /* Do a translation from raw id's to indexes in the distance matrix */
            double get(unsigned row, unsigned col) {
                if (row == col) { return std::numeric_limits<double>::max(); }
                translate(row, col);
                return mat[row][col];
            }

            // Testing purposes only
            double get_min_dist(const unsigned row) {
                double best = std::numeric_limits<double>::max();
                for (unsigned col = 0; col < rows+1; col++) {
                    if (col != row) {
                        double val = get(row, col);
                        if (val < best) best = val;
                    }
                }
                BOOST_VERIFY(best < std::numeric_limits<double>::max());
                return best;
            }

            void set(unsigned row, unsigned col, double val) {
                BOOST_VERIFY(row != col);
                translate(row, col);
                mat[row][col] = val;
            }

            const unsigned get_num_rows() { return rows; }

            void print() {
                for (unsigned row = 0; row < rows; row++) {
                    std::cout << row << " ==> ";
                    print_vector<double>(mat[row]);
                }
            }

            void compute_dist(std::vector<prune_cluster::ptr>& vcl, const unsigned num_clust) {
                if (num_clust <= 1) return;

                BOOST_VERIFY(get_num_rows() == vcl.size()-1); // -1 since the last item has no row

                for (unsigned i = 0; i < num_clust; i++) {
                    vcl[i]->reset_s_val();
                }

                //#pragma omp parallel for collapse(2) // FIXME: Opt Coalese perhaps
                for (unsigned i = 0; i < num_clust; i++) {
                    for (unsigned j = i+1; j < num_clust; j++) {
                        double dist = eucl_dist<std::vector<double>>(&(vcl[i]->get_mean()),
                                &(vcl[j]->get_mean())) / 2.0;
                        set(i,j, dist);

                        // Set s(x) for each cluster
                        if (dist < vcl[i]->get_s_val()) {
                            vcl[i]->set_s_val(dist);
                        }

                        if (dist < vcl[j]->get_s_val()) {
                            vcl[j]->set_s_val(dist);
                        }
                    }
                }
#if KM_TEST
                for (unsigned cl = 0; cl < num_clust; cl++)
                    BOOST_VERIFY(vcl[cl]->get_s_val() == get_min_dist(cl));
#endif
            }
    };
}
#endif
