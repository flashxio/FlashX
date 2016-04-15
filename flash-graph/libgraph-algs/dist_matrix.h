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

#include <memory>
#include <limits>

#include "log.h"

namespace km {
    class prune_clusters;
}

namespace prune {
    // NOTE: Creates a matrix like this e.g for K = 5
    /* - Don't store full matrix, don't store dist to myself -> space: (k*k-1)/2
       0 ==> 1 2 3 4
       1 ==> 2 3 4
       2 ==> 3 4
       3 ==> 4
       (4 ==> not needed)
       */
    class dist_matrix {
        private:
            std::vector<std::vector<double>> mat;
            unsigned rows;

            dist_matrix(const unsigned rows);
            void translate(unsigned& row, unsigned& col);

        public:
            typedef typename std::shared_ptr<dist_matrix> ptr;

            static ptr create(const unsigned rows) {
                return ptr(new dist_matrix(rows));
            }

            const unsigned get_num_rows() { return rows; }

            /* Do a translation from raw id's to indexes in the distance matrix */
            double get(unsigned row, unsigned col);
            // Testing purposes only
            double get_min_dist(const unsigned row);
            void set(unsigned row, unsigned col, double val);

            void print();
            void compute_dist(std::shared_ptr<km::prune_clusters> cls, const unsigned ncol);
    };
}
#endif
