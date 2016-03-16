#ifndef __THD_SAFE_BOOL_VECTOR_H__
#define __THD_SAFE_BOOL_VECTOR_H__

/*
 * Copyright 2014 Open Connectome Project (http://openconnecto.me)
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

#include <string.h>
#include <vector>
#include <limits>
#include <memory>
#include <iostream>
#include <boost/assert.hpp>
#include "log.h"
#include "../libgraph-algs/sem_kmeans_util.h"

namespace {
    constexpr unsigned LEN = 2;

    class _bool{
        private:
            char _ [LEN];
        public:
            _bool() { }
            _bool(char c) {
                _[0] = c;
            }

            friend std::ostream& operator<<(std::ostream& os,
                    const _bool& b);

            operator bool() const {
                if (strcmp(_, "1"))
                    return false;
                return true;
            }
    };

    std::ostream& operator<<(std::ostream& os,
            const _bool& b) {
        os << b._;
        return os;
    }

    /**
      * \brief This class uses a char[2] vector to represent a
      *     boolean vector since it's represented as a bit vector
      *     which leads to data races for concurrent writers.
    */
    class thd_safe_bool_vector {
        private:
            std::vector<_bool> data;

            thd_safe_bool_vector(const size_t len) {
                data.resize(len); 
            }

            thd_safe_bool_vector(const size_t len, const bool init) : 
                thd_safe_bool_vector(len) {
                    if (init) {
                        for (unsigned i = 0; i < data.size(); i++)
                            data[i] = _bool('1');
                    } else {
                        for (unsigned i = 0; i < data.size(); i++)
                            data[i] = _bool('0');
                    }
                }
        public:
            typedef std::shared_ptr<thd_safe_bool_vector> ptr;
            static ptr create(const size_t len) {
                return ptr(new thd_safe_bool_vector(len));
            }

            static ptr create(const size_t len, const bool init) {
                return ptr(new thd_safe_bool_vector(len, init));
            }

            const bool get(const unsigned idx) const {
                return data[idx];
            }

            void set(const unsigned idx, const bool val) {
                if (val) {
                    data[idx] = _bool('1');
                } else {
                    data[idx] = _bool('0');
                }
            }

            unsigned size() const {
                return data.size();
            }

            void print() const {
                print_vector<_bool>(data);
            }
    };
}
#endif
