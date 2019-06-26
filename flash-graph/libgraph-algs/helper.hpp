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

#ifndef FG_LIBGRAPH_HELPER_HPP__
#define FG_LIBGRAPH_HELPER_HPP__

#include <map>
#include <unordered_map>
#include <vector>
#include <atomic>

namespace fg {
// Unordered Map
template <typename K, typename V>
void print(const std::unordered_map<K,V>& map) {
#ifndef BIND
    for (auto const& kv : map) {
        std::cout << "k: " << kv.first << ", v: " << kv.second << std::endl;
    }
    std::cout << "\n";
#endif
}

// Map
template <typename K, typename V>
void print(const std::map<K,V>& map) {
#ifndef BIND
    for (auto const& kv : map) {
        std::cout << "k: " << kv.first << ", v: " << kv.second << std::endl;
    }
    std::cout << "\n";
#endif
}

// Vector
template <typename T>
void print(const typename std::vector<T>& v, size_t max_print=100) {
#ifndef BIND
    auto print_len = v.size() > max_print ? max_print : v.size();
    std::cout << "[";
    for (auto const& val : v)
        std::cout << " "<< val;

    if (v.size() > print_len) std::cout << " ...";
    std::cout <<  " ]\n";
#endif
}

// Array
template <typename T>
void print(const T* arr, const unsigned len) {
#ifndef BIND
    printf("[ ");
    for (unsigned i = 0; i < len; i++) {
        std::cout << arr[i] << " ";
    }
    printf("]\n");
#endif
}

class barrier {
    private:
        std::atomic<unsigned> ncomplete;
        unsigned nmembers;
    public:
        typedef std::shared_ptr<barrier> ptr;

        barrier(const unsigned nmembers) {
            ncomplete.store(0);
            set_nmembers(nmembers);
        }

        static ptr create(const unsigned nmembers) {
            return ptr(new barrier(nmembers));
        }

        void set_nmembers(const unsigned nmembers) {
            this->nmembers = nmembers;
        }

        const unsigned get_nmembers() const {
            return nmembers;
        }

        bool ping() {
            ncomplete++;
            bool complete = nmembers == ncomplete.load();
            if (complete) {
                ncomplete.store(0); // Always reset
            }
            return complete;
        }
};

template <typename T>
class atomicwrapper
{
    private:
        std::atomic<T> _a;

    public:
        atomicwrapper() :_a() {}
        atomicwrapper(const std::atomic<T> &a) :_a(a.load()) {}
        atomicwrapper(const atomicwrapper &other) :_a(other._a.load()) {}
        atomicwrapper &operator=(const atomicwrapper &other) {
            _a.store(other._a.load());
            return *this;
        }

        T get() {
            return _a;
        }

        void assign(T val) {
            _a = val;
        }

        void minus_eq(T val) {
            _a -= val;
        }

        void plus_eq(T val) {
            _a += val;
        }
};

template <typename T>
static void print_atomicwrapper_v(std::vector<atomicwrapper<T>>& v ) {

#ifndef BIND
    std::cout << "[ ";
    for (typename std::vector<atomicwrapper<T>>::iterator it =
            v.begin(); it != v.end(); ++it)  {
        std::cout << it->get() << " ";
    }
    std::cout << "]\n";
#endif
}

class bitmapvector {

    // TODO

};

} // End namespace fg
#endif
