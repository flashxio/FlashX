/**
 * Copyright 2014 Open Connectome Project (http://openconnecto.me)
 * Written by Disa Mhembere (disa)
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
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#include <fstream>
#include "worker_thread.h"
#include "graph_engine.h"

using namespace fg;

namespace {
    void make_index(std::string outfn, size_t num_vert, size_t num_cols) {
        std::cout << "Creating header " << std::endl;
        graph_header header(UNDIRECTED, num_vert, 
                (num_vert*num_cols), 8);

        std::cout << "Creating vertices vector" << std::endl;
        // TODO: vertices verify content and vertex_entry_type (currently double)
        std::vector<double> vertices;
        vertices.resize(num_vert);
        vertex_index::ptr index = vertex_index_temp<double>::create(
                header, vertices);

        std::cout << "Opening file" << std::endl;
        std::ofstream outfile;
        outfile.open(outfn);
        outfile.open(outfn, std::ios::binary | std::ios::trunc | std::ios::out);
        std::cout << "Writing " << sizeof(*index) << "bytes .." << std::endl;
        outfile.write((char*)index.get(), sizeof(*index));
        outfile.close();
        std::cout << "Closed file ..." << std::endl;
    }

}

int main(int argc, char* argv[]) {

    if (argc < 4) {
        fprintf(stderr, "usage: ./mkindex outfilename num_vertices num_cols\n");
        exit(EXIT_FAILURE);
    }

    make_index(argv[1], atol(argv[2]), atol(argv[3]));

    return EXIT_SUCCESS;
}
