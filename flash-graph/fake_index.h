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

#include "graph_engine.h"

using namespace fg;

namespace {
    static graph_header make_graph_header(const size_t num_vert, const size_t num_cols, 
            std::string outfn="") {
        std::cout << "Creating header " << std::endl;
        const short edge_data_size = 8;

        graph_header header(UNDIRECTED, num_vert, (num_vert*num_cols), edge_data_size);
        
        if (!outfn.empty()) {
            std::cout << "Opening file" << std::endl;
            
            FILE *f = fopen(outfn.c_str(), "w");
            assert(f);

            BOOST_VERIFY(fwrite((char*) &header, sizeof(header), 1, f));
            fclose(f);
        }
        return header;
    }

    static vertex_index::ptr make_index(const size_t num_vert, 
            const size_t num_cols, const std::string outfn="") {
        graph_header header = make_graph_header(num_vert, num_cols);

        off_t header_offset = sizeof(header);
        std::vector<vertex_offset> vertices;

        for (size_t row = 0; row < num_vert; row++) {
            vertices.push_back(header_offset + (row*(sizeof(double)*num_cols)));
        }
        vertices.push_back(header_offset + (num_vert*(sizeof(double)*num_cols)));
#if 0
        // Test print
        printf("Test print of vertices index ..\n[");
        for (std::vector<vertex_offset>::iterator it = vertices.begin();
                it != vertices.end(); ++it) {
            std::cout << " " << it->get_off();
        }
        std::cout << "]\n"; 
        // End Test print
#endif
        if (!outfn.empty()) {
            std::cout << "[WARNING]: Writing graph index to file!!\n";
            vertex_index_temp<vertex_offset>::dump(outfn, header, vertices);
        }

        return vertex_index_temp<vertex_offset>::create(header, vertices);
    }
}
