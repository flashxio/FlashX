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

#include "fake_index.h"

using namespace fg;

int main(int argc, char* argv[]) {

    if (argc < 4) {
        fprintf(stderr, "usage: ./make_page_row_header filename num_rows num_cols\n");
        return EXIT_FAILURE;
    }

    size_t nrow = atol(argv[2]);
    size_t ncol = atol(argv[3]);
    printf("Creating header for %lu X %lu..\n", nrow, ncol);
    graph_header header = make_graph_header(nrow, ncol, argv[1]);
    FILE *f = fopen(argv[1], "wb");
    fwrite(&header, sizeof(header), 1, f);
    fclose(f);
    printf("Writing header!\n");

    return EXIT_SUCCESS;
}
