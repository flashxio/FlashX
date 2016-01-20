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
        fprintf(stderr, "[ERROR]: usage => ./test_fake_index num_vert num_cols output_fn");
        exit(EXIT_FAILURE);
    }

    // TODO: Do something more test-y
    printf("Testing fake index ...\n");
    make_index(atol(argv[2]), atol(argv[3]), argv[1]);
    printf("Fake index tested correctly!\n");
    return EXIT_SUCCESS;
}
