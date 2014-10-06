/**
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
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#include "FGlib.h"
#include "graph.h"

FG_vector<float>::ptr compute_transitivity(FG_graph::ptr fg) {
	printf("Transitivity starts ...\n");
	struct timeval start, end;
	gettimeofday(&start, NULL);

	FG_vector<vsize_t>::ptr deg_v = get_degree(fg, BOTH_EDGES);
	FG_vector<size_t>::ptr ss_v = compute_local_scan(fg);
	FG_vector<float>::ptr trans_v = FG_vector<float>::create(deg_v->get_size());

#pragma omp parallel for
	for (vsize_t i = 0; i < deg_v->get_size(); i++) {
		size_t curr_deg = deg_v->get(i);
		curr_deg > 0 ? trans_v->set(i, (ss_v->get(i)/(float)((curr_deg+1) * (curr_deg)))):
					   trans_v->set(i, 0);
	}
	
	gettimeofday(&end, NULL);
	printf("Transitivity takes %f seconds\n", time_diff(start, end));
	return trans_v; 
}
