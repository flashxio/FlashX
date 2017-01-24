#include "FGlib.h"
#include "in_mem_storage.h"

#include "sparse_matrix.h"
#include "fg_utils.h"

using namespace fm;

int main(int argc, char *argv[])
{
	if (argc < 5) {
		fprintf(stderr, "fg_lcc conf_file graph_file index_file new_graph\n");
		exit(1);
	}

	std::string conf_file = argv[1];
	std::string graph_file = argv[2];
	std::string index_file = argv[3];
	std::string new_graph_name = argv[4];


	config_map::ptr configs = config_map::create(conf_file);
	fg::graph_engine::init_flash_graph(configs);
	init_flash_matrix(configs);

	fg::FG_graph::ptr g = fg::FG_graph::create(graph_file, index_file, configs);
	fg::vertex_index::ptr vindex = g->get_index_data();
	printf("The graph has %ld vertices and %ld edges\n",
			vindex->get_num_vertices(), vindex->get_graph_header().get_num_edges());

	// Get the vertices in the largest connected component.
	bool is_directed = vindex->get_graph_header().is_directed_graph();
	fg::FG_vector<fg::vertex_id_t>::ptr cc_ids
		= is_directed ? fg::compute_wcc(g) : fg::compute_cc(g);
	fg::count_map<fg::vertex_id_t> counts;
	cc_ids->count_unique(counts);
	std::pair<fg::vertex_id_t, size_t> max_cid = counts.get_max_count();
	std::vector<fg::vertex_id_t> lcc_vids;
	for (size_t i = 0; i < cc_ids->get_size(); i++)
		if (cc_ids->get(i) == max_cid.first)
			lcc_vids.push_back(i);

	fg::FG_graph::ptr lcc = fg::fetch_subgraph(g, lcc_vids, new_graph_name, false);
	if (lcc) {
		const fg::graph_header &header = lcc->get_graph_header();
		printf("lcc has %ld vertices and %ld edges\n", header.get_num_vertices(),
				header.get_num_edges());
		auto graph_data = lcc->get_graph_data();
		if (graph_data) {
			graph_data->dump(new_graph_name + ".adj");
			lcc->get_index_data()->dump(new_graph_name + ".index");
		}
	}

	destroy_flash_matrix();
	fg::graph_engine::destroy_flash_graph();
}
