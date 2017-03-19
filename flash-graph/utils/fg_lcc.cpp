#include "FGlib.h"
#include "in_mem_storage.h"

#include "sparse_matrix.h"
#include "fg_utils.h"

#include "col_vec.h"

using namespace fm;

std::pair<fg::vertex_id_t, size_t> get_max_cid(fm::vector::ptr cc_ids)
{
	fm::col_vec::ptr vec = fm::col_vec::create(cc_ids);
	fm::bulk_operate::const_ptr add = bulk_operate::conv2ptr(
			get_scalar_type<size_t>().get_basic_ops().get_add());
	fm::agg_operate::const_ptr sum = fm::agg_operate::create(add);;

	fm::data_frame::ptr counts = vec->groupby(sum, true);
	fm::detail::mem_vec_store::const_ptr ids
		= std::dynamic_pointer_cast<const fm::detail::mem_vec_store>(
				counts->get_vec(0));
	fm::detail::mem_vec_store::const_ptr cnts
		= std::dynamic_pointer_cast<const fm::detail::mem_vec_store>(
				counts->get_vec(1));
	assert(cnts->get_type() == get_scalar_type<size_t>());
	size_t idx = 0;
	size_t max_size = cnts->get<size_t>(0);
	for (size_t i = 1; i < cnts->get_length(); i++)
		if (cnts->get<size_t>(i) > max_size) {
			max_size = cnts->get<size_t>(i);
			idx = i;
		}
	fg::vertex_id_t id = ids->get<fg::vertex_id_t>(idx);
	return std::pair<fg::vertex_id_t, size_t>(id, max_size);
}

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
	fm::vector::ptr cc_ids = is_directed ? fg::compute_wcc(g) : fg::compute_cc(g);

	std::pair<fg::vertex_id_t, size_t> max_cid = get_max_cid(cc_ids);
	fm::detail::mem_vec_store::const_ptr id_store
		= std::dynamic_pointer_cast<const fm::detail::mem_vec_store>(
				cc_ids->get_raw_store());
	// we can use which().
	std::vector<fg::vertex_id_t> lcc_vids;
	for (size_t i = 0; i < cc_ids->get_length(); i++)
		if (id_store->get<fg::vertex_id_t>(i) == max_cid.first)
			lcc_vids.push_back(i);

	fg::FG_graph::ptr lcc = fg::fetch_subgraph(g, lcc_vids, new_graph_name, true);
	if (lcc) {
		const fg::graph_header &header = lcc->get_graph_header();
		printf("lcc has %ld vertices and %ld edges\n", header.get_num_vertices(),
				header.get_num_edges());
		auto graph_data = lcc->get_graph_data();
		graph_data->dump(new_graph_name + ".adj");
		lcc->get_index_data()->dump(new_graph_name + ".index");
	}

	destroy_flash_matrix();
	fg::graph_engine::destroy_flash_graph();
}
