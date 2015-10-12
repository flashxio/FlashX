#include "FG2Tpetra.h"

size_t edge_data_size = 0;

ArrayRCP<size_t> getNumEntriesPerRow(fg::vertex_index::ptr index,
		size_t first_row, size_t num_local_rows, fg::edge_type type)
{
	fg::in_mem_query_vertex_index::ptr query_index
		= fg::in_mem_query_vertex_index::create(index, false);

	ArrayRCP<size_t> numEntries(num_local_rows);
	for (size_t i = 0; i < num_local_rows; i++)
		numEntries[i] = query_index->get_num_edges(i + first_row, type);
	return numEntries;
}

RCP<crs_matrix_type> create_crs(const std::string &graph_file,
		fg::vertex_index::ptr index, fg::edge_type type,
		RCP<map_type> map, int my_rank)
{
	const size_t numMyElements = map->getNodeNumElements ();
	const global_ordinal_type gblRow0 = map->getGlobalElement(0);
	const size_t numRows = index->get_num_vertices();
	printf("first vertex in process %d: %ld\n", my_rank, gblRow0);

	ArrayRCP<const size_t> numEntriesPerRow = getNumEntriesPerRow(index,
			gblRow0, numMyElements, type);
	size_t tot_size = 0;
	for (size_t i = 0; i < numMyElements; i++)
		tot_size += fg::ext_mem_undirected_vertex::num_edges2vsize(
				numEntriesPerRow[i], edge_data_size);;

	off_t first_off;
	if (index->get_graph_header().is_directed_graph()) {
		fg::in_mem_cdirected_vertex_index::ptr cu_vindex
			= fg::in_mem_cdirected_vertex_index::create(*index);
		if (type == fg::edge_type::IN_EDGE)
			first_off = cu_vindex->get_vertex(gblRow0).get_in_off();
		else
			first_off = cu_vindex->get_vertex(gblRow0).get_out_off();
	}
	else {
		fg::in_mem_cundirected_vertex_index::ptr cu_vindex
			= fg::in_mem_cundirected_vertex_index::create(*index);
		first_off = cu_vindex->get_vertex(gblRow0).get_off();
	}

	// Read the portion of the graph image that belongs to the current process.
	printf("read %ld bytes of the graph\n", tot_size);
	std::unique_ptr<char[]> graph_data(new char[tot_size]);
	FILE *f = fopen(graph_file.c_str(), "r");
	assert(f);
	int seek_ret = fseek(f, first_off, SEEK_SET);
	assert(seek_ret == 0);
	size_t read_ret = fread(graph_data.get(), tot_size, 1, f);
	assert(read_ret == 1);
	fclose(f);

	printf("fill the CRS matrix\n");
	// Create a Tpetra sparse matrix whose rows have distribution given by the Map.
	RCP<crs_matrix_type> A (new crs_matrix_type (map, numEntriesPerRow,
				Tpetra::ProfileType::StaticProfile));
	std::vector<global_ordinal_type> cols;
	std::vector<double> vals;
	off_t local_off = 0;
	// Fill the sparse matrix, one row at a time.
	for (local_ordinal_type lclRow = 0;
			lclRow < static_cast<local_ordinal_type> (numMyElements); ++lclRow) {
		const global_ordinal_type gblRow = map->getGlobalElement (lclRow);
		const fg::ext_mem_undirected_vertex *v
			= (const fg::ext_mem_undirected_vertex *) (graph_data.get() + local_off);
		assert(v->get_id() == gblRow);
		fg::vsize_t num_edges = v->get_num_edges();
		cols.resize(num_edges);
		vals.resize(num_edges);
		for (fg::vsize_t i = 0; i < num_edges; i++) {
			cols[i] = v->get_neighbor(i);
			assert(cols[i] < numRows);
			vals[i] = 1;
		}
		assert(num_edges == numEntriesPerRow[lclRow]);
		assert(gblRow < numRows);
		A->insertGlobalValues (gblRow, cols, vals);
		local_off += fg::ext_mem_undirected_vertex::num_edges2vsize(num_edges,
				edge_data_size);
	}
	assert((size_t) local_off == tot_size);

	// Tell the sparse matrix that we are done adding entries to it.
	A->fillComplete ();
	return A;
}
