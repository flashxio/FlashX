#include "sparse_matrix.h"
#include "vertex.h"
#include "in_mem_storage.h"
#include "io_interface.h"

#include "Tpetra_Map.hpp"
#include "Tpetra_MultiVector.hpp"
#include "Tpetra_DefaultPlatform.hpp"
#include "Tpetra_CrsMatrix_decl.hpp"
#include "Tpetra_CrsMatrix_def.hpp"

using namespace fm;

using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::ArrayRCP;
using std::cerr;
using std::cout;
using std::endl;

typedef fg::vsize_t local_ordinal_type;
typedef fg::vsize_t global_ordinal_type;

typedef Tpetra::DefaultPlatform::DefaultPlatformType Platform;
typedef Tpetra::DefaultPlatform::DefaultPlatformType::NodeType  Node;
typedef Tpetra::Map<local_ordinal_type, global_ordinal_type, Node> map_type;
typedef Tpetra::CrsMatrix<double, local_ordinal_type, global_ordinal_type> crs_matrix_type;

size_t edge_data_size = 0;

ArrayRCP<size_t> getNumEntriesPerRow(fg::vertex_index::ptr index,
		size_t &num_rows)
{
	fg::in_mem_query_vertex_index::ptr query_index
		= fg::in_mem_query_vertex_index::create(index, false);

	ArrayRCP<size_t> numEntries(index->get_num_vertices());
	for (size_t i = 0; i < index->get_num_vertices(); i++)
		numEntries[i] = query_index->get_num_edges(i, fg::edge_type::IN_EDGE);

	num_rows = index->get_num_vertices();
	return numEntries;
}

RCP<crs_matrix_type> create_crs(fg::in_mem_graph::ptr g,
		fg::vertex_index::ptr index, RCP<map_type> map)
{
	safs::io_interface::ptr io = create_io(g->create_io_factory(),
			thread::get_curr_thread());
	const size_t numMyElements = map->getNodeNumElements ();

	printf("Get #entries per row\n");
	size_t numRows;
	ArrayRCP<const size_t> numEntriesPerRow = getNumEntriesPerRow(index, numRows);
	assert(numRows == numMyElements);
	printf("allocate CRS matrix\n");
	// Create a Tpetra sparse matrix whose rows have distribution given by the Map.
	RCP<crs_matrix_type> A (new crs_matrix_type (map, numEntriesPerRow,
				Tpetra::ProfileType::StaticProfile));

	printf("fill the CRS matrix\n");
	const size_t vheader_size = fg::ext_mem_undirected_vertex::get_header_size();
	char vheader_buf[vheader_size];
	const fg::ext_mem_undirected_vertex *vheader
		= (const fg::ext_mem_undirected_vertex *) vheader_buf;
	off_t off = fg::graph_header::get_header_size();
	std::vector<global_ordinal_type> cols;
	std::vector<double> vals;
	// Fill the sparse matrix, one row at a time.
	for (local_ordinal_type lclRow = 0;
			lclRow < static_cast<local_ordinal_type> (numMyElements); ++lclRow) {
		io->access(vheader_buf, off, vheader_size, READ);
		fg::vsize_t num_edges = vheader->get_num_edges();
		assert(num_edges > 0);
		size_t size = fg::ext_mem_undirected_vertex::num_edges2vsize(num_edges,
				edge_data_size);
		std::unique_ptr<char[]> buf(new char[size]);
		io->access(buf.get(), off, size, READ);

		const fg::ext_mem_undirected_vertex *v
			= (const fg::ext_mem_undirected_vertex *) buf.get();
		cols.resize(num_edges);
		vals.resize(num_edges);
		for (fg::vsize_t i = 0; i < num_edges; i++) {
			cols[i] = v->get_neighbor(i);
			assert(cols[i] < numRows);
			vals[i] = 1;
		}
		assert(num_edges == numEntriesPerRow[lclRow]);
		const global_ordinal_type gblRow = map->getGlobalElement (lclRow);
		assert(gblRow < numRows);
		A->insertGlobalValues (gblRow, cols, vals);
		off += size;
	}

	// Tell the sparse matrix that we are done adding entries to it.
	A->fillComplete ();
	return A;
}

void test_tpetra(fg::FG_graph::ptr fg)
{
	typedef Tpetra::MultiVector<double, global_ordinal_type, global_ordinal_type, Node> MV;

	Platform &platform = Tpetra::DefaultPlatform::getDefaultPlatform();
	RCP<const Teuchos::Comm<int> > comm = platform.getComm();
	RCP<map_type> Map = rcp (new map_type(
				fg->get_graph_header().get_num_vertices(), 0, comm));

	printf("start to convert FG format to CSR format\n");
	struct timeval start, end;
	gettimeofday(&start, NULL);
	RCP<crs_matrix_type> A = create_crs(fg->get_graph_data(), fg->get_index_data(), Map);
	gettimeofday(&end, NULL);
	printf("conversion takes %.3f seconds\n", time_diff(start, end));

	for (size_t num_vecs = 1; num_vecs <= 16; num_vecs *= 2) {
		printf("SpMM: #vecs: %ld\n", num_vecs);
		RCP<MV> evecs = rcp(new MV(Map, num_vecs));
		RCP<MV> res = rcp(new MV(Map, num_vecs));
		evecs->randomize ();

		printf("start SpMM, A has %ld rows and %ld cols\n", A->getGlobalNumRows(), A->getGlobalNumCols());
		struct timeval start, end;
		gettimeofday(&start, NULL);
		A->apply (*evecs, *res, Teuchos::NO_TRANS);
		gettimeofday(&end, NULL);
		printf("SpMM(%ld) takes %.3fs\n", num_vecs, time_diff(start, end));
	}
}

int main(int argc, char *argv[])
{
	if (argc < 4) {
		fprintf(stderr, "eigensolver conf_file graph_file index_file\n");
		exit(1);
	}

	std::string conf_file = argv[1];
	std::string graph_file = argv[2];
	std::string index_file = argv[3];

	config_map::ptr configs = config_map::create(conf_file);
	init_flash_matrix(configs);

	fg::FG_graph::ptr fg = fg::FG_graph::create(graph_file, index_file, configs);
	test_tpetra(fg);
}
