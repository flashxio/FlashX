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

typedef size_t local_ordinal_type;
typedef size_t global_ordinal_type;

typedef Tpetra::DefaultPlatform::DefaultPlatformType Platform;
typedef Tpetra::DefaultPlatform::DefaultPlatformType::NodeType  Node;
typedef Tpetra::Map<local_ordinal_type, global_ordinal_type, Node> map_type;
typedef Tpetra::CrsMatrix<double, local_ordinal_type, global_ordinal_type> crs_matrix_type;

size_t edge_data_size = 0;

ArrayRCP<size_t> getNumEntriesPerRow(fg::vertex_index::ptr index,
		size_t first_row, size_t num_local_rows)
{
	fg::in_mem_query_vertex_index::ptr query_index
		= fg::in_mem_query_vertex_index::create(index, false);

	ArrayRCP<size_t> numEntries(num_local_rows);
	for (size_t i = 0; i < num_local_rows; i++)
		numEntries[i] = query_index->get_num_edges(i + first_row,
				fg::edge_type::IN_EDGE);
	return numEntries;
}

RCP<crs_matrix_type> create_crs(fg::in_mem_graph::ptr g,
		fg::vertex_index::ptr index, RCP<map_type> map, int my_rank)
{
	std::shared_ptr<const char> graph_data = g->get_raw_data();
	const size_t numMyElements = map->getNodeNumElements ();
	const global_ordinal_type gblRow0 = map->getGlobalElement(0);
	const size_t numRows = index->get_num_vertices();

	printf("Get #entries per row\n");
	ArrayRCP<const size_t> numEntriesPerRow = getNumEntriesPerRow(index,
			gblRow0, numMyElements);
	printf("allocate CRS matrix\n");
	// Create a Tpetra sparse matrix whose rows have distribution given by the Map.
	RCP<crs_matrix_type> A (new crs_matrix_type (map, numEntriesPerRow,
				Tpetra::ProfileType::StaticProfile));

	printf("fill the CRS matrix\n");
	fg::in_mem_cundirected_vertex_index::ptr cu_vindex
		= fg::in_mem_cundirected_vertex_index::create(*index);
	printf("first vertex in process %d: %ld\n", my_rank, gblRow0);
	off_t off = cu_vindex->get_vertex(gblRow0).get_off();
	printf("The offset of the first vertex in the image is %ld\n", off);
	std::vector<global_ordinal_type> cols;
	std::vector<double> vals;
	// Fill the sparse matrix, one row at a time.
	for (local_ordinal_type lclRow = 0;
			lclRow < static_cast<local_ordinal_type> (numMyElements); ++lclRow) {
		const global_ordinal_type gblRow = map->getGlobalElement (lclRow);
		const fg::ext_mem_undirected_vertex *v
			= (const fg::ext_mem_undirected_vertex *) (graph_data.get() + off);
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
		off += fg::ext_mem_undirected_vertex::num_edges2vsize(num_edges,
				edge_data_size);
	}

	// Tell the sparse matrix that we are done adding entries to it.
	A->fillComplete ();
	return A;
}

void test_tpetra(fg::in_mem_graph::ptr g, fg::vertex_index::ptr index)
{
	typedef Tpetra::MultiVector<double, global_ordinal_type, global_ordinal_type, Node> MV;

	Platform &platform = Tpetra::DefaultPlatform::getDefaultPlatform();
	RCP<const Teuchos::Comm<int> > comm = platform.getComm();
	const int myRank = comm->getRank ();
	RCP<map_type> Map = rcp (new map_type(
				index->get_graph_header().get_num_vertices(), 0, comm));

	printf("start to convert FG format to CSR format\n");
	struct timeval start, end;
	gettimeofday(&start, NULL);
	RCP<crs_matrix_type> A = create_crs(g, index, Map, myRank);
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
	if (argc < 3) {
		fprintf(stderr, "multiply graph_file index_file\n");
		exit(1);
	}

	std::string graph_file = argv[1];
	std::string index_file = argv[2];
	fg::in_mem_graph::ptr g = fg::in_mem_graph::load_graph(graph_file);
	fg::vertex_index::ptr index = fg::vertex_index::load(index_file);

	Teuchos::GlobalMPISession mpiSession (&argc, &argv);
	printf("MPI initialized: %d, #procs: %d, rank: %d\n",
			Teuchos::GlobalMPISession::mpiIsInitialized(),
			Teuchos::GlobalMPISession::getNProc(),
			Teuchos::GlobalMPISession::getRank());

	test_tpetra(g, index);
}
