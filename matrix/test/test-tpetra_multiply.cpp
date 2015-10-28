#include "FG2Tpetra.h"

using namespace fm;

void test_tpetra(const std::string graph_file, fg::vertex_index::ptr index)
{
	typedef Tpetra::DefaultPlatform::DefaultPlatformType Platform;

	Platform &platform = Tpetra::DefaultPlatform::getDefaultPlatform();
	RCP<const Teuchos::Comm<int> > comm = platform.getComm();
	const int myRank = comm->getRank ();
	RCP<map_type> Map = rcp (new map_type(
				index->get_graph_header().get_num_vertices(), 0, comm));

	struct timeval start, end;
	gettimeofday(&start, NULL);
	RCP<crs_matrix_type> A = create_crs(graph_file, index,
			fg::edge_type::OUT_EDGE, Map, myRank,
			index->get_graph_header().get_edge_data_size());
	gettimeofday(&end, NULL);
	printf("conversion takes %.3f seconds\n", time_diff(start, end));

	for (size_t num_vecs = 1; num_vecs <= 16; num_vecs *= 2) {
		if (myRank == 0)
			printf("SpMM: #vecs: %ld\n", num_vecs);
		RCP<MV> evecs = rcp(new MV(Map, num_vecs));
		RCP<MV> res = rcp(new MV(Map, num_vecs));
		evecs->randomize ();

		printf("start SpMM, A on proc %d has %ld rows and %ld cols\n",
				myRank, A->getNodeNumRows(), A->getNodeNumCols());
		struct timeval start, end;
		gettimeofday(&start, NULL);
		A->apply (*evecs, *res, Teuchos::NO_TRANS);
		gettimeofday(&end, NULL);
		printf("SpMM(%ld) takes %.3fs\n", num_vecs, time_diff(start, end));
		comm->barrier();
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
	fg::vertex_index::ptr index = fg::vertex_index::load(index_file);

	Teuchos::GlobalMPISession mpiSession (&argc, &argv);
	printf("MPI initialized: %d, #procs: %d, rank: %d\n",
			Teuchos::GlobalMPISession::mpiIsInitialized(),
			Teuchos::GlobalMPISession::getNProc(),
			Teuchos::GlobalMPISession::getRank());

	test_tpetra(graph_file, index);
}
