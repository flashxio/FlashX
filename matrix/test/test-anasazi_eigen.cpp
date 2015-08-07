// This example computes the eigenvalues of largest magnitude of an
// eigenvalue problem $A x = \lambda x$, using Anasazi's
// implementation of the Block Davidson method.

#include "sparse_matrix.h"
#include "vertex.h"
#include "vertex_index.h"
#include "in_mem_storage.h"
#include "io_interface.h"

// Include header for block Davidson eigensolver
#include "AnasaziBlockDavidsonSolMgr.hpp"
// Include header for LOBPCG eigensolver
#include "AnasaziLOBPCGSolMgr.hpp"
// Include header for block Davidson eigensolver
#include "AnasaziBlockKrylovSchurSolMgr.hpp"
// Include header to define eigenproblem Ax = \lambda*x
#include "AnasaziBasicEigenproblem.hpp"
#include "AnasaziOperator.hpp"
#include "AnasaziTpetraAdapter.hpp"
#include "Tpetra_Map.hpp"
#include "Tpetra_MultiVector.hpp"
#include "Tpetra_DefaultPlatform.hpp"
#include "Tpetra_CrsMatrix_decl.hpp"
#include "Tpetra_CrsMatrix_def.hpp"

// ****************************************************************************
// BEGIN MAIN ROUTINE
// ****************************************************************************

using namespace fm;

using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::ArrayRCP;
using std::cerr;
using std::cout;
using std::endl;

std::atomic<size_t> iter_no;

typedef size_t local_ordinal_type;
typedef size_t global_ordinal_type;

typedef Tpetra::DefaultPlatform::DefaultPlatformType::NodeType  Node;
typedef Tpetra::MultiVector<double, global_ordinal_type, global_ordinal_type, Node> MV;
typedef Tpetra::Operator<double, global_ordinal_type, global_ordinal_type, Node> OP;
typedef Anasazi::OperatorTraits<double, MV, OP> OPT;

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

RCP<crs_matrix_type> create_crs(const std::string &graph_file,
		fg::vertex_index::ptr index, RCP<map_type> map, int my_rank)
{
	fg::in_mem_cundirected_vertex_index::ptr cu_vindex
		= fg::in_mem_cundirected_vertex_index::create(*index);

	const size_t numMyElements = map->getNodeNumElements ();
	const global_ordinal_type gblRow0 = map->getGlobalElement(0);
	const size_t numRows = index->get_num_vertices();
	printf("first vertex in process %d: %ld\n", my_rank, gblRow0);

	ArrayRCP<const size_t> numEntriesPerRow = getNumEntriesPerRow(index,
			gblRow0, numMyElements);
	size_t tot_size = 0;
	for (size_t i = 0; i < numMyElements; i++)
		tot_size += fg::ext_mem_undirected_vertex::num_edges2vsize(
				numEntriesPerRow[i], edge_data_size);;
	off_t off = cu_vindex->get_vertex(gblRow0).get_off();

	// Read the portion of the graph image that belongs to the current process.
	std::unique_ptr<char[]> graph_data(new char[tot_size]);
	FILE *f = fopen(graph_file.c_str(), "r");
	assert(f);
	int seek_ret = fseek(f, off, SEEK_SET);
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

RCP<map_type> Map;

void compute_eigen(RCP<crs_matrix_type> A, int nev, const std::string &solver,
		int blockSize, int numBlocks, double tol, int my_rank)
{
	// Here, Scalar is double, MV is Tpetra::MultiVector, and OP is FMTp_Operator.
	typedef Anasazi::MultiVecTraits<double, MV> MVT;

	// Set eigensolver parameters.
	const int maxRestarts = 100; // maximum number of restart cycles
	const int maxIters = 500;

	if (my_rank == 0) {
		printf("solver: %s\n", solver.c_str());
		printf("block size: %d\n", blockSize);
		printf("#blocks: %d\n", numBlocks);
		printf("tol: %g\n", tol);
	}

	// Create a set of initial vectors to start the eigensolver.
	// This needs to have the same number of columns as the block size.
	RCP<MV> ivec = rcp (new MV (Map, (size_t) blockSize));
	ivec->randomize ();

	// Create the eigenproblem.  This object holds all the stuff about
	// your problem that Anasazi will see.  In this case, it knows about
	// the matrix A and the inital vectors.
	RCP<Anasazi::BasicEigenproblem<double, MV, OP> > problem =
		rcp (new Anasazi::BasicEigenproblem<double, MV, OP> (A, ivec));

	// Tell the eigenproblem that the operator A is symmetric.
	problem->setHermitian (true);

	// Set the number of eigenvalues requested
	problem->setNEV (nev);

	// Tell the eigenproblem that you are finishing passing it information.
	const bool boolret = problem->setProblem();
	if (boolret != true) {
		cerr << "Anasazi::BasicEigenproblem::setProblem() returned an error." << endl;
		return;
	}

	// Create a ParameterList, to pass parameters into the Block
	// Davidson eigensolver.
	Teuchos::ParameterList anasaziPL;
	anasaziPL.set ("Which", "LM");
	anasaziPL.set ("Block Size", blockSize);
	anasaziPL.set ("Convergence Tolerance", tol);
	anasaziPL.set ("Verbosity", Anasazi::Errors + Anasazi::Warnings +
			Anasazi::TimingDetails + Anasazi::FinalSummary);

	Anasazi::ReturnType returnCode;
	// Create the Block Davidson eigensolver.
	if (solver == "Davidson") {
		anasaziPL.set ("Num Blocks", numBlocks);
		anasaziPL.set ("Maximum Restarts", maxRestarts);

		Anasazi::BlockDavidsonSolMgr<double, MV, OP> anasaziSolver (problem,
				anasaziPL);
		// Solve the eigenvalue problem.
		returnCode = anasaziSolver.solve ();
	}
	else if (solver == "KrylovSchur") {
		anasaziPL.set ("Num Blocks", numBlocks);
		anasaziPL.set ("Maximum Restarts", maxRestarts);

		Anasazi::BlockKrylovSchurSolMgr<double, MV, OP> anasaziSolver (problem,
				anasaziPL);
		// Solve the eigenvalue problem.
		returnCode = anasaziSolver.solve ();
	}
	else if (solver == "LOBPCG") {
		anasaziPL.set ("Maximum Iterations", maxIters);
		anasaziPL.set ("Full Ortho", true);
		anasaziPL.set ("Use Locking", true);

		// Create the LOBPCG eigensolver.
		Anasazi::LOBPCGSolMgr<double, MV, OP> anasaziSolver (problem, anasaziPL);
		// Solve the eigenvalue problem.
		returnCode = anasaziSolver.solve ();
	}
	else
		assert(0);

	if (returnCode != Anasazi::Converged && my_rank == 0) {
		cout << "Anasazi eigensolver did not converge." << endl;
	}
	// Get the eigenvalues and eigenvectors from the eigenproblem.
	Anasazi::Eigensolution<double,MV> sol = problem->getSolution ();

	// Anasazi returns eigenvalues as Anasazi::Value, so that if
	// Anasazi's Scalar type is real-valued (as it is in this case), but
	// some eigenvalues are complex, you can still access the
	// eigenvalues correctly.  In this case, there are no complex
	// eigenvalues, since the matrix pencil is symmetric.
	std::vector<Anasazi::Value<double> > evals = sol.Evals;
	RCP<MV> evecs = sol.Evecs;

	// Compute residuals.
	std::vector<double> normR (sol.numVecs);
	if (sol.numVecs > 0) {
		Teuchos::SerialDenseMatrix<int,double> T (sol.numVecs, sol.numVecs);
		MV tempAevec (Map, sol.numVecs);
		T.putScalar (0.0);
		for (int i=0; i<sol.numVecs; ++i) {
			T(i,i) = evals[i].realpart;
		}
		OPT::Apply( *A, *evecs, tempAevec ); 
		MVT::MvTimesMatAddMv (-1.0, *evecs, T, 1.0, tempAevec);
		MVT::MvNorm (tempAevec, normR);
	}

	if (my_rank == 0) {
		// Print the results on MPI process 0.
		cout << "Solver manager returned "
			<< (returnCode == Anasazi::Converged ? "converged." : "unconverged.")
			<< endl << endl
			<< "------------------------------------------------------" << endl
			<< std::setw(16) << "Eigenvalue"
			<< std::setw(18) << "Direct Residual"
			<< endl
			<< "------------------------------------------------------" << endl;
		for (int i=0; i<sol.numVecs; ++i) {
			cout << std::setw(16) << evals[i].realpart
				<< std::setw(18) << normR[i] / evals[i].realpart
				<< endl;
		}
		cout << "------------------------------------------------------" << endl;
		cout << iter_no << endl;
	}
}

void print_usage()
{
	fprintf(stderr, "eigensolver matrix_file index_file nev [options]\n");
	fprintf(stderr, "-b block_size_start\n");
	fprintf(stderr, "-B block_size_end\n");
	fprintf(stderr, "-n num_blocks_start\n");
	fprintf(stderr, "-N num_blocks_end\n");
	fprintf(stderr, "-s solver: Davidson, KrylovSchur, LOBPCG\n");
	fprintf(stderr, "-t tolerance\n");
	fprintf(stderr, "-S: run SVD\n");
}

int main (int argc, char *argv[])
{
	int opt;
	int num_opts = 0;
	size_t blockSizeStart = 4;
	size_t blockSizeEnd = 0;
	size_t numBlocksStart = 8;
	size_t numBlocksEnd = 0;
	std::string solver = "LOBPCG";
	double tol = 1e-8;
	while ((opt = getopt(argc, argv, "b:B:n:N:s:t:")) != -1) {
		num_opts++;
		switch (opt) {
			case 'b':
				blockSizeStart = atoi(optarg);
				num_opts++;
				break;
			case 'B':
				blockSizeEnd = atoi(optarg);
				num_opts++;
				break;
			case 'n':
				numBlocksStart = atoi(optarg);
				num_opts++;
				break;
			case 'N':
				numBlocksEnd = atoi(optarg);
				num_opts++;
				break;
			case 's':
				solver = optarg;
				num_opts++;
				break;
			case 't':
				tol = atof(optarg);
				num_opts++;
				break;
			default:
				print_usage();
				abort();
		}
	}
	if (blockSizeEnd == 0)
		blockSizeEnd = blockSizeStart;
	if (numBlocksEnd == 0)
		numBlocksEnd = numBlocksStart;

	argv += 1 + num_opts;
	argc -= 1 + num_opts;
	if (argc < 3) {
		print_usage();
		exit(1);
	}

	std::string graph_file = argv[0];
	std::string index_file = argv[1];
	int nev = atoi(argv[3]); // number of eigenvalues for which to solve;

	fg::vertex_index::ptr index = fg::vertex_index::load(index_file);
	Teuchos::GlobalMPISession mpiSession (&argc, &argv);

	//
	// Set up the test problem.
	//
	typedef Tpetra::DefaultPlatform::DefaultPlatformType Platform;
	Platform &platform = Tpetra::DefaultPlatform::getDefaultPlatform();
	RCP<const Teuchos::Comm<int> > comm = platform.getComm();
	Map = rcp (new map_type(index->get_graph_header().get_num_vertices(),
				0, comm));
	int my_rank = comm->getRank();
	RCP<crs_matrix_type> A = create_crs(graph_file, index, Map, my_rank);
	for (size_t blockSize = blockSizeStart; blockSize <= blockSizeEnd;
			blockSize *= 2) {
		for (size_t numBlocks = numBlocksStart; numBlocks <= numBlocksEnd;
				numBlocks *= 2)
			compute_eigen(A, nev, solver, blockSize, numBlocks, tol, my_rank);
	}

	return 0;
}
