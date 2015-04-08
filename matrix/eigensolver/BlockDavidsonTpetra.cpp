// This example computes the eigenvalues of largest magnitude of an
// eigenvalue problem $A x = \lambda x$, using Anasazi's
// implementation of the Block Davidson method.

#include "sparse_matrix.h"
#include "vertex.h"
#include "in_mem_storage.h"
#include "io_interface.h"

// Include header for block Davidson eigensolver
#include "AnasaziBlockDavidsonSolMgr.hpp"
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
using std::cerr;
using std::cout;
using std::endl;

std::atomic<size_t> iter_no;

typedef fg::vsize_t local_ordinal_type;
typedef fg::vsize_t global_ordinal_type;

typedef Tpetra::DefaultPlatform::DefaultPlatformType::NodeType  Node;
typedef Tpetra::MultiVector<double, global_ordinal_type, global_ordinal_type, Node> MV;
typedef Tpetra::Operator<double, global_ordinal_type, global_ordinal_type, Node> OP;
typedef Anasazi::OperatorTraits<double, MV, OP> OPT;

typedef Tpetra::Map<local_ordinal_type, global_ordinal_type, Node> map_type;
typedef Tpetra::CrsMatrix<double, local_ordinal_type, global_ordinal_type> crs_matrix_type;

size_t edge_data_size = 0;

RCP<crs_matrix_type> create_csr(fg::in_mem_graph::ptr g, RCP<map_type> map)
{
	safs::io_interface::ptr io = create_io(g->create_io_factory(),
			thread::get_curr_thread());
	const size_t numMyElements = map->getNodeNumElements ();

	// Create a Tpetra sparse matrix whose rows have distribution given by the Map.
	RCP<crs_matrix_type> A (new crs_matrix_type (map, 0));

	const size_t vheader_size = fg::ext_mem_undirected_vertex::get_header_size();
	char vheader_buf[vheader_size];
	const fg::ext_mem_undirected_vertex *vheader
		= (const fg::ext_mem_undirected_vertex *) vheader_buf;
	off_t off = fg::graph_header::get_header_size();
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
		std::vector<global_ordinal_type> cols(num_edges);
		std::vector<double> vals(num_edges);
		for (fg::vsize_t i = 0; i < num_edges; i++) {
			cols[i] = v->get_neighbor(i);
			vals[i] = 1;
		}
		const global_ordinal_type gblRow = map->getGlobalElement (lclRow);
		A->insertGlobalValues (gblRow, cols, vals);
		off += size;
	}

	// Tell the sparse matrix that we are done adding entries to it.
	A->fillComplete ();
	return A;
}

int main (int argc, char *argv[])
{

	if (argc < 5) {
		fprintf(stderr, "eigensolver conf_file graph_file index_file nev\n");
		exit(1);
	}

	std::string conf_file = argv[1];
	std::string graph_file = argv[2];
	std::string index_file = argv[3];
	int nev = atoi(argv[4]); // number of eigenvalues for which to solve

	// Anasazi solvers have the following template parameters:
	//
	//   - Scalar: The type of dot product results.
	//   - MV: The type of (multi)vectors.
	//   - OP: The type of operators (functions from multivector to
	//     multivector).  A matrix (like Epetra_CrsMatrix) is an example
	//     of an operator; an Ifpack preconditioner is another example.
	//
	// Here, Scalar is double, MV is Tpetra::MultiVector, and OP is FMTp_Operator.
	typedef Anasazi::MultiVecTraits<double, MV> MVT;
	typedef Tpetra::DefaultPlatform::DefaultPlatformType Platform;

	//
	// Set up the test problem.
	//
	config_map::ptr configs = config_map::create(conf_file);
	init_flash_matrix(configs);

	fg::FG_graph::ptr fg = fg::FG_graph::create(graph_file, index_file, configs);

	// Set eigensolver parameters.
	const double tol = 1.0e-12; // convergence tolerance
	const int blockSize = nev + 1; // block size (number of eigenvectors processed at once)
	const int numBlocks = 3; // restart length
	const int maxRestarts = 100; // maximum number of restart cycles

	Platform &platform = Tpetra::DefaultPlatform::getDefaultPlatform();
	RCP<const Teuchos::Comm<int> > comm = platform.getComm();
	RCP<map_type> Map = rcp (new map_type(
				fg->get_graph_header().get_num_vertices(), 0, comm));
	RCP<crs_matrix_type> A = create_csr(fg->get_graph_data(), Map);
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
		return -1;
	}

	// Create a ParameterList, to pass parameters into the Block
	// Davidson eigensolver.
	Teuchos::ParameterList anasaziPL;
	anasaziPL.set ("Which", "LM");
	anasaziPL.set ("Block Size", blockSize);
	anasaziPL.set ("Num Blocks", numBlocks);
	anasaziPL.set ("Maximum Restarts", maxRestarts);
	anasaziPL.set ("Convergence Tolerance", tol);
	anasaziPL.set ("Verbosity", Anasazi::Errors + Anasazi::Warnings +
			Anasazi::TimingDetails + Anasazi::FinalSummary);

	// Create the Block Davidson eigensolver.
	Anasazi::BlockDavidsonSolMgr<double, MV, OP> anasaziSolver (problem, anasaziPL);

	// Solve the eigenvalue problem.
	//
	// Note that creating the eigensolver is separate from solving it.
	// After creating the eigensolver, you may call solve() multiple
	// times with different parameters or initial vectors.  This lets
	// you reuse intermediate state, like allocated basis vectors.
	Anasazi::ReturnType returnCode = anasaziSolver.solve ();
	if (returnCode != Anasazi::Converged) {
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

	return 0;
}
