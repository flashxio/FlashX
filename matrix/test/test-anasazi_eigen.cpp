// This example computes the eigenvalues of largest magnitude of an
// eigenvalue problem $A x = \lambda x$, using Anasazi's
// implementation of the Block Davidson method.

#include "FG2Tpetra.h"

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

// ****************************************************************************
// BEGIN MAIN ROUTINE
// ****************************************************************************

using namespace fm;

std::atomic<size_t> iter_no;

RCP<map_type> Map;

class spm_function
{
public:
	typedef std::shared_ptr<const spm_function> const_ptr;

	virtual void run(const MV& x, MV& y) const = 0;
};

class eigen_operator: public Tpetra::Operator<double, local_ordinal_type, global_ordinal_type, Node>
{
	RCP<crs_matrix_type> A;
public:
	eigen_operator(RCP<crs_matrix_type> A) {
		this->A = A;
	}

	Teuchos::RCP<const map_type> getDomainMap() const {
		return A->getDomainMap();
	}

	Teuchos::RCP<const map_type> getRangeMap() const {
		return A->getRangeMap();
	}

	void apply (const MV &X, MV &Y, Teuchos::ETransp mode, double alpha,
			double beta) const {
		struct timeval start, end;
		gettimeofday(&start, NULL);
		A->apply(X, Y);
		gettimeofday(&end, NULL);
		printf("SpMM takes %.3f seconds\n", time_diff(start, end));
	}
};

class SVD_operator: public Tpetra::Operator<double, local_ordinal_type, global_ordinal_type, Node>
{
	RCP<crs_matrix_type> A;
	RCP<crs_matrix_type> tA;
public:
	SVD_operator(RCP<crs_matrix_type> A, RCP<crs_matrix_type> tA) {
		this->A = A;
		this->tA = tA;
	}

	Teuchos::RCP<const map_type> getDomainMap() const {
		return A->getDomainMap();
	}

	Teuchos::RCP<const map_type> getRangeMap() const {
		return A->getRangeMap();
	}

	void apply (const MV &X, MV &Y, Teuchos::ETransp mode, double alpha,
			double beta) const {
		struct timeval start, end;
		gettimeofday(&start, NULL);
		RCP<MV> tmp = rcp (new MV (Map, X.getNumVectors()));
		A->apply(X, *tmp);
		tA->apply(*tmp, Y);
		gettimeofday(&end, NULL);
		printf("SVD SpMM takes %.3f seconds\n", time_diff(start, end));
	}
};

typedef Tpetra::Operator<double, global_ordinal_type, global_ordinal_type, Node> OP;
typedef Anasazi::OperatorTraits<double, MV, OP> OPT;

void compute_eigen(RCP<OP> A, int nev, const std::string &solver,
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
#if 0
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
#endif
}

void print_usage()
{
	fprintf(stderr, "eigensolver matrix_file index_file [options]\n");
	fprintf(stderr, "-b block_size_start\n");
	fprintf(stderr, "-B block_size_end\n");
	fprintf(stderr, "-n num_blocks_start\n");
	fprintf(stderr, "-N num_blocks_end\n");
	fprintf(stderr, "-e nev_start\n");
	fprintf(stderr, "-E nev_end\n");
	fprintf(stderr, "-s solver: Davidson, KrylovSchur, LOBPCG\n");
	fprintf(stderr, "-t tolerance\n");
}

int main (int argc, char *argv[])
{
	int opt;
	int num_opts = 0;
	size_t blockSizeStart = 4;
	size_t blockSizeEnd = 0;
	size_t numBlocksStart = 8;
	size_t numBlocksEnd = 0;
	size_t nevStart = 8;
	size_t nevEnd = 0;
	std::string solver = "LOBPCG";
	double tol = 1e-8;
	while ((opt = getopt(argc, argv, "b:B:n:N:s:t:e:E:")) != -1) {
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
			case 'e':
				nevStart = atoi(optarg);
				num_opts++;
				break;
			case 'E':
				nevEnd = atoi(optarg);
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
	if (nevEnd == 0)
		nevEnd = nevStart;

	argv += 1 + num_opts;
	argc -= 1 + num_opts;
	if (argc < 2) {
		print_usage();
		exit(1);
	}

	std::string graph_file = argv[0];
	std::string index_file = argv[1];

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

	RCP<OP> A;
	if (index->get_graph_header().is_directed_graph()) {
		RCP<crs_matrix_type> mat = create_crs(graph_file, index,
				fg::edge_type::OUT_EDGE, Map, my_rank);
		RCP<crs_matrix_type> t_mat = create_crs(graph_file, index,
				fg::edge_type::IN_EDGE, Map, my_rank);
		A = rcp(new SVD_operator(mat, t_mat));
	}
	else {
		RCP<crs_matrix_type> mat = create_crs(graph_file, index,
				fg::edge_type::OUT_EDGE, Map, my_rank);
		A = mat;
//		A = rcp(new eigen_operator(mat));
	}

	for (size_t nev = nevStart; nev <= nevEnd; nev *= 2) {
		for (size_t blockSize = blockSizeStart; blockSize <= blockSizeEnd;
				blockSize *= 2) {
			for (size_t numBlocks = numBlocksStart; numBlocks <= numBlocksEnd;
					numBlocks *= 2)
				compute_eigen(A, nev, solver, blockSize, numBlocks, tol, my_rank);
		}
	}

	return 0;
}
