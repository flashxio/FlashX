// Include header for block Davidson eigensolver
#include "AnasaziBlockDavidsonSolMgr.hpp"
// Include header for LOBPCG eigensolver
#include "AnasaziLOBPCGSolMgr.hpp"
// Include header for block Davidson eigensolver
#include "AnasaziBlockKrylovSchurSolMgr.hpp"

// Include header to define eigenproblem Ax = \lambda*x
#include "AnasaziBasicEigenproblem.hpp"
#include "AnasaziOperator.hpp"
#include "FM_MultiVector.h"

#include "sparse_matrix.h"
#include "matrix_stats.h"

#include "eigensolver.h"
#include "block_dense_matrix.h"

// ****************************************************************************
// BEGIN MAIN ROUTINE
// ****************************************************************************

using namespace fm;

using std::cerr;
using std::cout;
using std::endl;

using namespace fm::eigen;

namespace fm
{
namespace eigen
{
extern size_t num_cached_mats;
}
}

namespace Anasazi
{

template<>
class OperatorTraits <double, FM_MultiVector<double>, spm_function>
{
public:
	/*! \brief This method takes the FM_MultiVector \c x and
	 * applies the spm_function \c Op to it resulting in the FM_MultiVector \c y.
	 */
	static void Apply ( const spm_function& Op,
			const FM_MultiVector<double>& x, FM_MultiVector<double>& y ) {
		BOOST_LOG_TRIVIAL(info) << get_curr_time_str() << ":" << boost::format(
				"SpMM: y(%1%) = A * x(%2%)") % y.get_name() % x.get_name();
		struct timeval start, end;
		gettimeofday(&start, NULL);
		bool out_mat_in_mem;
		if (x.get_data()->is_in_mem())
			out_mat_in_mem = true;
		else
			out_mat_in_mem = num_cached_mats > 0;
		block_multi_vector::sparse_matrix_multiply(Op, *x.get_data(), *y.get_data(),
				out_mat_in_mem);
		gettimeofday(&end, NULL);
		BOOST_LOG_TRIVIAL(info) << "SpMM takes " << time_diff(start, end)
			<< " seconds";
		y.sync_fm2ep();
	}

};

}

namespace fm
{

namespace eigen
{

eigen_options::eigen_options(int nev, std::string solver)
{
	tol = 1.0e-8;
	max_restarts = 100;
	max_iters = 500;
	this->nev = nev;
	this->solver = solver;
	which="LM";
	in_mem = true;

	if (solver == "KrylovSchur") {
		block_size = 1;
		// The KrylovSchur solver wants the number of blocks to be at least 3.
		num_blocks = std::max(nev * 2, 3);
	}
	else if (solver == "Davidson") {
		block_size = nev;
		num_blocks = 4;
	}
	else if (solver == "LOBPCG") {
		block_size = 4;
		num_blocks = 10;
	}
	else {
		BOOST_LOG_TRIVIAL(error) << "Unknown solver: " << solver;
		exit(1);
	}
}

eigen_res compute_eigen(spm_function *func, bool sym,
		struct eigen_options &_opts)
{
	using Teuchos::RCP;
	using Teuchos::rcp;

	// Anasazi solvers have the following template parameters:
	//
	//   - Scalar: The type of dot product results.
	//   - MV: The type of (multi)vectors.
	//   - OP: The type of operators (functions from multivector to
	//     multivector).  A matrix (like Epetra_CrsMatrix) is an example
	//     of an operator; an Ifpack preconditioner is another example.
	//
	// Here, Scalar is double, MV is FM_MultiVector, and OP is spm_function.
	typedef FM_MultiVector<double> MV;
	typedef spm_function OP;
	typedef Anasazi::MultiVecTraits<double, MV> MVT;

	RCP<spm_function> A = rcp(func);

	struct eigen_options opts = _opts;
	if (opts.block_size == 0)
		opts.block_size = 2 * opts.nev;
	// Set eigensolver parameters.
	const double tol = opts.tol; // convergence tolerance
	const int numBlocks = opts.num_blocks; // restart length
	const int maxRestarts = opts.max_restarts; // maximum number of restart cycles
	const int maxIters = opts.max_iters; // maximum number of iterations
	int blockSize = opts.block_size; // block size (number of eigenvectors processed at once)
	std::string solver = opts.solver;
	int nev = opts.nev;
	std::string which = opts.which;

	// Create a set of initial vectors to start the eigensolver.
	// This needs to have the same number of columns as the block size.
	RCP<MV> ivec = rcp (new MV (A->get_num_cols(), blockSize, blockSize,
				numBlocks * blockSize, opts.in_mem, opts.solver));
	ivec->Random ();

	// Create the eigenproblem.  This object holds all the stuff about
	// your problem that Anasazi will see.  In this case, it knows about
	// the matrix A and the inital vectors.
	RCP<Anasazi::BasicEigenproblem<double, MV, OP> > problem =
		rcp (new Anasazi::BasicEigenproblem<double, MV, OP> (A, ivec));

	// Tell the eigenproblem that the operator A is symmetric.
	problem->setHermitian (sym);

	// Set the number of eigenvalues requested
	problem->setNEV (nev);

	// Tell the eigenproblem that you are finishing passing it information.
	const bool boolret = problem->setProblem();
	if (boolret != true) {
		BOOST_LOG_TRIVIAL(error)
			<< "Anasazi::BasicEigenproblem::setProblem() returned an error.";
		return eigen_res();
	}

	// Create a ParameterList, to pass parameters into the Block
	// Davidson eigensolver.
	Teuchos::ParameterList anasaziPL;
	anasaziPL.set ("Which", which.c_str());
	anasaziPL.set ("Block Size", blockSize);
	anasaziPL.set ("Convergence Tolerance", tol);
	anasaziPL.set ("Verbosity", Anasazi::Errors + Anasazi::Warnings +
			Anasazi::TimingDetails + Anasazi::FinalSummary);

	Anasazi::ReturnType returnCode;
	if (solver == "Davidson") {
		anasaziPL.set ("Num Blocks", numBlocks);
		anasaziPL.set ("Maximum Restarts", maxRestarts);

		// Create the Block Davidson eigensolver.
		Anasazi::BlockDavidsonSolMgr<double, MV, OP> anasaziSolver (problem, anasaziPL);

		// Solve the eigenvalue problem.
		//
		// Note that creating the eigensolver is separate from solving it.
		// After creating the eigensolver, you may call solve() multiple
		// times with different parameters or initial vectors.  This lets
		// you reuse intermediate state, like allocated basis vectors.
		returnCode = anasaziSolver.solve ();
	}
	else if (solver == "KrylovSchur") {
		anasaziPL.set ("Num Blocks", numBlocks);
		anasaziPL.set ("Maximum Restarts", maxRestarts);

		// Create the Block Davidson eigensolver.
		Anasazi::BlockKrylovSchurSolMgr<double, MV, OP> anasaziSolver (problem, anasaziPL);

		// Solve the eigenvalue problem.
		//
		// Note that creating the eigensolver is separate from solving it.
		// After creating the eigensolver, you may call solve() multiple
		// times with different parameters or initial vectors.  This lets
		// you reuse intermediate state, like allocated basis vectors.
		returnCode = anasaziSolver.solve ();
	}
	else if (solver == "LOBPCG") {
		anasaziPL.set ("Maximum Iterations", maxIters);
		anasaziPL.set ("Full Ortho", true);
		anasaziPL.set ("Use Locking", true);

		// Create the LOBPCG eigensolver.
		Anasazi::LOBPCGSolMgr<double, MV, OP> anasaziSolver (problem, anasaziPL);

		// Solve the eigenvalue problem.
		//
		// Note that creating the eigensolver is separate from solving it.
		// After creating the eigensolver, you may call solve() multiple
		// times with different parameters or initial vectors.  This lets
		// you reuse intermediate state, like allocated basis vectors.
		returnCode = anasaziSolver.solve ();
	}
	else {
		BOOST_LOG_TRIVIAL(error) << "a wrong solver: " << solver;
		return eigen_res();
	}

	if (returnCode != Anasazi::Converged) {
		BOOST_LOG_TRIVIAL(error) << "Anasazi eigensolver did not converge.";
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
		MV tempAevec (A->get_num_rows(), sol.numVecs, evecs->get_block_size(),
				numBlocks * blockSize, opts.in_mem, opts.solver);
		T.putScalar (0.0);
		for (int i=0; i<sol.numVecs; ++i) {
			T(i,i) = evals[i].realpart;
		}
		block_multi_vector::sparse_matrix_multiply(*A, *evecs->get_data(),
				*tempAevec.get_data(), evecs->get_data()->is_in_mem());
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

	struct eigen_res res;
	res.vecs = evecs->get_data()->conv2matrix();
	res.vals.resize(sol.numVecs);
	for (int i=0; i<sol.numVecs; ++i) {
		res.vals[i] = evals[i].realpart;
		cout << std::setw(16) << evals[i].realpart
			<< std::setw(18) << normR[i] / evals[i].realpart
			<< endl;
	}
	cout << "------------------------------------------------------" << endl;
	cout << "#col writes: " << num_col_writes << endl;
	cout << "#col reads: " << num_col_reads_concept << " in concept" << endl;
	cout << "#col writes: " << num_col_writes_concept << " in concept" << endl;
	cout << "#multiply: " << num_multiply_concept << " in concept" << endl;
	cout << "#mem read bytes: " << detail::matrix_stats.get_read_bytes(true)
		<< endl;
	cout << "#mem write bytes: " << detail::matrix_stats.get_write_bytes(true)
		<< endl;
	cout << "#EM read bytes: " << detail::matrix_stats.get_read_bytes(false)
		<< endl;
	cout << "#EM write bytes: " << detail::matrix_stats.get_write_bytes(false)
		<< endl;
	cout << "#double float-point multiplies: "
		<< detail::matrix_stats.get_multiplies() << endl;
	cached_mats.clear();

	return res;
}

}

}
