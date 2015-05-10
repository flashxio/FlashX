// This example computes the eigenvalues of largest magnitude of an
// eigenvalue problem $A x = \lambda x$, using Anasazi's
// implementation of the Block Davidson method.

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

#include "safs_file.h"

#include "sparse_matrix.h"

// ****************************************************************************
// BEGIN MAIN ROUTINE
// ****************************************************************************

using namespace fm;

std::atomic<size_t> iter_no;

class FM_Operator//: public Anasazi::Operator<double>
{
	sparse_matrix::ptr mat;
public:
	FM_Operator(sparse_matrix::ptr mat) {
		this->mat = mat;
	}

	virtual void Apply(const FM_MultiVector<double>& x,
			FM_MultiVector<double>& y) const {
		printf("SpMM: y(%s) = A * x(%s)\n", y.get_name().c_str(), x.get_name().c_str());
		block_multi_vector::sparse_matrix_multiply<double>(*mat, *x.get_data(),
				*y.get_data());
		y.sync_fm2ep();

//		assert((size_t) x.GetGlobalLength() == mat->get_num_cols());
//		assert((size_t) y.GetGlobalLength() == mat->get_num_rows());
//		mem_vector::ptr in = mem_vector::create(mat->get_num_cols(),
//				get_scalar_type<double>());
//		mem_vector::ptr out = mem_vector::create(mat->get_num_rows(),
//				get_scalar_type<double>());
//		for (int i = 0; i < x.GetNumberVecs(); i++) {
//			memcpy(in->get_raw_arr(), x.get_ep_mv()[i], in->get_length() * sizeof(double));
//			out->reset_data();
//			mat->multiply<double>(*in, *out);
//			memcpy(y.get_ep_mv()[i], out->get_raw_arr(), out->get_length() * sizeof(double));
//		}
//		y.sync_ep2fm();
	}

	size_t get_num_cols() const {
		return mat->get_num_cols();
	}

	size_t get_num_rows() const {
		return mat->get_num_rows();
	}
};

namespace Anasazi
{

template<>
class OperatorTraits <double, FM_MultiVector<double>, FM_Operator>
{
public:
	/*! \brief This method takes the FM_MultiVector \c x and
	 * applies the FM_Operator \c Op to it resulting in the FM_MultiVector \c y.
	 */
	static void Apply ( const FM_Operator& Op,
			const FM_MultiVector<double>& x, FM_MultiVector<double>& y ) {
		iter_no++;
		Op.Apply(x,y);
	}

};

}

int
main (int argc, char *argv[])
{
	using Teuchos::RCP;
	using Teuchos::rcp;
	using std::cerr;
	using std::cout;
	using std::endl;

	if (argc < 6) {
		fprintf(stderr, "eigensolver conf_file matrix_file index_file nev solver\n");
		exit(1);
	}

	std::string conf_file = argv[1];
	std::string matrix_file = argv[2];
	std::string index_file = argv[3];
	int nev = atoi(argv[4]); // number of eigenvalues for which to solve
	std::string solver = argv[5];

	// Anasazi solvers have the following template parameters:
	//
	//   - Scalar: The type of dot product results.
	//   - MV: The type of (multi)vectors.
	//   - OP: The type of operators (functions from multivector to
	//     multivector).  A matrix (like Epetra_CrsMatrix) is an example
	//     of an operator; an Ifpack preconditioner is another example.
	//
	// Here, Scalar is double, MV is FM_MultiVector, and OP is FM_Operator.
	typedef FM_MultiVector<double> MV;
	typedef FM_Operator OP;
	typedef Anasazi::MultiVecTraits<double, MV> MVT;

	//
	// Set up the test problem.
	//
	config_map::ptr configs = config_map::create(conf_file);
	init_flash_matrix(configs);

	// Load index.
	SpM_2d_index::ptr index;
	safs::safs_file idx_f(safs::get_sys_RAID_conf(), index_file);
	if (idx_f.exist())
		index = SpM_2d_index::safs_load(index_file);
	else
		index = SpM_2d_index::load(index_file);

	// Load matrix.
	sparse_matrix::ptr mat;
	safs::safs_file mat_f(safs::get_sys_RAID_conf(), matrix_file);
	if (mat_f.exist())
		mat = sparse_matrix::create(index, safs::create_io_factory(
					matrix_file, safs::REMOTE_ACCESS));
	else
		mat = sparse_matrix::create(index,
				SpM_2d_storage::load(matrix_file, index));

	RCP<FM_Operator> A = rcp(new FM_Operator(mat));

	// Set eigensolver parameters.
	const double tol = 1.0e-8; // convergence tolerance
	const int numBlocks = 8; // restart length
	const int maxRestarts = 100; // maximum number of restart cycles
	const int maxIters = 500; // maximum number of iterations
	int blockSize = -1; // block size (number of eigenvectors processed at once)

	if (solver == "Davidson" || solver == "KrylovSchur") {
		blockSize = nev + 1;
	}
	else if (solver == "LOBPCG") {
		blockSize = 4;
	}

	// Create a set of initial vectors to start the eigensolver.
	// This needs to have the same number of columns as the block size.
	RCP<MV> ivec = rcp (new MV (A->get_num_cols(), blockSize, blockSize));
	ivec->Random ();

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
		cerr << "a wrong solver: " << solver << endl;
		exit(1);
	}

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
		MV tempAevec (A->get_num_rows(), sol.numVecs, evecs->get_block_size());
		T.putScalar (0.0);
		for (int i=0; i<sol.numVecs; ++i) {
			T(i,i) = evals[i].realpart;
		}
		A->Apply (*evecs, tempAevec);
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
	cout << "#iterations: " << iter_no << endl;
	cout << "#col writes: " << num_col_writes << endl;

	return 0;
}
