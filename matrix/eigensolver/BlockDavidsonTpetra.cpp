// This example computes the eigenvalues of largest magnitude of an
// eigenvalue problem $A x = \lambda x$, using Anasazi's
// implementation of the Block Davidson method.

#include "sparse_matrix.h"

// Include header for block Davidson eigensolver
#include "AnasaziBlockDavidsonSolMgr.hpp"
// Include header to define eigenproblem Ax = \lambda*x
#include "AnasaziBasicEigenproblem.hpp"
#include "AnasaziOperator.hpp"
#include "AnasaziTpetraAdapter.hpp"
#include "Tpetra_Map.hpp"
#include "Tpetra_MultiVector.hpp"
#include "Tpetra_DefaultPlatform.hpp"

// ****************************************************************************
// BEGIN MAIN ROUTINE
// ****************************************************************************

using namespace fm;

std::atomic<size_t> iter_no;

typedef Tpetra::DefaultPlatform::DefaultPlatformType::NodeType  Node;
typedef Tpetra::MultiVector<double, size_t, size_t, Node> MV;

class FMTp_Operator//: public Anasazi::Operator<double>
{
	sparse_matrix::ptr mat;
public:
	FMTp_Operator(sparse_matrix::ptr mat) {
		this->mat = mat;
	}

	virtual void Apply(const MV &x, MV &y) const;

	size_t get_num_cols() const {
		return mat->get_num_cols();
	}

	size_t get_num_rows() const {
		return mat->get_num_rows();
	}
};

void FMTp_Operator::Apply(const MV &x, MV &y) const
{
	printf("SpMM %ld columns\n", x.getNumVectors());
	assert(x.getGlobalLength() == mat->get_num_cols());
	assert(y.getGlobalLength() == mat->get_num_rows());
	mem_vector::ptr in = mem_vector::create(mat->get_num_cols(),
			get_scalar_type<double>());
	mem_vector::ptr out = mem_vector::create(mat->get_num_rows(),
			get_scalar_type<double>());
	for (size_t i = 0; i < x.getNumVectors(); i++) {
		memcpy(in->get_raw_arr(), x.getData(i).getRawPtr(),
				in->get_length() * sizeof(double));
		mat->multiply<double>(*in, *out);
		memcpy(y.getDataNonConst(i).getRawPtr(), out->get_raw_arr(),
				out->get_length() * sizeof(double));
	}
}

namespace Anasazi
{

template<>
class OperatorTraits <double, MV, FMTp_Operator>
{
public:
	static void Apply(const FMTp_Operator &Op, const MV &x, MV &y) {
		iter_no++;
		Op.Apply(x, y);
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
	typedef FMTp_Operator OP;
	typedef Anasazi::MultiVecTraits<double, MV> MVT;
	typedef Tpetra::DefaultPlatform::DefaultPlatformType Platform;

	//
	// Set up the test problem.
	//
	config_map::ptr configs = config_map::create(conf_file);
	init_flash_matrix(configs);

	fg::FG_graph::ptr fg = fg::FG_graph::create(graph_file, index_file, configs);
	sparse_matrix::ptr m = sparse_matrix::create(fg);
	RCP<FMTp_Operator> A = rcp(new FMTp_Operator(m));

	// Set eigensolver parameters.
	const double tol = 1.0e-12; // convergence tolerance
	const int blockSize = nev + 1; // block size (number of eigenvectors processed at once)
	const int numBlocks = 3; // restart length
	const int maxRestarts = 100; // maximum number of restart cycles

	Platform &platform = Tpetra::DefaultPlatform::getDefaultPlatform();
	RCP<const Teuchos::Comm<int> > comm = platform.getComm();
	RCP<Tpetra::Map<size_t, size_t, Node> > Map
		= rcp (new Tpetra::Map<size_t, size_t, Node>(m->get_num_cols(), 0, comm));
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
	cout << iter_no << endl;

	return 0;
}
