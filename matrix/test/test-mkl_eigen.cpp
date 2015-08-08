// This example computes the eigenvalues of largest magnitude of an
// eigenvalue problem $A x = \lambda x$, using Anasazi's
// implementation of the Block Davidson method.

#include "sparse_matrix.h"
#include "vertex.h"
#include "vertex_index.h"
#include "in_mem_storage.h"
#include "io_interface.h"

#include "crs_header.h"

#include <mkl.h>

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

#if 0
class MKL_MultiVector: public Anasazi::MultiVec<double>
{
#ifdef FM_VERIFY
	std::shared_ptr<Anasazi::EpetraMultiVec> ep_mat;
#endif
	std::vector<double> data;
	double *mat;
	size_t num_rows;
	size_t num_cols;

	MKL_MultiVector(size_t num_rows, size_t num_cols, double *mat) {
		this->num_rows = num_rows;
		this->num_cols = num_cols;
		this->mat = mat;
	}

public:
	MKL_MultiVector(size_t num_rows, size_t num_cols) {
		this->num_rows = num_rows;
		this->num_cols = num_cols;
		data.resize(num_rows * num_cols);
		mat = data.data();
#ifdef FM_VERIFY
		Epetra_SerialComm Comm;
		Teuchos::RCP<Epetra_Map> Map = Teuchos::rcp (
				new Epetra_Map((int) num_rows, 0, Comm));
		ep_mat = std::shared_ptr<Anasazi::EpetraMultiVec>(
				new Anasazi::EpetraMultiVec(*Map, num_cols));
#endif
	}

	void verify() const {
#ifdef FM_VERIFY
		int len = ep_mat->GlobalLength();
		assert((size_t) len == mat->get_num_rows());
		assert(mat->get_num_cols() == (size_t) ep_mat->NumVectors());
		size_t ep_col_idx = 0;
		for (size_t i = 0; i < mat->get_num_blocks(); i++) {
			dense_matrix::ptr block = mat->get_block(i);
			detail::local_matrix_store::const_ptr portion
				= block->get_data().get_portion(0);
			assert(block->get_num_rows() == portion->get_num_rows());
			assert(block->get_num_cols() == portion->get_num_cols());
			for (size_t j = 0; j < portion->get_num_cols(); j++) {
				for (size_t k = 0; k < portion->get_num_rows(); k++) {
					double v1 = (*ep_mat)[ep_col_idx][k];
					double v2 = portion->get<ScalarType>(k, j);
					if (v1 != v2)
						printf("v1: %g, v2: %g, diff: %g\n", v1, v2, v1 - v2);
					assert(abs(v1 - v2) == 0);
				}
				ep_col_idx++;
			}
		}
#endif
	}

	void Random () {
		std::mt19937_64 generator;
		std::uniform_real_distribution<double> dist;
		for (size_t i = 0; i < data.size(); i++)
			data[i] = dist(generator);
	}

	void sync_ep2fm() {
#ifdef FM_VERIFY
		assert(mat->get_num_cols() == (size_t) ep_mat->NumVectors());
		// This is a temporary solution.
		for (size_t i = 0; i < mat->get_num_blocks(); i++) {
			fm::detail::mem_matrix_store::ptr store
				= fm::detail::mem_matrix_store::create(mat->get_num_rows(),
						mat->get_block_size(), fm::matrix_layout_t::L_COL,
						mat->get_type(), mat->get_num_nodes());
			fm::detail::mem_col_matrix_store::ptr col_store
				= std::static_pointer_cast<fm::detail::mem_col_matrix_store>(store);
			for (size_t j = 0; j < mat->get_block_size(); j++) {
				size_t col_idx = i * mat->get_block_size() + j;
				memcpy(col_store->get_col(j), (*ep_mat)[col_idx],
						mat->get_num_rows() * mat->get_entry_size());
			}
			mat->set_block(i, fm::dense_matrix::create(col_store));
		}
#endif
	}

	void sync_fm2ep() {
#ifdef FM_VERIFY
		printf("There are %ld blocks and each block has %ld cols\n",
				mat->get_num_blocks(), mat->get_block(0)->get_num_cols());
		size_t ep_col_idx = 0;
		for (size_t i = 0; i < mat->get_num_blocks(); i++) {
			dense_matrix::ptr block = mat->get_block(i);
			detail::local_matrix_store::const_ptr portion
				= block->get_data().get_portion(0);
			assert(block->get_num_rows() == portion->get_num_rows());
			assert(block->get_num_cols() == portion->get_num_cols());
			for (size_t j = 0; j < portion->get_num_cols(); j++) {
				for (size_t k = 0; k < portion->get_num_rows(); k++)
					(*ep_mat)[ep_col_idx][k] = portion->get<ScalarType>(k, j);
				ep_col_idx++;
			}
		}
#endif
	}

	static void par_copy(void *dst, void *src, size_t num_bytes) {
		memcpy(dst, src, num_bytes);
	}

	static bool is_contig(const std::vector<int> &index) {
		if (index.size() == 1)
			return true;
		int prev = index[0];
		for (size_t i = 1; i < index.size(); i++) {
			if (prev + 1 != index[i])
				return false;
			prev++;
		}
		return true;
	}

	static void conv2row(const double *mat, size_t num_rows, size_t num_cols,
			double *res) {
		for (size_t i = 0; i < num_rows; i++)
			for (size_t j = 0; j < num_cols; j++)
				res[i * num_cols + j] = mat[j * num_rows + i];
	}

	static void conv2col(const double *mat, size_t num_rows, size_t num_cols,
			double *res) {
		for (size_t i = 0; i < num_rows; i++)
			for (size_t j = 0; j < num_cols; j++)
				res[j * num_rows + i] = mat[i * num_cols + j];
	}

	size_t get_num_rows() const {
		return num_rows;
	}

	size_t get_num_cols() const {
		return num_cols;
	}

	double *get_col(int off) {
		return mat + off * num_rows;
	}

	const double *get_col(int off) const {
		return mat + off * num_rows;
	}

	double *get_data() {
		return mat;
	}

	const double *get_data() const {
		return mat;
	}

	//! @name Creation methods

	/// \brief Create a new MultiVec with \c numvecs columns.
	/// \return Pointer to the new multivector with uninitialized values.
	virtual Anasazi::MultiVec<double> * Clone(const int numvecs) const {
		return new MKL_MultiVector(num_rows, numvecs);
	}

	/// \brief Create a new MultiVec and copy contents of \c *this into it (deep copy).
	/// \return Pointer to the new multivector	
	virtual Anasazi::MultiVec<double> * CloneCopy () const {
		MKL_MultiVector *ret = new MKL_MultiVector(num_rows, num_cols);
		par_copy(ret->mat, mat, num_rows * num_cols * sizeof(double));
		return ret;
#ifdef FM_VERIFY
		ret->ep_mat = std::shared_ptr<Anasazi::EpetraMultiVec>(
				dynamic_cast<Anasazi::EpetraMultiVec *>(ep_mat->CloneCopy()));
#endif
		verify();
		ret->verify();
		return ret;
	}

	/*! \brief Creates a new Anasazi::MultiVec and copies the selected contents of \c *this 
	  into the new vector (deep copy).  The copied 
	  vectors from \c *this are indicated by the \c index.size() indices in \c index.

	  \return Pointer to the new multivector	
	  */
	virtual Anasazi::MultiVec<double> * CloneCopy (
			const std::vector<int>& index) const {
		assert(is_contig(index));
		MKL_MultiVector *ret = new MKL_MultiVector(num_rows, index.size());
		par_copy(mat, ret->get_col(index[0]),
				num_rows * index.size() * sizeof(double));
#ifdef FM_VERIFY
		ret->ep_mat = std::shared_ptr<Anasazi::EpetraMultiVec>(
				dynamic_cast<Anasazi::EpetraMultiVec *>(ep_mat->CloneCopy(index)));
#endif
		verify();
		ret->verify();
		return ret;
	}

	/*! \brief Creates a new Anasazi::MultiVec that shares the selected contents of \c *this.
	  The index of the \c numvecs vectors shallow copied from \c *this are indicated by the
	  indices given in \c index.

	  \return Pointer to the new multivector	
	  */
	virtual Anasazi::MultiVec<double> * CloneViewNonConst (
			const std::vector<int>& index) {
		assert(is_contig(index));
#ifdef FM_VERIFY
		ret->ep_mat = std::shared_ptr<Anasazi::EpetraMultiVec>(
				dynamic_cast<Anasazi::EpetraMultiVec *>(ep_mat->CloneViewNonConst(index)));
#endif
//		verify();
//		ret->verify();
		return new MKL_MultiVector(num_rows, index.size(), get_col(index[0]));
	}

	/*! \brief Creates a new Anasazi::MultiVec that shares the selected contents of \c *this.
	  The index of the \c numvecs vectors shallow copied from \c *this are indicated by the
	  indices given in \c index.

	  \return Pointer to the new multivector	
	  */
	virtual const Anasazi::MultiVec<double> * CloneView (
			const std::vector<int>& index) const {
		assert(is_contig(index));
		const MKL_MultiVector *ret = new MKL_MultiVector(num_rows, index.size(),
				const_cast<double *>(get_col(index[0])));
#ifdef FM_VERIFY
		ret->ep_mat = std::shared_ptr<Anasazi::EpetraMultiVec>(
				dynamic_cast<Anasazi::EpetraMultiVec *>(ep_mat->CloneViewNonConst(index)));
#endif
//		verify();
		ret->verify();
		return ret;
	}

	//! @name Dimension information methods	

	//! The number of rows in the multivector.
	ANASAZI_DEPRECATED virtual int GetVecLength () const {
		return num_rows;
	}

	//! The number of rows in the multivector.
	//! \note This method supersedes GetVecLength, which will be deprecated.
	virtual ptrdiff_t GetGlobalLength () const {
		return num_rows;
	}

	//! The number of vectors (i.e., columns) in the multivector.
	virtual int GetNumberVecs () const {
		return num_cols;
	}

	//! @name Update methods

	//! Update \c *this with \c alpha * \c A * \c B + \c beta * (\c *this).
	virtual void MvTimesMatAddMv (double alpha, 
			const Anasazi::MultiVec<double>& A, 
			const Teuchos::SerialDenseMatrix<int, double>& B, double beta) {
		sync_fm2ep();
		const MKL_MultiVector &fm_A = dynamic_cast<const MKL_MultiVector &>(A);

		std::vector<double> Bmat(B.numRows() * B.numCols());
		for (int j = 0; j < B.numCols(); j++)
			for (int i = 0; i < B.numRows(); i++)
				Bmat[B.numRows() * j + i] = B(i, j);
		cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
				fm_A.num_rows, B.numCols(), fm_A.num_cols, alpha, fm_A.mat,
				fm_A.num_rows, Bmat.data(), B.numRows(), beta, mat,
				num_rows);
#ifdef FM_VERIFY
		this->ep_mat->MvTimesMatAddMv(alpha, *fm_A.ep_mat, B, beta);
#endif
		fm_A.verify();
		this->verify();
	}

	//! Replace \c *this with \c alpha * \c A + \c beta * \c B.
	virtual void MvAddMv (double alpha, const Anasazi::MultiVec<double>& A,
			double beta, const Anasazi::MultiVec<double>& B ) {
		sync_fm2ep();
		const MKL_MultiVector &fm_A = dynamic_cast<const MKL_MultiVector &>(A);
		const MKL_MultiVector &fm_B = dynamic_cast<const MKL_MultiVector &>(B);

		if (alpha == 1 && beta == 0)
			par_copy(mat, fm_A.mat, num_rows * num_cols * sizeof(double));
		else if (alpha == 0 && beta == 1)
			par_copy(mat, fm_B.mat, num_rows * num_cols * sizeof(double));
		else {
			long double lalpha = alpha;
			long double lbeta = beta;
			for (size_t i = 0; i < num_rows * num_cols; i++)
				mat[i] = fm_A.mat[i] * lalpha + fm_B.mat[i] * lbeta;
		}
#ifdef FM_VERIFY
		this->ep_mat->MvAddMv(alpha, *fm_A.ep_mat, beta, *fm_B.ep_mat);
#endif
		fm_A.verify();
		fm_B.verify();
		this->verify();
	}

	//! Scale each element of the vectors in \c *this with \c alpha.
	virtual void MvScale (double alpha) {
		sync_fm2ep();
		long double lalpha = alpha;
		for (size_t i = 0; i < num_rows * num_cols; i++)
			mat[i] *= lalpha;
#ifdef FM_VERIFY
		ep_mat->MvScale(alpha);
#endif
		this->verify();
	}

	//! Scale each element of the <tt>i</tt>-th vector in \c *this with <tt>alpha[i]</tt>.
	virtual void MvScale ( const std::vector<double>& alpha ) {
		sync_fm2ep();
		std::vector<long double> lalpha(alpha.begin(), alpha.end());
		for (size_t i = 0; i < num_cols; i++)
			for (size_t j = 0; j < num_rows; j++)
				get_col(i)[j] *= lalpha[i];
#ifdef FM_VERIFY
		ep_mat->MvScale(alpha);
#endif
		this->verify();
	}

	/*! \brief Compute a dense matrix \c B through the matrix-matrix multiply 
	  \c alpha * \c A^T * (\c *this).
	  */
	virtual void MvTransMv (double alpha, const Anasazi::MultiVec<double>& A,
			Teuchos::SerialDenseMatrix<int, double>& B
#ifdef HAVE_ANASAZI_EXPERIMENTAL
			, ConjType conj = Anasazi::CONJ
#endif				 
			) const {
		const MKL_MultiVector &fm_A = dynamic_cast<const MKL_MultiVector &>(A);

		assert((size_t) B.numRows() == fm_A.num_cols);
		assert((size_t) B.numCols() == this->num_cols);
		assert(fm_A.num_rows == this->num_rows);
		this->verify();

		std::vector<double> tA(fm_A.num_rows * fm_A.num_cols);
		conv2row(fm_A.mat, fm_A.num_rows, fm_A.num_cols, tA.data());
		std::vector<double> res(B.numRows() * B.numCols());
		cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
				fm_A.num_cols, num_cols, fm_A.num_rows, alpha, tA.data(),
				fm_A.num_cols, mat, num_rows, 0, res.data(), B.numRows());
#if 0
		ep_mat->MvTransMv(alpha, *fm_A.ep_mat, B);
		for (int i = 0; i < B.numRows(); i++)
			for (int j = 0; j < B.numCols(); j++) {
				printf("%g\n", B(i, j) - (double) (res->get<ScalarType>(i, j) * lalpha));
				assert(B(i, j) == (double) (res->get<ScalarType>(i, j) * lalpha));
			}
#endif
		for (int i = 0; i < B.numRows(); i++)
			for (int j = 0; j < B.numCols(); j++)
				B(i, j) = res[B.numRows() * j + i];
	}

	/// \brief Compute the dot product of each column of *this with the corresponding column of A.
	///
	/// Compute a vector \c b whose entries are the individual
	/// dot-products.  That is, <tt>b[i] = A[i]^H * (*this)[i]</tt>
	/// where <tt>A[i]</tt> is the i-th column of A.
	virtual void MvDot ( const Anasazi::MultiVec<double>& A,
			std::vector<double> & b 
#ifdef HAVE_ANASAZI_EXPERIMENTAL
			, ConjType conj = Anasazi::CONJ
#endif
			) const {
		// TODO
		assert(0);
	}

	//! @name Norm method

	/// \brief Compute the 2-norm of each vector in \c *this.  
	///
	/// \param normvec [out] On output, normvec[i] holds the 2-norm of the
	///   \c i-th vector of \c *this.
	virtual void MvNorm (
			std::vector<typename Teuchos::ScalarTraits<double>::magnitudeType> & normvec) const {
		verify();
		for (size_t i = 0; i < normvec.size(); i++) {
			long double sq_sum = 0;
			for (size_t j = 0; j < num_rows; j++)
				sq_sum += ((long double) get_col(i)[j]) * get_col(i)[j];
			normvec[i] = sqrt(sq_sum);
		}
#ifdef FM_VERIFY
		std::vector<typename Teuchos::ScalarTraits<double>::magnitudeType> normvec1(
				normvec.size()); 
		ep_mat->MvNorm(normvec1);
		for (size_t i = 0; i < normvec.size(); i++) {
			printf("%ld: %g\n", i, normvec[i] - normvec1[i]);
			assert(normvec[i] == normvec1[i]);
		}
#endif
	}

	//! @name Initialization methods

	/// \brief Copy the vectors in \c A to a set of vectors in \c *this.  
	///
	/// The \c numvecs vectors in \c A are copied to a subset of vectors
	/// in \c *this indicated by the indices given in \c index.
	virtual void SetBlock (const Anasazi::MultiVec<double>& A,
			const std::vector<int>& index) {
		sync_fm2ep();
		assert(is_contig(index));
		assert((size_t) A.GetNumberVecs() == index.size());
		const MKL_MultiVector &fm_A = dynamic_cast<const MKL_MultiVector &>(A);
		par_copy(get_col(index[0]), fm_A.mat,
				fm_A.num_rows * fm_A.num_cols * sizeof(double));
#ifdef FM_VERIFY
		this->ep_mat->SetBlock(*fm_A.ep_mat, index);
#endif
		fm_A.verify();
		this->verify();
	}

	//! Fill all the vectors in \c *this with random numbers.
	virtual void MvRandom () {
		Random();
		sync_fm2ep();
		this->verify();
	}

	//! Replace each element of the vectors in \c *this with \c alpha.
	virtual void MvInit (double alpha) {
		// TODO
		assert(0);
	}

	//! @name Print method

	//! Print \c *this multivector to the \c os output stream.
	virtual void MvPrint ( std::ostream& os ) const {
		// TODO
		assert(0);
	}

#ifdef HAVE_ANASAZI_TSQR
	//! @name TSQR-related methods

	/// \brief Compute the QR factorization *this = QR, using TSQR.
	///
	/// The *this multivector on input is the multivector A to factor.
	/// It is overwritten with garbage on output.
	///
	/// \param Q [out] On input: a multivector with the same number of
	///   rows and columns as A (the *this multivector).  Its contents
	///   are overwritten on output with the (explicitly stored) Q
	///   factor in the QR factorization of A.
	///
	/// \param R [out] On output: the R factor in the QR factorization
	///   of the (input) multivector A.
	///
	/// \param forceNonnegativeDiagonal [in] If true, then (if
	///   necessary) do extra work (modifying both the Q and R
	///   factors) in order to force the R factor to have a
	///   nonnegative diagonal.
	///
	/// For syntax's sake, we provide a default implementation of this
	/// method that throws std::logic_error.  You should implement this
	/// method if you intend to use TsqrOrthoManager or
	/// TsqrMatOrthoManager with your subclass of MultiVec.
	virtual void factorExplicit (Anasazi::MultiVec<ScalarType>& Q, 
			Teuchos::SerialDenseMatrix<int, ScalarType>& R,
			const bool forceNonnegativeDiagonal=false) {
		TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "The Anasazi::MultiVec<" 
				<< Teuchos::TypeNameTraits<ScalarType>::name() << "> subclass which you "
				"are using does not implement the TSQR-related method factorExplicit().");
	}

	/// \brief Use result of factorExplicit() to compute rank-revealing decomposition.
	///
	/// When calling this method, the *this multivector should be the Q
	/// factor output of factorExplicit().  Using that Q factor and the
	/// R factor from factorExplicit(), compute the singular value
	/// decomposition (SVD) of R (\f$R = U \Sigma V^*\f$).  If R is full
	/// rank (with respect to the given relative tolerance tol), don't
	/// change Q (= *this) or R.  Otherwise, compute \f$Q := Q \cdot
	/// U\f$ and \f$R := \Sigma V^*\f$ in place (the latter may be no
	/// longer upper triangular).
	///
	/// The *this multivector on input must be the explicit Q factor
	/// output of a previous call to factorExplicit().  On output: On
	/// output: If R is of full numerical rank with respect to the
	/// tolerance tol, Q is unmodified.  Otherwise, Q is updated so that
	/// the first rank columns of Q are a basis for the column space of
	/// A (the original matrix whose QR factorization was computed by
	/// factorExplicit()).  The remaining columns of Q are a basis for
	/// the null space of A.
	///
	/// \param R [in/out] On input: N by N upper triangular matrix with
	///   leading dimension LDR >= N.  On output: if input is full rank,
	///   R is unchanged on output.  Otherwise, if \f$R = U \Sigma
	///   V^*\f$ is the SVD of R, on output R is overwritten with
	///   \f$\Sigma \cdot V^*\f$.  This is also an N by N matrix, but
	///   may not necessarily be upper triangular.
	///
	/// \param tol [in] Relative tolerance for computing the numerical
	///   rank of the matrix R.
	///
	/// For syntax's sake, we provide a default implementation of this
	/// method that throws std::logic_error.  You should implement this
	/// method if you intend to use TsqrOrthoManager or
	/// TsqrMatOrthoManager with your subclass of MultiVec.
	virtual int revealRank (Teuchos::SerialDenseMatrix<int, ScalarType>& R,
			const typename Teuchos::ScalarTraits<ScalarType>::magnitudeType& tol) {
		TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "The Anasazi::MultiVec<" 
				<< Teuchos::TypeNameTraits<ScalarType>::name() << "> subclass which you "
				"are using does not implement the TSQR-related method revealRank().");
	}

#endif // HAVE_ANASAZI_TSQR
};

namespace Anasazi
{

/// \brief Specialization of MultiVecTraits for Belos::MultiVec.
///
/// Anasazi interfaces to every multivector implementation through a
/// specialization of MultiVecTraits.  Thus, we provide a
/// specialization of MultiVecTraits for the MultiVec run-time
/// polymorphic interface above.
///
/// \tparam ScalarType The type of entries in the multivector; the
///   template parameter of MultiVec.
template<>
class MultiVecTraits<double, MKL_MultiVector> {
public:
	//! @name Creation methods
	//@{ 

	/// \brief Create a new empty \c MultiVec containing \c numvecs columns.
	/// \return Reference-counted pointer to the new \c MultiVec.
	static Teuchos::RCP<MKL_MultiVector> Clone (const MKL_MultiVector& mv,
			const int numvecs) {
		return Teuchos::rcp ((MKL_MultiVector *) mv.Clone(numvecs)); 
	}

	/*!
	 * \brief Creates a new \c Anasazi::MultiVec and copies contents of \c mv
	 * into the new vector (deep copy).
	  \return Reference-counted pointer to the new \c Anasazi::MultiVec.
	  */
	static Teuchos::RCP<MKL_MultiVector> CloneCopy(const MKL_MultiVector& mv) {
		return Teuchos::rcp((MKL_MultiVector *) mv.CloneCopy());
	}

	/*! \brief Creates a new \c Anasazi::MultiVec and copies the selected
	 * contents of \c mv into the new vector (deep copy).  

	  The copied vectors from \c mv are indicated by the \c index.size() indices in \c index.      
	  \return Reference-counted pointer to the new \c Anasazi::MultiVec.
	  */
	static Teuchos::RCP<MKL_MultiVector> CloneCopy(
			const MKL_MultiVector& mv, const std::vector<int>& index) {
		return Teuchos::rcp((MKL_MultiVector *) mv.CloneCopy(index));
	}

	/*! \brief Creates a new \c Anasazi::MultiVec that shares the selected
	 * contents of \c mv (shallow copy).

	  The index of the \c numvecs vectors shallow copied from \c mv are indicated
	  by the indices given in \c index.
	  \return Reference-counted pointer to the new \c Anasazi::MultiVec.
	  */    
	static Teuchos::RCP<MKL_MultiVector> CloneViewNonConst(
			MKL_MultiVector& mv, const std::vector<int>& index) {
		return Teuchos::rcp((MKL_MultiVector *) mv.CloneViewNonConst(index));
	}

	/*! \brief Creates a new const \c Anasazi::MultiVec that shares
	 * the selected contents of \c mv (shallow copy).

	  The index of the \c numvecs vectors shallow copied from \c mv are
	  indicated by the indices given in \c index.
	  \return Reference-counted pointer to the new const \c Anasazi::MultiVec.
	  */      
	static Teuchos::RCP<const MKL_MultiVector> CloneView(
			const MKL_MultiVector& mv, const std::vector<int>& index) {
		return Teuchos::rcp((MKL_MultiVector *) mv.CloneView(index));
	}

	//@}

	//! @name Attribute methods
	//@{ 

	//! Obtain the vector length of \c mv.
	ANASAZI_DEPRECATED static int GetVecLength(const MKL_MultiVector& mv) {
		return mv.GetGlobalLength();
	}

	//! Obtain the number of vectors in \c mv
	static int GetNumberVecs(const MKL_MultiVector& mv) {
		return mv.GetNumberVecs();
	}

	//! Obtain the vector length of \c mv.
	//! \note This method supersedes GetVecLength, which will be deprecated.
	static ptrdiff_t GetGlobalLength(const MKL_MultiVector& mv) {
		return mv.GetGlobalLength();
	}

	//@}

	//! @name Update methods
	//@{ 

	/*! \brief Update \c mv with \f$ \alpha AB + \beta mv \f$.
	*/
	static void MvTimesMatAddMv(double alpha, const MKL_MultiVector& A, 
			const Teuchos::SerialDenseMatrix<int, double>& B, 
			double beta, MKL_MultiVector& mv)
	{ mv.MvTimesMatAddMv(alpha, A, B, beta); }

	/*! \brief Replace \c mv with \f$\alpha A + \beta B\f$.
	*/
	static void MvAddMv(double alpha, const MKL_MultiVector& A,
			double beta, const MKL_MultiVector& B,
			MKL_MultiVector& mv)
	{ mv.MvAddMv(alpha, A, beta, B); }

	/*! \brief Compute a dense matrix \c B through the matrix-matrix multiply \f$ \alpha A^Tmv \f$.
	*/
	static void MvTransMv(double alpha, const MKL_MultiVector& A,
			const MKL_MultiVector& mv, Teuchos::SerialDenseMatrix<int, double>& B
#ifdef HAVE_ANASAZI_EXPERIMENTAL
			, ConjType conj = Anasazi::CONJ
#endif
			)
	{ mv.MvTransMv(alpha, A, B
#ifdef HAVE_ANASAZI_EXPERIMENTAL
			, conj
#endif
			); }

	/*! \brief Compute a vector \c b where the components are the individual dot-products of the \c i-th columns of \c A and \c mv, i.e.\f$b[i] = A[i]^H mv[i]\f$.
	*/
	static void MvDot( const MKL_MultiVector& mv,
			const MKL_MultiVector& A, std::vector<double> &b
#ifdef HAVE_ANASAZI_EXPERIMENTAL
			, ConjType conj = Anasazi::CONJ
#endif
			)
	{ mv.MvDot( A, b
#ifdef HAVE_ANASAZI_EXPERIMENTAL
			, conj
#endif
			); }

	//! Scale each element of the vectors in \c *this with \c alpha.
	static void MvScale (MKL_MultiVector& mv, double alpha)
	{ mv.MvScale( alpha ); }

	//! Scale each element of the \c i-th vector in \c *this with \c alpha[i].
	static void MvScale (MKL_MultiVector& mv, const std::vector<double>& alpha)
	{ mv.MvScale( alpha ); }

	//@}
	//! @name Norm method
	//@{ 

	/*! \brief Compute the 2-norm of each individual vector of \c mv.  
	  Upon return, \c normvec[i] holds the value of \f$||mv_i||_2\f$, the \c i-th column of \c mv.
	  */
	static void MvNorm(const MKL_MultiVector& mv,
			std::vector<typename Teuchos::ScalarTraits<double>::magnitudeType> & normvec)
	{ mv.MvNorm(normvec); }

	//@}
	//! @name Initialization methods
	//@{ 
	/*! \brief Copy the vectors in \c A to a set of vectors in \c mv indicated by the indices given in \c index.

	  The \c numvecs vectors in \c A are copied to a subset of vectors in \c mv indicated by the indices given in \c index,
	  i.e.<tt> mv[index[i]] = A[i]</tt>.
	  */
	static void SetBlock(const MKL_MultiVector& A, const std::vector<int>& index,
			MKL_MultiVector& mv )
	{ mv.SetBlock(A, index); }

	/*! \brief Replace the vectors in \c mv with random vectors.
	*/
	static void MvRandom(MKL_MultiVector& mv)
	{ mv.MvRandom(); }

	/*! \brief Replace each element of the vectors in \c mv with \c alpha.
	*/
	static void MvInit(MKL_MultiVector& mv,
			double alpha = Teuchos::ScalarTraits<double>::zero())
	{ mv.MvInit(alpha); }

	//@}
	//! @name Print method
	//@{ 

	//! Print the \c mv multi-vector to the \c os output stream.
	static void MvPrint(const MKL_MultiVector& mv, std::ostream& os)
	{ mv.MvPrint(os); }

	//@}

#ifdef HAVE_ANASAZI_TSQR
	/// \typedef tsqr_adaptor_type
	/// \brief TSQR adapter for MultiVec.
	///
	/// Our TSQR adapter for MultiVec calls MultiVec's virtual
	/// methods.  If you want to use TSQR with your MultiVec subclass,
	/// you must implement these methods yourself, as the default
	/// implementations throw std::logic_error.
	typedef details::MultiVecTsqrAdapter<double> tsqr_adaptor_type;
#endif // HAVE_ANASAZI_TSQR
};

}
#endif

std::atomic<size_t> iter_no;

typedef int local_ordinal_type;
typedef int global_ordinal_type;

typedef Tpetra::DefaultPlatform::DefaultPlatformType::NodeType  Node;
typedef Tpetra::MultiVector<double, global_ordinal_type, global_ordinal_type, Node> MV;

size_t edge_data_size = 0;

static void read_index(const std::string &crs_file, std::vector<MKL_INT> &row_idxs,
		std::vector<MKL_INT> &col_vec)
{
	FILE *f = fopen(crs_file.c_str(), "r");
	if (f == NULL) {
		fprintf(stderr, "can't open %s: %s\n", crs_file.c_str(), strerror(errno));
		exit(1);
	}

	crs_header header;
	size_t ret = fread(&header, sizeof(header), 1, f);
	if (ret == 0) {
		fprintf(stderr, "can't read header: %s\n", strerror(errno));
		exit(1);
	}

	size_t nrowA = header.get_num_rows();
	size_t ncolA = header.get_num_cols();
	size_t nnz = header.get_num_non_zeros();
	printf("There are %ld rows, %ld cols and %ld nnz\n", nrowA, ncolA, nnz);

	std::vector<crs_idx_t> row_idxs64(nrowA + 1);
	std::vector<crs_idx_t> col_vec64(nnz);
	printf("read CRS\n");
	ret = fread(row_idxs64.data(), row_idxs64.size() * sizeof(row_idxs64[0]),
			1, f);
	if (ret == 0) {
		fprintf(stderr, "can't read row idxs: %s\n", strerror(errno));
		exit(1);
	}

	ret = fread(col_vec64.data(), col_vec64.size() * sizeof(col_vec64[0]), 1, f);
	if (ret == 0) {
		fprintf(stderr, "can't read col idxs: %s\n", strerror(errno));
		exit(1);
	}

	for (size_t i = 0; i < col_vec64.size(); i++)
		assert(col_vec64[i] < ncolA);

	assert(row_idxs.empty());
	row_idxs.insert(row_idxs.end(), row_idxs64.begin(), row_idxs64.end());
	assert(col_vec.empty());
	col_vec.insert(col_vec.end(), col_vec64.begin(), col_vec64.end());
}

static void conv2row(const MV &src, double *res)
{
	size_t col_length = src.getGlobalLength();
	for (size_t j = 0; j < src.getNumVectors(); j++) {
		Teuchos::ArrayRCP<const double> col = src.getData(j);
		for (size_t i = 0; i < col_length; i++)
			res[i * src.getNumVectors() + j] = col[i];
	}
}

void conv2col(const double *mat, MV &res)
{
	size_t col_length = res.getGlobalLength();
	size_t num_cols = res.getNumVectors();
	for (size_t j = 0; j < num_cols; j++) {
		Teuchos::ArrayRCP<double> col = res.getDataNonConst(j);
		for (size_t i = 0; i < col_length; i++)
			col[i] = mat[i * num_cols + j];
	}
}

class spm_function
{
public:
	typedef std::shared_ptr<const spm_function> const_ptr;
	virtual void run(const MV &in, MV &out) const = 0;
	virtual size_t get_num_rows() const = 0;
	virtual size_t get_num_cols() const = 0;
};

class eigen_function: public spm_function
{
	std::vector<MKL_INT> row_idxs;
	std::vector<MKL_INT> col_vec;
	std::vector<double> val_vec;
public:
	eigen_function(const std::string &crs_file);

	virtual void run(const MV &in, MV &out) const;
	virtual size_t get_num_rows() const {
		return row_idxs.size() - 1;
	}
	virtual size_t get_num_cols() const {
		return row_idxs.size() - 1;
	}
};

class SVD_function: public spm_function
{
	std::vector<MKL_INT> row_idxs;
	std::vector<MKL_INT> col_vec;
	std::vector<MKL_INT> t_row_idxs;
	std::vector<MKL_INT> t_col_vec;
	// Both matrices share the same value array because both matrices have
	// the same number of nnz and all nnz have the same value.
	std::vector<double> val_vec;
public:
	SVD_function(const std::string &crs_file, const std::string t_crs_file);

	virtual void run(const MV &in, MV &out) const;
	virtual size_t get_num_rows() const {
		return row_idxs.size() - 1;
	}
	virtual size_t get_num_cols() const {
		return row_idxs.size() - 1;
	}
};

eigen_function::eigen_function(const std::string &crs_file)
{
	read_index(crs_file, row_idxs, col_vec);
	size_t nnz = col_vec.size();
	val_vec.resize(nnz);
	printf("Create value vector of %ld entries for the adj matrix\n", val_vec.size());
#pragma omp parallel for
	for (size_t i = 0; i < val_vec.size(); i++)
		val_vec[i] = 1;
}

void eigen_function::run(const MV &in, MV &out) const
{
	std::vector<double> in_data(in.getGlobalLength() * in.getNumVectors());
	std::vector<double> out_data(out.getGlobalLength() * out.getNumVectors());
	MKL_INT ncolC = out.getNumVectors();
	MKL_INT nrowA = get_num_rows();
	MKL_INT ncolA = get_num_cols();
	double alpha = 1;
	double beta = 0;
	// Copy data from `in' (in column major) to `in_data' (in row major).
	conv2row(in, in_data.data());
	mkl_dcsrmm("N", &nrowA, &ncolC, &ncolA, &alpha, "G  C", val_vec.data(),
			col_vec.data(), row_idxs.data(), row_idxs.data() + 1, in_data.data(),
			&ncolC, &beta, out_data.data(), &ncolC);
	// Copy data from `out_data' (in row major) to `out' (in column major).
	conv2col(out_data.data(), out);
}

SVD_function::SVD_function(const std::string &crs_file,
		const std::string t_crs_file)
{
	read_index(crs_file, row_idxs, col_vec);
	read_index(t_crs_file, t_row_idxs, t_col_vec);
	size_t nnz = col_vec.size();
	val_vec.resize(nnz);
	printf("Create value vector of %ld entries for the adj matrix\n", val_vec.size());
#pragma omp parallel for
	for (size_t i = 0; i < val_vec.size(); i++)
		val_vec[i] = 1;
}

void SVD_function::run(const MV &in, MV &out) const
{
	std::vector<double> in_data(in.getGlobalLength() * in.getNumVectors());
	std::vector<double> tmp_data(out.getGlobalLength() * out.getNumVectors());
	std::vector<double> out_data(out.getGlobalLength() * out.getNumVectors());
	MKL_INT ncolC = out.getNumVectors();
	MKL_INT nrowA = get_num_rows();
	MKL_INT ncolA = get_num_cols();
	double alpha = 1;
	double beta = 0;
	// Copy data from `in' (in column major) to `in_data' (in row major).
	conv2row(in, in_data.data());
	mkl_dcsrmm("N", &nrowA, &ncolC, &ncolA, &alpha, "G  C", val_vec.data(),
			col_vec.data(), row_idxs.data(), row_idxs.data() + 1, in_data.data(),
			&ncolC, &beta, tmp_data.data(), &ncolC);

	mkl_dcsrmm("N", &nrowA, &ncolC, &ncolA, &alpha, "G  C", val_vec.data(),
			t_col_vec.data(), t_row_idxs.data(), t_row_idxs.data() + 1,
			tmp_data.data(), &ncolC, &beta, out_data.data(), &ncolC);
	// Copy data from `out_data' (in row major) to `out' (in column major).
	conv2col(out_data.data(), out);
}

namespace Anasazi
{

template<>
class OperatorTraits <double, MV, spm_function>
{
public:
	static void Apply(const spm_function& Op, const MV& x, MV& y) {
		printf("run sparse matrix multiplication\n");
		Op.run(x, y);
	}
};

}

typedef spm_function OP;
typedef Anasazi::OperatorTraits<double, MV, OP> OPT;
typedef Tpetra::Map<local_ordinal_type, global_ordinal_type, Node> map_type;

RCP<map_type> Map;

void compute_eigen(spm_function *func, int nev, const std::string &solver,
		int blockSize, int numBlocks, double tol)
{
	// Here, Scalar is double, MV is Tpetra::MultiVector, and OP is FMTp_Operator.
	typedef Anasazi::MultiVecTraits<double, MV> MVT;

	// Set eigensolver parameters.
	const int maxRestarts = 100; // maximum number of restart cycles
	const int maxIters = 500;

	printf("solver: %s\n", solver.c_str());
	printf("block size: %d\n", blockSize);
	printf("#blocks: %d\n", numBlocks);
	printf("tol: %g\n", tol);

	RCP<OP> A = rcp(func);
	// Create a set of initial vectors to start the eigensolver.
	// This needs to have the same number of columns as the block size.
	RCP<MV> ivec = rcp (new MV (Map, (size_t) blockSize));
	ivec->randomize();

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
}

void print_usage()
{
	fprintf(stderr, "eigensolver csr_file nev [options]\n");
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
	std::string type = "eigen";
	std::string solver = "LOBPCG";
	double tol = 1e-8;
	while ((opt = getopt(argc, argv, "b:B:n:N:s:t:T:")) != -1) {
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
			case 'T':
				type = optarg;
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
	if ((argc < 2 && type != "SVD") || (argc < 3 && type == "SVD")) {
		print_usage();
		exit(1);
	}

	std::string csr_file = argv[0];
	std::string t_csr_file;
	int nev;
	if (type == "SVD") {
		t_csr_file = argv[1];
		nev = atoi(argv[2]);
	}
	else
		nev = atoi(argv[1]); // number of eigenvalues for which to solve;

	// Anasazi solvers have the following template parameters:
	//
	//   - Scalar: The type of dot product results.
	//   - MV: The type of (multi)vectors.
	//   - OP: The type of operators (functions from multivector to
	//     multivector).  A matrix (like Epetra_CrsMatrix) is an example
	//     of an operator; an Ifpack preconditioner is another example.
	//
	typedef Tpetra::DefaultPlatform::DefaultPlatformType Platform;

	//
	// Set up the test problem.
	//
	spm_function *A;
	if (type == "eigen")
		A = new eigen_function(csr_file);
	else if (type == "SVD")
		A = new SVD_function(csr_file, t_csr_file);
	else {
		fprintf(stderr, "wrong type: %s\n", type.c_str());
		exit(1);
	}
	assert(A->get_num_rows() == A->get_num_cols());

	Platform &platform = Tpetra::DefaultPlatform::getDefaultPlatform();
	RCP<const Teuchos::Comm<int> > comm = platform.getComm();
	Map = rcp (new map_type(A->get_num_rows(), 0, comm));
	for (size_t blockSize = blockSizeStart; blockSize <= blockSizeEnd;
			blockSize *= 2) {
		for (size_t numBlocks = numBlocksStart; numBlocks <= numBlocksEnd;
				numBlocks *= 2)
			compute_eigen(A, nev, solver, blockSize, numBlocks, tol);
	}

	return 0;
}
