#ifndef __FM_MULTIVECTOR_H__
#define __FM_MULTIVECTOR_H__

/*
 * Copyright 2014 Open Connectome Project (http://openconnecto.me)
 * Written by Da Zheng (zhengda1936@gmail.com)
 *
 * This file is part of FlashMatrix.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#include <stdlib.h>

#include "AnasaziMultiVec.hpp"

#include "NUMA_vector.h"
#include "mem_dense_matrix.h"
#include "matrix_config.h"
#include "generic_type.h"

#include "Epetra_MultiVector.h"
#include "Epetra_SerialComm.h"
#include "Epetra_Map.h"
#include "AnasaziEpetraAdapter.hpp"

static int MV_id;

template<class ScalarType>
class FM_MultiVector: public Anasazi::MultiVec<ScalarType>
{
	class rand_init: public fm::type_set_operate<ScalarType>
	{
	public:
		virtual void set(ScalarType *arr, size_t num_eles, off_t row_idx,
				off_t col_idx) const {
			for (size_t i = 0; i < num_eles; i++)
				arr[i] = random();
		}
	};

	std::string name;
	fm::mem_col_dense_matrix::ptr mat;
	Teuchos::RCP<Epetra_Map> Map;
	std::shared_ptr<Anasazi::EpetraMultiVec> ep_mat;

	static void copy_ep2fm(const Anasazi::EpetraMultiVec &ep_mat,
			fm::mem_col_dense_matrix &mat) {
		size_t len = ep_mat.GlobalLength();
		assert((size_t) ep_mat.NumVectors() == mat.get_num_cols());
		assert(len == mat.get_num_rows());
		for (int i = 0; i < ep_mat.NumVectors(); i++)
			memcpy(mat.get_col(i), ep_mat[i], len * sizeof(double));
	}

	static void verify(const fm::mem_col_dense_matrix &mat1,
			const fm::mem_col_dense_matrix &mat2) {
		assert(mat1.get_num_cols() == mat1.get_num_cols());
		assert(mat1.get_num_rows() == mat1.get_num_rows());
		for (size_t i = 0; i < mat1.get_num_cols(); i++) {
			assert(memcmp(mat1.get_col(i), mat2.get_col(i),
						mat1.get_num_rows() * mat1.get_entry_size()) == 0);
		}
	}

	void verify() const {
//		int len = ep_mat->GlobalLength();
//		assert((size_t) ep_mat->NumVectors() == mat->get_num_cols());
//		assert((size_t) len == mat->get_num_rows());
//		double max = 0;
//		for (int i = 0; i < ep_mat->NumVectors(); i++) {
//			for (int j = 0; j < len; j++) {
//				double v1 = (*ep_mat)[i][j];
//				double v2 = *(double *) mat->get(j, i);
//				if (abs(v1 - v2) > 0) {
//					printf("%g\n", *(double *) mat->get(0, 0));
//					fprintf(stderr, "%s: %d, %d: %g\n", name.c_str(), i, j, abs(v1 - v2));
//				}
//				assert(abs(v1 - v2) == 0);
//				max = std::max(max, abs(v1 - v2));
//			}
//		}
//		printf("max diff: %g\n", max);
	}

	FM_MultiVector(const std::string &extra) {
		char name_buf[128];
		snprintf(name_buf, sizeof(name_buf), "MV-%d", MV_id++);
		this->name = std::string(name_buf) + " " + extra;
	}
public:
	FM_MultiVector(size_t num_rows, size_t num_cols) {
		mat = fm::mem_col_dense_matrix::create(num_rows, num_cols,
				fm::get_scalar_type<ScalarType>());

		Epetra_SerialComm Comm;
		Map = Teuchos::rcp (new Epetra_Map((int) num_rows, 0, Comm));
		ep_mat = std::shared_ptr<Anasazi::EpetraMultiVec>(
				new Anasazi::EpetraMultiVec(*Map, num_cols));

		char name_buf[128];
		snprintf(name_buf, sizeof(name_buf), "MV-%d", MV_id++);
		this->name = name_buf;
	}

	fm::mem_col_dense_matrix &get_data() {
		return *mat;
	}

	const fm::mem_col_dense_matrix &get_data() const {
		return *mat;
	}

	void Random () {
		ep_mat->Random();
		copy_ep2fm(*ep_mat, *mat);
		verify();
//		mat->set_data(rand_init());
	}

	void sync_ep2fm() {
		printf("sync ep to fm %s\n", get_name().c_str());
		for (int i = 0; i < ep_mat->NumVectors(); i++)
			memcpy(mat->get_col(i), (*ep_mat)[i],
					mat->get_num_rows() * mat->get_entry_size());
	}

	void sync_fm2ep() {
		printf("sync fm to ep %s\n", get_name().c_str());
		for (int i = 0; i < ep_mat->NumVectors(); i++)
			memcpy((*ep_mat)[i], mat->get_col(i),
					mat->get_num_rows() * mat->get_entry_size());
	}

	std::string get_name() const {
		return name;
	}

	//! @name Creation methods

	/// \brief Create a new MultiVec with \c numvecs columns.
	/// \return Pointer to the new multivector with uninitialized values.
	virtual Anasazi::MultiVec<ScalarType> * Clone(const int numvecs) const {
		FM_MultiVector<ScalarType> *ret = new FM_MultiVector<ScalarType>(mat->get_num_rows(), numvecs);
		printf("create new %s(#cols: %d)\n", ret->get_name().c_str(), numvecs);
		return ret;
	}

	/// \brief Create a new MultiVec and copy contents of \c *this into it (deep copy).
	/// \return Pointer to the new multivector	
	virtual Anasazi::MultiVec<ScalarType> * CloneCopy () const {
		printf("deep copy %s(#cols: %ld)\n", name.c_str(), mat->get_num_cols());
		std::string extra = std::string("(deep copy from ") + get_name() + ")";
		FM_MultiVector<ScalarType> *ret = new FM_MultiVector<ScalarType>(extra);
		ret->mat = fm::mem_col_dense_matrix::cast(this->mat->deep_copy());
		ret->ep_mat = std::shared_ptr<Anasazi::EpetraMultiVec>(
				dynamic_cast<Anasazi::EpetraMultiVec *>(ep_mat->CloneCopy()));
		verify();
		ret->verify();
		return ret;
	}

	/*! \brief Creates a new Anasazi::MultiVec and copies the selected contents of \c *this 
	  into the new vector (deep copy).  The copied 
	  vectors from \c *this are indicated by the \c index.size() indices in \c index.

	  \return Pointer to the new multivector	
	  */
	virtual Anasazi::MultiVec<ScalarType> * CloneCopy (
			const std::vector<int>& index) const {
		printf("deep copy sub %s (#cols: %ld)\n", name.c_str(), index.size());
		std::string extra = std::string("(deep copy from sub ") + get_name() + ")";
		FM_MultiVector<ScalarType> *ret = new FM_MultiVector<ScalarType>(extra);
		std::vector<off_t> offs(index.begin(), index.end());
		ret->mat = fm::mem_col_dense_matrix::cast(
				fm::mem_col_dense_matrix::cast(this->mat->get_cols(offs))->deep_copy());
		ret->ep_mat = std::shared_ptr<Anasazi::EpetraMultiVec>(
				dynamic_cast<Anasazi::EpetraMultiVec *>(ep_mat->CloneCopy(index)));
		verify();
		ret->verify();
		return ret;
	}

	/*! \brief Creates a new Anasazi::MultiVec that shares the selected contents of \c *this.
	  The index of the \c numvecs vectors shallow copied from \c *this are indicated by the
	  indices given in \c index.

	  \return Pointer to the new multivector	
	  */
	virtual Anasazi::MultiVec<ScalarType> * CloneViewNonConst (
			const std::vector<int>& index) {
		verify();
		std::string extra = std::string("(view of ") + get_name() + ")";
		FM_MultiVector<ScalarType> *ret = new FM_MultiVector<ScalarType>(extra);
		std::vector<off_t> offs(index.begin(), index.end());
		// TODO this deep copy might be unnecessary.
		ret->mat = fm::mem_col_dense_matrix::cast(this->mat->get_cols(offs));
		ret->ep_mat = std::shared_ptr<Anasazi::EpetraMultiVec>(
				dynamic_cast<Anasazi::EpetraMultiVec *>(ep_mat->CloneViewNonConst(index)));
		verify((fm::mem_col_dense_matrix &) *mat->get_cols(offs), *ret->mat);
		verify();
		ret->verify();
		return ret;
	}

	/*! \brief Creates a new Anasazi::MultiVec that shares the selected contents of \c *this.
	  The index of the \c numvecs vectors shallow copied from \c *this are indicated by the
	  indices given in \c index.

	  \return Pointer to the new multivector	
	  */
	virtual const Anasazi::MultiVec<ScalarType> * CloneView (
			const std::vector<int>& index) const {
		verify();
		std::string extra = std::string("(const view of ") + get_name() + ")";
		FM_MultiVector<ScalarType> *ret = new FM_MultiVector<ScalarType>(extra);
		std::vector<off_t> offs(index.begin(), index.end());
		// TODO this deep copy might be unnecessary.
		ret->mat = fm::mem_col_dense_matrix::cast(this->mat->get_cols(offs));
		ret->ep_mat = std::shared_ptr<Anasazi::EpetraMultiVec>(
				// TODO can this be all right?
				dynamic_cast<Anasazi::EpetraMultiVec *>(ep_mat->CloneViewNonConst(index)));
		verify((fm::mem_col_dense_matrix &) *mat->get_cols(offs), *ret->mat);
		verify();
		ret->verify();
		return ret;
	}

	//! @name Dimension information methods	

	//! The number of rows in the multivector.
	ANASAZI_DEPRECATED virtual int GetVecLength () const {
		return mat->get_num_rows();
	}

	//! The number of rows in the multivector.
	//! \note This method supersedes GetVecLength, which will be deprecated.
	virtual ptrdiff_t GetGlobalLength () const {
		return mat->get_num_rows();
	}

	//! The number of vectors (i.e., columns) in the multivector.
	virtual int GetNumberVecs () const {
		return mat->get_num_cols();
	}

	//! @name Update methods

	//! Update \c *this with \c alpha * \c A * \c B + \c beta * (\c *this).
	virtual void MvTimesMatAddMv (ScalarType alpha, 
			const Anasazi::MultiVec<ScalarType>& A, 
			const Teuchos::SerialDenseMatrix<int,ScalarType>& B, ScalarType beta) {
		const FM_MultiVector &fm_A = dynamic_cast<const FM_MultiVector &>(A);
		this->verify();
		fm::mem_row_dense_matrix::ptr aB = fm::mem_row_dense_matrix::create(
				B.numRows(), B.numCols(), fm::get_scalar_type<ScalarType>());
		for (int i = 0; i < B.numRows(); i++) {
			for (int j = 0; j < B.numCols(); j++) {
				if (alpha == 1)
					aB->set<ScalarType>(i, j, B(i, j));
				else
					aB->set<ScalarType>(i, j, B(i, j) * alpha);
			}
		}
		fm::dense_matrix::ptr aAB = fm_A.mat->multiply(*aB);
		if (beta == 0)
			this->mat->copy_from(*aAB);
		else if (beta == 1)
			this->mat->copy_from(*aAB->add(*this->mat));
		else {
			fm::dense_matrix::ptr bThis = this->mat->multiply_scalar(beta);
			this->mat->copy_from(*aAB->add(*bThis));
		}
		this->ep_mat->MvTimesMatAddMv(alpha, *fm_A.ep_mat, B, beta);
		this->sync_ep2fm();
		fm_A.verify();
		this->verify();
	}

	//! Replace \c *this with \c alpha * \c A + \c beta * \c B.
	virtual void MvAddMv ( ScalarType alpha, const Anasazi::MultiVec<ScalarType>& A,
			ScalarType beta, const Anasazi::MultiVec<ScalarType>& B ) {
		const FM_MultiVector &fm_A = dynamic_cast<const FM_MultiVector &>(A);
		const FM_MultiVector &fm_B = dynamic_cast<const FM_MultiVector &>(B);
		fm::dense_matrix::ptr aA = fm_A.mat->multiply_scalar(alpha);
		fm::dense_matrix::ptr bB = fm_B.mat->multiply_scalar(beta);
		this->mat->copy_from(*aA->add(*bB));
		this->ep_mat->MvAddMv(alpha, *fm_A.ep_mat, beta, *fm_B.ep_mat);
		fm_A.verify();
		fm_B.verify();
		this->verify();
	}

	//! Scale each element of the vectors in \c *this with \c alpha.
	virtual void MvScale ( ScalarType alpha ) {
		mat->copy_from(*mat->multiply_scalar<ScalarType>(alpha));
		ep_mat->MvScale(alpha);
		this->verify();
	}

	//! Scale each element of the <tt>i</tt>-th vector in \c *this with <tt>alpha[i]</tt>.
	virtual void MvScale ( const std::vector<ScalarType>& alpha ) {
		mat->scale_cols<ScalarType>(alpha);
		ep_mat->MvScale(alpha);
		this->verify();
	}

	/*! \brief Compute a dense matrix \c B through the matrix-matrix multiply 
	  \c alpha * \c A^T * (\c *this).
	  */
	virtual void MvTransMv ( ScalarType alpha, const Anasazi::MultiVec<ScalarType>& A,
			Teuchos::SerialDenseMatrix<int,ScalarType>& B
#ifdef HAVE_ANASAZI_EXPERIMENTAL
			, ConjType conj = Anasazi::CONJ
#endif				 
			) const {
		const FM_MultiVector &fm_A = dynamic_cast<const FM_MultiVector &>(A);
		assert((size_t) B.numRows() == fm_A.mat->get_num_cols());
		assert((size_t) B.numCols() == this->mat->get_num_cols());
		assert(fm_A.mat->get_num_rows() == this->mat->get_num_rows());
		fm::mem_dense_matrix::ptr tA = fm::mem_dense_matrix::cast(fm_A.mat->transpose());
		fm::mem_dense_matrix::ptr res = fm::mem_dense_matrix::cast(tA->multiply(*this->mat));
		for (int i = 0; i < B.numRows(); i++) {
			for (int j = 0; j < B.numCols(); j++) {
				B(i, j) = res->get<ScalarType>(i, j) * alpha;
			}
		}
	}

	/// \brief Compute the dot product of each column of *this with the corresponding column of A.
	///
	/// Compute a vector \c b whose entries are the individual
	/// dot-products.  That is, <tt>b[i] = A[i]^H * (*this)[i]</tt>
	/// where <tt>A[i]</tt> is the i-th column of A.
	virtual void MvDot ( const Anasazi::MultiVec<ScalarType>& A,
			std::vector<ScalarType> & b 
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
			std::vector<typename Teuchos::ScalarTraits<ScalarType>::magnitudeType> & normvec) const {
		printf("norm(%s)(#cols: %ld)\n", name.c_str(), normvec.size());
		verify();
		std::vector<off_t> col_idx(1);
		for (size_t i = 0; i < mat->get_num_cols(); i++) {
			col_idx[0] = i;
			normvec[i] = mat->get_cols(col_idx)->norm2();
		}
	}

	//! @name Initialization methods

	/// \brief Copy the vectors in \c A to a set of vectors in \c *this.  
	///
	/// The \c numvecs vectors in \c A are copied to a subset of vectors
	/// in \c *this indicated by the indices given in \c index.
	virtual void SetBlock (const Anasazi::MultiVec<ScalarType>& A,
			const std::vector<int>& index) {
		assert((size_t) A.GetNumberVecs() == index.size());
		const FM_MultiVector &fm_A = dynamic_cast<const FM_MultiVector &>(A);
		std::vector<off_t> offs(index.begin(), index.end());
		this->mat->set_cols(*fm_A.mat, offs);
		this->ep_mat->SetBlock(*fm_A.ep_mat, index);
		fm_A.verify();
		this->verify();
	}

	//! Fill all the vectors in \c *this with random numbers.
	virtual void MvRandom () {
		ep_mat->MvRandom();
		copy_ep2fm(*ep_mat, *mat);
//		mat->set_data(rand_init());
		this->verify();
	}

	//! Replace each element of the vectors in \c *this with \c alpha.
	virtual void MvInit ( ScalarType alpha ) {
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
template<class ScalarType>
class MultiVecTraits<ScalarType, FM_MultiVector<ScalarType> > {
public:
	//! @name Creation methods
	//@{ 

	/// \brief Create a new empty \c MultiVec containing \c numvecs columns.
	/// \return Reference-counted pointer to the new \c MultiVec.
	static Teuchos::RCP<FM_MultiVector<ScalarType> > Clone (
			const FM_MultiVector<ScalarType>& mv, const int numvecs) {
		return Teuchos::rcp ((FM_MultiVector<ScalarType> *)
				const_cast<FM_MultiVector<ScalarType>&> (mv).Clone(numvecs)); 
	}

	/*!
	 * \brief Creates a new \c Anasazi::MultiVec and copies contents of \c mv
	 * into the new vector (deep copy).
	  \return Reference-counted pointer to the new \c Anasazi::MultiVec.
	  */
	static Teuchos::RCP<FM_MultiVector<ScalarType> > CloneCopy(
			const FM_MultiVector<ScalarType>& mv) {
		return Teuchos::rcp((FM_MultiVector<ScalarType> *)
				const_cast<FM_MultiVector<ScalarType>&>(mv).CloneCopy());
	}

	/*! \brief Creates a new \c Anasazi::MultiVec and copies the selected
	 * contents of \c mv into the new vector (deep copy).  

	  The copied vectors from \c mv are indicated by the \c index.size() indices in \c index.      
	  \return Reference-counted pointer to the new \c Anasazi::MultiVec.
	  */
	static Teuchos::RCP<FM_MultiVector<ScalarType> > CloneCopy(
			const FM_MultiVector<ScalarType>& mv, const std::vector<int>& index) {
		return Teuchos::rcp((FM_MultiVector<ScalarType> *)
				const_cast<FM_MultiVector<ScalarType>&>(mv).CloneCopy(index));
	}

	/*! \brief Creates a new \c Anasazi::MultiVec that shares the selected
	 * contents of \c mv (shallow copy).

	  The index of the \c numvecs vectors shallow copied from \c mv are indicated
	  by the indices given in \c index.
	  \return Reference-counted pointer to the new \c Anasazi::MultiVec.
	  */    
	static Teuchos::RCP<FM_MultiVector<ScalarType> > CloneViewNonConst(
			FM_MultiVector<ScalarType>& mv, const std::vector<int>& index) {
		return Teuchos::rcp((FM_MultiVector<ScalarType> *)
				mv.CloneViewNonConst(index));
	}

	/*! \brief Creates a new const \c Anasazi::MultiVec that shares
	 * the selected contents of \c mv (shallow copy).

	  The index of the \c numvecs vectors shallow copied from \c mv are
	  indicated by the indices given in \c index.
	  \return Reference-counted pointer to the new const \c Anasazi::MultiVec.
	  */      
	static Teuchos::RCP<const FM_MultiVector<ScalarType> > CloneView(
			const FM_MultiVector<ScalarType>& mv, const std::vector<int>& index) {
		return Teuchos::rcp((FM_MultiVector<ScalarType> *)
				const_cast<FM_MultiVector<ScalarType>&>(mv).CloneView(index));
	}

	//@}

	//! @name Attribute methods
	//@{ 

	//! Obtain the vector length of \c mv.
	ANASAZI_DEPRECATED static int GetVecLength(const FM_MultiVector<ScalarType>& mv) {
		return mv.GetGlobalLength();
	}

	//! Obtain the number of vectors in \c mv
	static int GetNumberVecs(const FM_MultiVector<ScalarType>& mv) {
		return mv.GetNumberVecs();
	}

	//@}

	//! @name Update methods
	//@{ 

	/*! \brief Update \c mv with \f$ \alpha AB + \beta mv \f$.
	*/
	static void MvTimesMatAddMv( ScalarType alpha, const FM_MultiVector<ScalarType>& A, 
			const Teuchos::SerialDenseMatrix<int,ScalarType>& B, 
			ScalarType beta, FM_MultiVector<ScalarType>& mv )
	{ mv.MvTimesMatAddMv(alpha, A, B, beta); }

	/*! \brief Replace \c mv with \f$\alpha A + \beta B\f$.
	*/
	static void MvAddMv( ScalarType alpha, const FM_MultiVector<ScalarType>& A,
			ScalarType beta, const FM_MultiVector<ScalarType>& B,
			FM_MultiVector<ScalarType>& mv)
	{ mv.MvAddMv(alpha, A, beta, B); }

	/*! \brief Compute a dense matrix \c B through the matrix-matrix multiply \f$ \alpha A^Tmv \f$.
	*/
	static void MvTransMv( ScalarType alpha, const FM_MultiVector<ScalarType>& A,
			const FM_MultiVector<ScalarType>& mv,
			Teuchos::SerialDenseMatrix<int, ScalarType>& B
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
	static void MvDot( const FM_MultiVector<ScalarType>& mv,
			const FM_MultiVector<ScalarType>& A, std::vector<ScalarType> & b
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
	static void MvScale ( FM_MultiVector<ScalarType>& mv, ScalarType alpha )
	{ mv.MvScale( alpha ); }

	//! Scale each element of the \c i-th vector in \c *this with \c alpha[i].
	static void MvScale ( FM_MultiVector<ScalarType>& mv,
			const std::vector<ScalarType>& alpha )
	{ mv.MvScale( alpha ); }

	//@}
	//! @name Norm method
	//@{ 

	/*! \brief Compute the 2-norm of each individual vector of \c mv.  
	  Upon return, \c normvec[i] holds the value of \f$||mv_i||_2\f$, the \c i-th column of \c mv.
	  */
	static void MvNorm( const FM_MultiVector<ScalarType>& mv,
			std::vector<typename Teuchos::ScalarTraits<ScalarType>::magnitudeType> & normvec)
	{ mv.MvNorm(normvec); }

	//@}
	//! @name Initialization methods
	//@{ 
	/*! \brief Copy the vectors in \c A to a set of vectors in \c mv indicated by the indices given in \c index.

	  The \c numvecs vectors in \c A are copied to a subset of vectors in \c mv indicated by the indices given in \c index,
	  i.e.<tt> mv[index[i]] = A[i]</tt>.
	  */
	static void SetBlock( const FM_MultiVector<ScalarType>& A,
			const std::vector<int>& index, FM_MultiVector<ScalarType>& mv )
	{ mv.SetBlock(A, index); }

	/*! \brief Replace the vectors in \c mv with random vectors.
	*/
	static void MvRandom( FM_MultiVector<ScalarType>& mv )
	{ mv.MvRandom(); }

	/*! \brief Replace each element of the vectors in \c mv with \c alpha.
	*/
	static void MvInit( FM_MultiVector<ScalarType>& mv,
			ScalarType alpha = Teuchos::ScalarTraits<ScalarType>::zero() )
	{ mv.MvInit(alpha); }

	//@}
	//! @name Print method
	//@{ 

	//! Print the \c mv multi-vector to the \c os output stream.
	static void MvPrint( const FM_MultiVector<ScalarType>& mv, std::ostream& os )
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
	typedef details::MultiVecTsqrAdapter<ScalarType> tsqr_adaptor_type;
#endif // HAVE_ANASAZI_TSQR
};

}

#endif
