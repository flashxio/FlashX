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

#include <boost/format.hpp>

#include "log.h"

#include "AnasaziMultiVec.hpp"

#include "NUMA_vector.h"
#include "matrix_config.h"
#include "generic_type.h"

#ifdef FM_VERIFY
#include "Epetra_MultiVector.h"
#include "Epetra_SerialComm.h"
#include "Epetra_Map.h"
#include "AnasaziEpetraAdapter.hpp"
#endif
#include "block_dense_matrix.h"

#include "mem_matrix_store.h"
#include "dense_matrix.h"
#include "dotp_matrix_store.h"
#include "matrix_stats.h"

static int MV_id;

namespace fm
{

namespace eigen
{

template<class ScalarType>
class FM_MultiVector: public Anasazi::MultiVec<ScalarType>
{
	std::string solver;
	bool in_mem;
	std::string name;
	block_multi_vector::ptr mat;
#ifdef FM_VERIFY
	std::shared_ptr<Anasazi::EpetraMultiVec> ep_mat;
#endif

	FM_MultiVector(const std::string &extra, bool in_mem,
			const std::string &solver) {
		char name_buf[128];
		snprintf(name_buf, sizeof(name_buf), "MV-%d", MV_id++);
		this->name = std::string(name_buf) + " " + extra;
		this->in_mem = in_mem;
		this->solver = solver;
	}
public:
	FM_MultiVector(size_t num_rows, size_t num_cols, size_t block_size,
			bool in_mem, const std::string &solver) {
		this->solver = solver;
		this->in_mem = in_mem;
		// We don't materialize the column matrix.
		mat = block_multi_vector::create(num_rows, num_cols, block_size,
				fm::get_scalar_type<ScalarType>(), in_mem);

#ifdef FM_VERIFY
		Epetra_SerialComm Comm;
		Teuchos::RCP<Epetra_Map> Map = Teuchos::rcp (
				new Epetra_Map((int) num_rows, 0, Comm));
		ep_mat = std::shared_ptr<Anasazi::EpetraMultiVec>(
				new Anasazi::EpetraMultiVec(*Map, num_cols));
#endif

		char name_buf[128];
		snprintf(name_buf, sizeof(name_buf), "MV-%d", MV_id++);
		this->name = name_buf;
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

#ifdef FM_VERIFY
	Anasazi::EpetraMultiVec &get_ep_mv() {
		return *ep_mat;
	}

	const Anasazi::EpetraMultiVec &get_ep_mv() const {
		return *ep_mat;
	}
#endif

	block_multi_vector::ptr get_data() {
		return mat;
	}

	block_multi_vector::ptr get_data() const {
		return mat;
	}

	size_t get_block_size() const {
		return mat->get_block_size();
	}

	void Random () {
		mat->init_rand<ScalarType>(-1, 1);
		sync_fm2ep();
		verify();
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

	std::string get_name() const {
		return name;
	}

	//! @name Creation methods

	/// \brief Create a new MultiVec with \c numvecs columns.
	/// \return Pointer to the new multivector with uninitialized values.
	virtual Anasazi::MultiVec<ScalarType> * Clone(const int numvecs) const {
		FM_MultiVector<ScalarType> *ret;
		if (numvecs % mat->get_block_size() == 0)
			ret = new FM_MultiVector<ScalarType>(mat->get_num_rows(), numvecs,
					mat->get_block_size(), in_mem, solver);
		else
			ret = new FM_MultiVector<ScalarType>(
					mat->get_num_rows(), numvecs, numvecs, in_mem, solver);
		BOOST_LOG_TRIVIAL(info) << boost::format("create new %1% (#cols: %2%)")
			% ret->get_name() % numvecs;
		return ret;
	}

	/// \brief Create a new MultiVec and copy contents of \c *this into it (deep copy).
	/// \return Pointer to the new multivector	
	virtual Anasazi::MultiVec<ScalarType> * CloneCopy () const {
		num_col_writes_concept += mat->get_num_cols();
		num_col_reads_concept += mat->get_num_cols();

		BOOST_LOG_TRIVIAL(info) << boost::format("deep copy %1% (#cols: %2%)")
			% name % mat->get_num_cols();
		std::string extra = std::string("(deep copy from ") + get_name() + ")";
		FM_MultiVector<ScalarType> *ret = new FM_MultiVector<ScalarType>(extra,
				in_mem, solver);
		ret->mat = this->mat->clone();
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
	virtual Anasazi::MultiVec<ScalarType> * CloneCopy (
			const std::vector<int>& index) const {
		num_col_writes_concept += index.size();
		num_col_reads_concept += index.size();

		BOOST_LOG_TRIVIAL(info) << boost::format("deep copy sub %1% (#cols: %2%)")
			% name % index.size();
		std::string extra = std::string("(deep copy from sub ") + get_name() + ")";
		FM_MultiVector<ScalarType> *ret = new FM_MultiVector<ScalarType>(extra,
				in_mem, solver);
		ret->mat = mat->get_cols(index);
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
	virtual Anasazi::MultiVec<ScalarType> * CloneViewNonConst (
			const std::vector<int>& index) {
		std::string extra = std::string("(") + get_name() + "[" + vec2str(index) + "])";
		FM_MultiVector<ScalarType> *ret = new FM_MultiVector<ScalarType>(extra,
				in_mem, solver);
		BOOST_LOG_TRIVIAL(info) << boost::format("view %1% (#cols: %2%)")
			% ret->name % index.size();
		ret->mat = mat->get_cols_mirror(index);
#ifdef FM_VERIFY
		ret->ep_mat = std::shared_ptr<Anasazi::EpetraMultiVec>(
				dynamic_cast<Anasazi::EpetraMultiVec *>(ep_mat->CloneViewNonConst(index)));
#endif
//		verify();
//		ret->verify();
		return ret;
	}

	static std::string vec2str(const std::vector<int> &index) {
		std::vector<int> copy = index;
		std::sort(copy.begin(), copy.end());
		typedef std::pair<int, int> range_t;
		std::vector<range_t> ranges;
		ranges.push_back(range_t(copy.front(), copy.front()));
		for (size_t i = 1; i < copy.size(); i++) {
			if (ranges.back().second + 1 == copy[i])
				ranges.back().second++;
			else
				ranges.push_back(range_t(copy[i], copy[i]));
		}
		std::string ret = "";
		for (size_t i = 0; i < ranges.size(); i++) {
			if (ranges[i].first == ranges[i].second)
				ret += itoa(ranges[i].first);
			else
				ret += std::string(std::string(itoa(ranges[i].first)) + ":") + itoa(ranges[i].second);
		}
		return ret;
	}

	/*! \brief Creates a new Anasazi::MultiVec that shares the selected contents of \c *this.
	  The index of the \c numvecs vectors shallow copied from \c *this are indicated by the
	  indices given in \c index.

	  \return Pointer to the new multivector	
	  */
	virtual const Anasazi::MultiVec<ScalarType> * CloneView (
			const std::vector<int>& index) const {
		std::string extra = std::string("(const ") + get_name() + "[" + vec2str(index) + "])";
		FM_MultiVector<ScalarType> *ret = new FM_MultiVector<ScalarType>(extra,
				in_mem, solver);
		if (index.size() == 1
				&& (solver == "KrylovSchur" || solver == "Davidson")) {
			size_t block_idx = index[0] / mat->get_block_size();
			dense_matrix::ptr block = mat->get_block(block_idx);
			const dotp_matrix_store *store
				= dynamic_cast<const dotp_matrix_store *>(block->get_raw_store().get());
			if (store == NULL) {
				dotp_matrix_store::ptr dotp = dotp_matrix_store::create(
						mat->get_block(block_idx)->get_raw_store());
				mat->set_block(block_idx, fm::dense_matrix::create(dotp));
			}
		}
		BOOST_LOG_TRIVIAL(info) << boost::format("const view %1% (#cols: %2%)")
			% ret->name % index.size();
		ret->mat = mat->get_cols(index);
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
		sync_fm2ep();
		const FM_MultiVector &fm_A = dynamic_cast<const FM_MultiVector &>(A);
		BOOST_LOG_TRIVIAL(info) << boost::format(
				"this(%1%) = %2% * A(%3%) * B(%4%x%5%) + %6% * this")
			% name % alpha % fm_A.name % B.numRows() % B.numCols() % beta;
		if (alpha != 0)
			num_col_reads_concept += fm_A.mat->get_num_cols();
		if (beta != 0)
			num_col_reads_concept += mat->get_num_cols();
		num_col_writes_concept += mat->get_num_cols();

		fm::detail::mem_col_matrix_store::ptr Bstore
			= fm::detail::mem_col_matrix_store::create(
					B.numRows(), B.numCols(), fm::get_scalar_type<ScalarType>());
		for (int i = 0; i < B.numRows(); i++) {
			for (int j = 0; j < B.numCols(); j++)
				Bstore->set<ScalarType>(i, j, B(i, j));
		}
		fm::scalar_variable_impl<ScalarType> alpha_var(alpha);
		fm::scalar_variable_impl<ScalarType> beta_var(beta);
		this->mat->assign(*this->mat->gemm(*fm_A.mat, Bstore, alpha_var, beta_var));
#ifdef FM_VERIFY
		this->ep_mat->MvTimesMatAddMv(alpha, *fm_A.ep_mat, B, beta);
#endif
		fm_A.verify();
		this->verify();
	}

	//! Replace \c *this with \c alpha * \c A + \c beta * \c B.
	virtual void MvAddMv ( ScalarType alpha, const Anasazi::MultiVec<ScalarType>& A,
			ScalarType beta, const Anasazi::MultiVec<ScalarType>& B ) {
		sync_fm2ep();
		const FM_MultiVector &fm_A = dynamic_cast<const FM_MultiVector &>(A);
		const FM_MultiVector &fm_B = dynamic_cast<const FM_MultiVector &>(B);
		BOOST_LOG_TRIVIAL(info) << boost::format(
				"this(%1%) = %2% * A(%3%) + %4% *  B(%5%)")
			% name % alpha % fm_A.name % beta % fm_B.name;
		num_col_reads_concept += fm_A.mat->get_num_cols() + fm_B.mat->get_num_cols();
		num_col_writes_concept += mat->get_num_cols();

		if (alpha == 1 && beta == 0)
			this->mat->assign(*fm_A.mat);
		else if (alpha == 0 && beta == 1)
			this->mat->assign(*fm_B.mat);
		else {
			block_multi_vector::ptr aA = fm_A.mat->multiply_scalar(alpha);
			block_multi_vector::ptr bB = fm_B.mat->multiply_scalar(beta);
			this->mat->assign(*aA->add(*bB));
		}
#ifdef FM_VERIFY
		this->ep_mat->MvAddMv(alpha, *fm_A.ep_mat, beta, *fm_B.ep_mat);
#endif
		fm_A.verify();
		fm_B.verify();
		this->verify();
	}

	//! Scale each element of the vectors in \c *this with \c alpha.
	virtual void MvScale ( ScalarType alpha ) {
		num_col_writes_concept += mat->get_num_cols();
		num_col_reads_concept += mat->get_num_cols();

		sync_fm2ep();
		BOOST_LOG_TRIVIAL(info) << boost::format("this(%1%) *= %2%") % name % alpha;
		mat->assign(*mat->multiply_scalar<ScalarType>(alpha));
#ifdef FM_VERIFY
		ep_mat->MvScale(alpha);
#endif
		this->verify();
	}

	//! Scale each element of the <tt>i</tt>-th vector in \c *this with <tt>alpha[i]</tt>.
	virtual void MvScale ( const std::vector<ScalarType>& alpha ) {
		num_col_writes_concept += mat->get_num_cols();
		num_col_reads_concept += mat->get_num_cols();

		sync_fm2ep();
		BOOST_LOG_TRIVIAL(info) << boost::format("this(%s) *= vec") % name;
		mat->assign(*mat->scale_cols<ScalarType>(alpha));
#ifdef FM_VERIFY
		ep_mat->MvScale(alpha);
#endif
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
		BOOST_LOG_TRIVIAL(info) << boost::format(
				"B(%1%x%2%) = %3% * A(%4%)^T * this(%5%)")
			% B.numRows() % B.numCols() % alpha % fm_A.name % name;
		num_col_reads_concept += fm_A.mat->get_num_cols() + mat->get_num_cols();

		assert((size_t) B.numRows() == fm_A.mat->get_num_cols());
		assert((size_t) B.numCols() == this->mat->get_num_cols());
		assert(fm_A.mat->get_num_rows() == this->mat->get_num_rows());
		this->verify();

		long double lalpha = alpha;
		// This is a special case: we compute the dot product of a col vector
		// with itself.
		if (fm_A.mat->get_num_cols() == 1 && mat->get_num_cols() == 1) {
			const sub_dotp_matrix_store *col1
				= dynamic_cast<const sub_dotp_matrix_store *>(
						fm_A.mat->get_block(0)->get_raw_store().get());
			const sub_dotp_matrix_store *col2
				= dynamic_cast<const sub_dotp_matrix_store *>(
						mat->get_block(0)->get_raw_store().get());
			if (col1 && col2
					&& col1->get_orig_store() == col2->get_orig_store()) {
				B(0, 0) = col1->get_col_dot(0) * lalpha;
				return;
			}
		}

		fm::dense_matrix::ptr res = mat->MvTransMv(*fm_A.mat);
		const_cast<FM_MultiVector *>(this)->sync_fm2ep();
#if 0
		ep_mat->MvTransMv(alpha, *fm_A.ep_mat, B);
		for (int i = 0; i < B.numRows(); i++)
			for (int j = 0; j < B.numCols(); j++) {
				printf("%g\n", B(i, j) - (double) (res->get<ScalarType>(i, j) * lalpha));
				assert(B(i, j) == (double) (res->get<ScalarType>(i, j) * lalpha));
			}
#endif
		const fm::detail::mem_matrix_store &mem_res
			= dynamic_cast<const fm::detail::mem_matrix_store &>(res->get_data());
		for (int i = 0; i < B.numRows(); i++) {
			for (int j = 0; j < B.numCols(); j++) {
				B(i, j) = mem_res.get<ScalarType>(i, j) * lalpha;
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
		num_col_reads_concept += mat->get_num_cols();

		BOOST_LOG_TRIVIAL(info) << boost::format("norm(%1%)(#cols: %2%)")
			% name % normvec.size();
		verify();
		fm::detail::matrix_stats_t orig_stats = fm::detail::matrix_stats;
		for (size_t i = 0; i < mat->get_num_blocks(); i++) {
			printf("materialize %s on the fly\n",
					mat->get_block(i)->get_data().get_name().c_str());
			dotp_matrix_store::ptr dotp
				= dotp_matrix_store::create(mat->get_block(i)->get_raw_store());
			mat->set_block(i, fm::dense_matrix::create(dotp));
			std::vector<ScalarType> col_dots = dotp->get_col_dot_prods();
			// TODO do I need long double here?
			for (size_t j = 0; j < col_dots.size(); j++)
				normvec[i * mat->get_block_size() + j] = std::sqrt(col_dots[j]);
		}
		fm::detail::matrix_stats.print_diff(orig_stats);
		const_cast<FM_MultiVector *>(this)->sync_fm2ep();
#ifdef FM_VERIFY
		std::vector<typename Teuchos::ScalarTraits<ScalarType>::magnitudeType> normvec1(
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
	virtual void SetBlock (const Anasazi::MultiVec<ScalarType>& A,
			const std::vector<int>& index) {
		num_col_reads_concept += index.size();
		num_col_writes_concept += index.size();

		sync_fm2ep();
		assert((size_t) A.GetNumberVecs() == index.size());
		const FM_MultiVector &fm_A = dynamic_cast<const FM_MultiVector &>(A);
		BOOST_LOG_TRIVIAL(info) << boost::format("this(%1%)[%2%] = A(%3%)")
			% name % vec2str(index) % fm_A.name;
		this->mat->set_block(*fm_A.mat, index);
#ifdef FM_VERIFY
		this->ep_mat->SetBlock(*fm_A.ep_mat, index);
#endif
		fm_A.verify();
		this->verify();
	}

	//! Fill all the vectors in \c *this with random numbers.
	virtual void MvRandom () {
		num_col_writes_concept += mat->get_num_cols();

		BOOST_LOG_TRIVIAL(info) << boost::format("this(%1%) = random") % name;
		mat->init_rand<ScalarType>(-1, 1);
//		ep_mat->MvRandom();
		sync_fm2ep();
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

}

}

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
class MultiVecTraits<ScalarType, fm::eigen::FM_MultiVector<ScalarType> > {
public:
	//! @name Creation methods
	//@{ 

	/// \brief Create a new empty \c MultiVec containing \c numvecs columns.
	/// \return Reference-counted pointer to the new \c MultiVec.
	static Teuchos::RCP<fm::eigen::FM_MultiVector<ScalarType> > Clone (
			const fm::eigen::FM_MultiVector<ScalarType>& mv, const int numvecs) {
		return Teuchos::rcp ((fm::eigen::FM_MultiVector<ScalarType> *)
				const_cast<fm::eigen::FM_MultiVector<ScalarType>&> (mv).Clone(numvecs)); 
	}

	/*!
	 * \brief Creates a new \c Anasazi::MultiVec and copies contents of \c mv
	 * into the new vector (deep copy).
	  \return Reference-counted pointer to the new \c Anasazi::MultiVec.
	  */
	static Teuchos::RCP<fm::eigen::FM_MultiVector<ScalarType> > CloneCopy(
			const fm::eigen::FM_MultiVector<ScalarType>& mv) {
		return Teuchos::rcp((fm::eigen::FM_MultiVector<ScalarType> *)
				const_cast<fm::eigen::FM_MultiVector<ScalarType>&>(mv).CloneCopy());
	}

	/*! \brief Creates a new \c Anasazi::MultiVec and copies the selected
	 * contents of \c mv into the new vector (deep copy).  

	  The copied vectors from \c mv are indicated by the \c index.size() indices in \c index.      
	  \return Reference-counted pointer to the new \c Anasazi::MultiVec.
	  */
	static Teuchos::RCP<fm::eigen::FM_MultiVector<ScalarType> > CloneCopy(
			const fm::eigen::FM_MultiVector<ScalarType>& mv,
			const std::vector<int>& index) {
		return Teuchos::rcp((fm::eigen::FM_MultiVector<ScalarType> *)
				const_cast<fm::eigen::FM_MultiVector<ScalarType>&>(
					mv).CloneCopy(index));
	}

	/*! \brief Creates a new \c Anasazi::MultiVec that shares the selected
	 * contents of \c mv (shallow copy).

	  The index of the \c numvecs vectors shallow copied from \c mv are indicated
	  by the indices given in \c index.
	  \return Reference-counted pointer to the new \c Anasazi::MultiVec.
	  */    
	static Teuchos::RCP<fm::eigen::FM_MultiVector<ScalarType> > CloneViewNonConst(
			fm::eigen::FM_MultiVector<ScalarType>& mv,
			const std::vector<int>& index) {
		return Teuchos::rcp((fm::eigen::FM_MultiVector<ScalarType> *)
				mv.CloneViewNonConst(index));
	}

	/*! \brief Creates a new const \c Anasazi::MultiVec that shares
	 * the selected contents of \c mv (shallow copy).

	  The index of the \c numvecs vectors shallow copied from \c mv are
	  indicated by the indices given in \c index.
	  \return Reference-counted pointer to the new const \c Anasazi::MultiVec.
	  */      
	static Teuchos::RCP<const fm::eigen::FM_MultiVector<ScalarType> > CloneView(
			const fm::eigen::FM_MultiVector<ScalarType>& mv,
			const std::vector<int>& index) {
		return Teuchos::rcp((fm::eigen::FM_MultiVector<ScalarType> *)
				const_cast<fm::eigen::FM_MultiVector<ScalarType>&>(
					mv).CloneView(index));
	}

	//@}

	//! @name Attribute methods
	//@{ 

	//! Obtain the vector length of \c mv.
	ANASAZI_DEPRECATED static int GetVecLength(
			const fm::eigen::FM_MultiVector<ScalarType>& mv) {
		return mv.GetGlobalLength();
	}

	//! Obtain the number of vectors in \c mv
	static int GetNumberVecs(const fm::eigen::FM_MultiVector<ScalarType>& mv) {
		return mv.GetNumberVecs();
	}

	//@}

	//! @name Update methods
	//@{ 

	/*! \brief Update \c mv with \f$ \alpha AB + \beta mv \f$.
	*/
	static void MvTimesMatAddMv( ScalarType alpha,
			const fm::eigen::FM_MultiVector<ScalarType>& A, 
			const Teuchos::SerialDenseMatrix<int,ScalarType>& B, 
			ScalarType beta, fm::eigen::FM_MultiVector<ScalarType>& mv )
	{ mv.MvTimesMatAddMv(alpha, A, B, beta); }

	/*! \brief Replace \c mv with \f$\alpha A + \beta B\f$.
	*/
	static void MvAddMv( ScalarType alpha,
			const fm::eigen::FM_MultiVector<ScalarType>& A,
			ScalarType beta, const fm::eigen::FM_MultiVector<ScalarType>& B,
			fm::eigen::FM_MultiVector<ScalarType>& mv)
	{ mv.MvAddMv(alpha, A, beta, B); }

	/*! \brief Compute a dense matrix \c B through the matrix-matrix multiply \f$ \alpha A^Tmv \f$.
	*/
	static void MvTransMv( ScalarType alpha,
			const fm::eigen::FM_MultiVector<ScalarType>& A,
			const fm::eigen::FM_MultiVector<ScalarType>& mv,
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
	static void MvDot( const fm::eigen::FM_MultiVector<ScalarType>& mv,
			const fm::eigen::FM_MultiVector<ScalarType>& A,
			std::vector<ScalarType> & b
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
	static void MvScale ( fm::eigen::FM_MultiVector<ScalarType>& mv, ScalarType alpha )
	{ mv.MvScale( alpha ); }

	//! Scale each element of the \c i-th vector in \c *this with \c alpha[i].
	static void MvScale ( fm::eigen::FM_MultiVector<ScalarType>& mv,
			const std::vector<ScalarType>& alpha )
	{ mv.MvScale( alpha ); }

	//@}
	//! @name Norm method
	//@{ 

	/*! \brief Compute the 2-norm of each individual vector of \c mv.  
	  Upon return, \c normvec[i] holds the value of \f$||mv_i||_2\f$, the \c i-th column of \c mv.
	  */
	static void MvNorm( const fm::eigen::FM_MultiVector<ScalarType>& mv,
			std::vector<typename Teuchos::ScalarTraits<ScalarType>::magnitudeType> & normvec)
	{ mv.MvNorm(normvec); }

	//@}
	//! @name Initialization methods
	//@{ 
	/*! \brief Copy the vectors in \c A to a set of vectors in \c mv indicated by the indices given in \c index.

	  The \c numvecs vectors in \c A are copied to a subset of vectors in \c mv indicated by the indices given in \c index,
	  i.e.<tt> mv[index[i]] = A[i]</tt>.
	  */
	static void SetBlock( const fm::eigen::FM_MultiVector<ScalarType>& A,
			const std::vector<int>& index,
			fm::eigen::FM_MultiVector<ScalarType>& mv )
	{ mv.SetBlock(A, index); }

	/*! \brief Replace the vectors in \c mv with random vectors.
	*/
	static void MvRandom( fm::eigen::FM_MultiVector<ScalarType>& mv )
	{ mv.MvRandom(); }

	/*! \brief Replace each element of the vectors in \c mv with \c alpha.
	*/
	static void MvInit( fm::eigen::FM_MultiVector<ScalarType>& mv,
			ScalarType alpha = Teuchos::ScalarTraits<ScalarType>::zero() )
	{ mv.MvInit(alpha); }

	//@}
	//! @name Print method
	//@{ 

	//! Print the \c mv multi-vector to the \c os output stream.
	static void MvPrint( const fm::eigen::FM_MultiVector<ScalarType>& mv,
			std::ostream& os )
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
