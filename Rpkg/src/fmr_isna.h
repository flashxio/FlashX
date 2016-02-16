// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; tab-width: 8 -*-
/* :tabSize=4:indentSize=4:noTabs=false:folding=explicit:collapseFolds=1: */
//
// na.h: Rcpp R/C++ interface class library -- optimized na checking
//
// Copyright (C) 2012-2014 Dirk Eddelbuettel, Romain Francois, Kevin Ushey
// and Da Zheng
//
// This file is part of Rcpp.
//
// Rcpp is free software: you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 2 of the License, or
// (at your option) any later version.
//
// Rcpp is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Rcpp.  If not, see <http://www.gnu.org/licenses/>.

#include <Rcpp.h>
#include <boost/static_assert.hpp>

#ifdef RCPP_HAS_LONG_LONG_TYPES

BOOST_STATIC_ASSERT_MSG(sizeof(rcpp_ulong_long_type) == sizeof(double),
		"unsigned long long and double have same size");

// motivation: on 32bit architectures, we only see 'LargeNA'
// as defined ahead; on 64bit architectures, R defaults to
// 'SmallNA' for R_NaReal, but this can get promoted to 'LargeNA'
// if a certain operation can create a 'signalling' NA, e.g. NA_real_+1
static const rcpp_ulong_long_type SmallNA = 0x7FF00000000007A2;
static const rcpp_ulong_long_type LargeNA = 0x7FF80000000007A2;

struct NACanChange {
	enum { value = sizeof(void*) == 8 };
};

template <bool NACanChange>
bool FM_IsNA__impl(const rcpp_ulong_long_type *x);

template <>
inline bool FM_IsNA__impl<true>(const rcpp_ulong_long_type *x) {
	return *x == SmallNA || *x == LargeNA;
}

template <>
inline bool FM_IsNA__impl<false>(const rcpp_ulong_long_type *x) {
	return *x == LargeNA;
}

inline bool FM_IsNA(const rcpp_ulong_long_type *x) {
	return FM_IsNA__impl< NACanChange::value >(x);
}

#endif
