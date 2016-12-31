# generated automatically by aclocal 1.14.1 -*- Autoconf -*-

# Copyright (C) 1996-2013 Free Software Foundation, Inc.

# This file is free software; the Free Software Foundation
# gives unlimited permission to copy and/or distribute it,
# with or without modifications, as long as this notice is preserved.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY, to the extent permitted by law; without
# even the implied warranty of MERCHANTABILITY or FITNESS FOR A
# PARTICULAR PURPOSE.

m4_ifndef([AC_CONFIG_MACRO_DIRS], [m4_defun([_AM_CONFIG_MACRO_DIRS], [])m4_defun([AC_CONFIG_MACRO_DIRS], [_AM_CONFIG_MACRO_DIRS($@)])])
# ===========================================================================
#     http://www.gnu.org/software/autoconf-archive/ax_check_library.html
# ===========================================================================
#
# SYNOPSIS
#
#   AX_CHECK_LIBRARY(VARIABLE-PREFIX, HEADER-FILE, LIBRARY-FILE,
#                    [ACTION-IF-FOUND], [ACTION-IF-NOT_FOUND])
#
# DESCRIPTION
#
#   Provides a generic test for a given library, similar in concept to the
#   PKG_CHECK_MODULES macro used by pkg-config.
#
#   Most simplest libraries can be checked against simply through the
#   presence of a header file and a library to link to. This macro allows to
#   wrap around the test so that it doesn't have to be recreated each time.
#
#   Rather than define --with-$LIBRARY arguments, it uses variables in the
#   same way that PKG_CHECK_MODULES does. It doesn't, though, use the same
#   names, since you shouldn't provide a value for LIBS or CFLAGS but rather
#   for LDFLAGS and CPPFLAGS, to tell the linker and compiler where to find
#   libraries and headers respectively.
#
#   If the library is find, HAVE_PREFIX is defined, and in all cases
#   PREFIX_LDFLAGS and PREFIX_CPPFLAGS are substituted.
#
#   Example:
#
#     AX_CHECK_LIBRARY([LIBEVENT], [event.h], [event], [],
#                      [AC_MSG_ERROR([Unable to find libevent])])
#
# LICENSE
#
#   Copyright (c) 2010 Diego Elio Petteno` <flameeyes@gmail.com>
#
#   This program is free software: you can redistribute it and/or modify it
#   under the terms of the GNU General Public License as published by the
#   Free Software Foundation, either version 3 of the License, or (at your
#   option) any later version.
#
#   This program is distributed in the hope that it will be useful, but
#   WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
#   Public License for more details.
#
#   You should have received a copy of the GNU General Public License along
#   with this program. If not, see <http://www.gnu.org/licenses/>.
#
#   As a special exception, the respective Autoconf Macro's copyright owner
#   gives unlimited permission to copy, distribute and modify the configure
#   scripts that are the output of Autoconf when processing the Macro. You
#   need not follow the terms of the GNU General Public License when using
#   or distributing such scripts, even though portions of the text of the
#   Macro appear in them. The GNU General Public License (GPL) does govern
#   all other use of the material that constitutes the Autoconf Macro.
#
#   This special exception to the GPL applies to versions of the Autoconf
#   Macro released by the Autoconf Archive. When you make and distribute a
#   modified version of the Autoconf Macro, you may extend this special
#   exception to the GPL to apply to your modified version as well.

#serial 4

AC_DEFUN([AX_CHECK_LIBRARY], [
  AC_ARG_VAR($1[_CPPFLAGS], [C preprocessor flags for ]$1[ headers])
  AC_ARG_VAR($1[_LDFLAGS], [linker flags for ]$1[ libraries])

  AC_CACHE_VAL(AS_TR_SH([ax_cv_have_]$1),
    [save_CPPFLAGS="$CPPFLAGS"
     save_LDFLAGS="$LDFLAGS"
     save_LIBS="$LIBS"

     AS_IF([test "x$]$1[_CPPFLAGS" != "x"],
       [CPPFLAGS="$CPPFLAGS $]$1[_CPPFLAGS"])

     AS_IF([test "x$]$1[_LDFLAGS" != "x"],
       [LDFLAGS="$LDFLAGS $]$1[_LDFLAGS"])

     AC_CHECK_HEADER($2, [
       AC_CHECK_LIB($3, [main],
         [AS_TR_SH([ax_cv_have_]$1)=yes],
         [AS_TR_SH([ax_cv_have_]$1)=no])
     ], [AS_TR_SH([ax_cv_have_]$1)=no])

     CPPFLAGS="$save_CPPFLAGS"
     LDFLAGS="$save_LDFLAGS"
     LIBS="$save_LIBS"
    ])

  AS_IF([test "$]AS_TR_SH([ax_cv_have_]$1)[" = "yes"],
    AC_DEFINE([HAVE_]$1, [1], [Define to 1 if ]$1[ is found])
    [$4],
    [$5])
])

# ============================================================================
#  http://www.gnu.org/software/autoconf-archive/ax_cxx_compile_stdcxx_0x.html
# ============================================================================
#
# SYNOPSIS
#
#   AX_CXX_COMPILE_STDCXX_0X
#
# DESCRIPTION
#
#   Check for baseline language coverage in the compiler for the C++0x
#   standard.
#
# LICENSE
#
#   Copyright (c) 2008 Benjamin Kosnik <bkoz@redhat.com>
#
#   Copying and distribution of this file, with or without modification, are
#   permitted in any medium without royalty provided the copyright notice
#   and this notice are preserved. This file is offered as-is, without any
#   warranty.

#serial 7

AU_ALIAS([AC_CXX_COMPILE_STDCXX_0X], [AX_CXX_COMPILE_STDCXX_0X])
AC_DEFUN([AX_CXX_COMPILE_STDCXX_0X], [
  AC_CACHE_CHECK(if g++ supports C++0x features without additional flags,
  ax_cv_cxx_compile_cxx0x_native,
  [AC_LANG_SAVE
  AC_LANG_CPLUSPLUS
  AC_TRY_COMPILE([
  template <typename T>
    struct check
    {
      static_assert(sizeof(int) <= sizeof(T), "not big enough");
    };

    typedef check<check<bool>> right_angle_brackets;

    int a;
    decltype(a) b;

    typedef check<int> check_type;
    check_type c;
    check_type&& cr = static_cast<check_type&&>(c);],,
  ax_cv_cxx_compile_cxx0x_native=yes, ax_cv_cxx_compile_cxx0x_native=no)
  AC_LANG_RESTORE
  ])

  AC_CACHE_CHECK(if g++ supports C++0x features with -std=c++0x,
  ax_cv_cxx_compile_cxx0x_cxx,
  [AC_LANG_SAVE
  AC_LANG_CPLUSPLUS
  ac_save_CXXFLAGS="$CXXFLAGS"
  CXXFLAGS="$CXXFLAGS -std=c++0x"
  AC_TRY_COMPILE([
  template <typename T>
    struct check
    {
      static_assert(sizeof(int) <= sizeof(T), "not big enough");
    };

    typedef check<check<bool>> right_angle_brackets;

    int a;
    decltype(a) b;

    typedef check<int> check_type;
    check_type c;
    check_type&& cr = static_cast<check_type&&>(c);],,
  ax_cv_cxx_compile_cxx0x_cxx=yes, ax_cv_cxx_compile_cxx0x_cxx=no)
  CXXFLAGS="$ac_save_CXXFLAGS"
  AC_LANG_RESTORE
  ])

  AC_CACHE_CHECK(if g++ supports C++0x features with -std=gnu++0x,
  ax_cv_cxx_compile_cxx0x_gxx,
  [AC_LANG_SAVE
  AC_LANG_CPLUSPLUS
  ac_save_CXXFLAGS="$CXXFLAGS"
  CXXFLAGS="$CXXFLAGS -std=gnu++0x"
  AC_TRY_COMPILE([
  template <typename T>
    struct check
    {
      static_assert(sizeof(int) <= sizeof(T), "not big enough");
    };

    typedef check<check<bool>> right_angle_brackets;

    int a;
    decltype(a) b;

    typedef check<int> check_type;
    check_type c;
    check_type&& cr = static_cast<check_type&&>(c);],,
  ax_cv_cxx_compile_cxx0x_gxx=yes, ax_cv_cxx_compile_cxx0x_gxx=no)
  CXXFLAGS="$ac_save_CXXFLAGS"
  AC_LANG_RESTORE
  ])

  if test "$ax_cv_cxx_compile_cxx0x_native" = yes ||
     test "$ax_cv_cxx_compile_cxx0x_cxx" = yes ||
     test "$ax_cv_cxx_compile_cxx0x_gxx" = yes; then
    AC_DEFINE(HAVE_STDCXX_0X,,[Define if g++ supports C++0x features. ])
  fi
])

