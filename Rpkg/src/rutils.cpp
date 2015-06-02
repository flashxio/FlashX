/*
 * Copyright 2014 Open Connectome Project (http://openconnecto.me)
 * Written by Da Zheng (zhengda1936@gmail.com)
 *
 * This file is part of FlashR.
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

#include <Rinternals.h>

#include "rutils.h"

/*
 * We need to implement these two functions here because Rinternals.h masses
 * up with other C++ header files.
 */

bool R_is_real(SEXP v)
{
	return isReal(v);
}

bool R_is_integer(SEXP v)
{
	return isInteger(v);
}

void R_gc()
{
	SEXP call = PROTECT(lang1(install("gc")));
	PROTECT(eval(call, R_GlobalEnv));
	UNPROTECT(2);
}
