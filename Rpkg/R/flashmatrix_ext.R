# Copyright 2014 Open Connectome Project (http://openconnecto.me)
# Written by Da Zheng (zhengda1936@gmail.com)
#
# This file is part of FlashR.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

# This file contains the FlashR functions implemented with
# the FlashR functions in flashmatrix.R.

setMethod("+", signature(e1 = "fm", e2 = "fm"), function(e1, e2)
		  fm.mapply2.fm(e1, e2, fm.bo.add))
setMethod("+", signature(e1 = "fmV", e2 = "fmV"), function(e1, e2)
		  fm.mapply2.fmV(e1, e2, fm.bo.add))
setMethod("+", signature(e1 = "fm", e2 = "fmV"), function(e1, e2)
		  fm.mapply.col(e1, e2, fm.bo.add))
setMethod("+", signature(e1 = "fmV", e2 = "fm"), function(e1, e2)
		  fm.mapply.col(e2, e1, fm.bo.add))
setMethod("+", signature(e1 = "fm", e2 = "matrix"), function(e1, e2)
		  fm.mapply2.fm.m(e1, e2, fm.bo.add))
setMethod("+", signature(e1 = "matrix", e2 = "fm"), function(e1, e2)
		  fm.mapply2.m.fm(e1, e2, fm.bo.add))
setMethod("+", signature(e1 = "fm", e2 = "ANY"), function(e1, e2)
		  fm.mapply2.fm.ANY(e1, e2, fm.bo.add))
setMethod("+", signature(e1 = "ANY", e2 = "fm"), function(e1, e2)
		  fm.mapply2.fm.ANY(e2, e1, fm.bo.add))
setMethod("+", signature(e1 = "fmV", e2 = "ANY"), function(e1, e2)
		  fm.mapply2.fmV.ANY(e1, e2, fm.bo.add))
setMethod("+", signature(e1 = "ANY", e2 = "fmV"), function(e1, e2)
		  fm.mapply2.ANY.fmV(e1, e2, fm.bo.add))

setMethod("-", signature(e1 = "fm", e2 = "fm"), function(e1, e2)
		  fm.mapply2.fm(e1, e2, fm.bo.sub))
setMethod("-", signature(e1 = "fmV", e2 = "fmV"), function(e1, e2)
		  fm.mapply2.fmV(e1, e2, fm.bo.sub))
setMethod("-", signature(e1 = "fm", e2 = "fmV"), function(e1, e2)
		  fm.mapply.col(e1, e2, fm.bo.sub))
setMethod("-", signature(e1 = "fmV", e2 = "fm"), function(e1, e2)
		  -fm.mapply.col(e2, e1, fm.bo.sub))
setMethod("-", signature(e1 = "fm", e2 = "matrix"), function(e1, e2)
		  fm.mapply2.fm.m(e1, e2, fm.bo.sub))
setMethod("-", signature(e1 = "matrix", e2 = "fm"), function(e1, e2)
		  fm.mapply2.m.fm(e1, e2, fm.bo.sub))
setMethod("-", signature(e1 = "fm", e2 = "ANY"), function(e1, e2)
		  fm.mapply2.fm.ANY(e1, e2, fm.bo.sub))
setMethod("-", signature(e1 = "ANY", e2 = "fm"), function(e1, e2)
		  -fm.mapply2.fm.ANY(e2, e1, fm.bo.sub))
setMethod("-", signature(e1 = "fmV", e2 = "ANY"), function(e1, e2)
		  fm.mapply2.fmV.ANY(e1, e2, fm.bo.sub))
setMethod("-", signature(e1 = "ANY", e2 = "fmV"), function(e1, e2)
		  fm.mapply2.ANY.fmV(e1, e2, fm.bo.sub))

setMethod("*", signature(e1 = "fm", e2 = "fm"), function(e1, e2)
		  fm.mapply2.fm(e1, e2, fm.bo.mul))
setMethod("*", signature(e1 = "fmV", e2 = "fmV"), function(e1, e2)
		  fm.mapply2.fmV(e1, e2, fm.bo.mul))
setMethod("*", signature(e1 = "fm", e2 = "fmV"), function(e1, e2)
		  fm.mapply.col(e1, e2, fm.bo.mul))
setMethod("*", signature(e1 = "fmV", e2 = "fm"), function(e1, e2)
		  fm.mapply.col(e2, e1, fm.bo.mul))
setMethod("*", signature(e1 = "fm", e2 = "matrix"), function(e1, e2)
		  fm.mapply2.fm.m(e1, e2, fm.bo.mul))
setMethod("*", signature(e1 = "matrix", e2 = "fm"), function(e1, e2)
		  fm.mapply2.m.fm(e1, e2, fm.bo.mul))
setMethod("*", signature(e1 = "fm", e2 = "ANY"), function(e1, e2)
		  fm.mapply2.fm.ANY(e1, e2, fm.bo.mul))
setMethod("*", signature(e1 = "ANY", e2 = "fm"), function(e1, e2)
		  fm.mapply2.fm.ANY(e2, e1, fm.bo.mul))
setMethod("*", signature(e1 = "fmV", e2 = "ANY"), function(e1, e2)
		  fm.mapply2.fmV.ANY(e1, e2, fm.bo.mul))
setMethod("*", signature(e1 = "ANY", e2 = "fmV"), function(e1, e2)
		  fm.mapply2.ANY.fmV(e1, e2, fm.bo.mul))

setMethod("/", signature(e1 = "fm", e2 = "fm"), function(e1, e2)
		  fm.mapply2.fm(e1, e2, fm.bo.div))
setMethod("/", signature(e1 = "fmV", e2 = "fmV"), function(e1, e2)
		  fm.mapply2.fmV(e1, e2, fm.bo.div))
setMethod("/", signature(e1 = "fm", e2 = "fmV"), function(e1, e2)
		  fm.mapply.col(e1, e2, fm.bo.div))
setMethod("/", signature(e1 = "fmV", e2 = "fm"), function(e1, e2) {
		  print("don't support fmV/fm")
		  NULL
		  })
setMethod("/", signature(e1 = "fm", e2 = "matrix"), function(e1, e2)
		  fm.mapply2.fm.m(e1, e2, fm.bo.div))
setMethod("/", signature(e1 = "matrix", e2 = "fm"), function(e1, e2)
		  fm.mapply2.m.fm(e1, e2, fm.bo.div))
setMethod("/", signature(e1 = "fm", e2 = "ANY"), function(e1, e2)
		  fm.mapply2.fm.ANY(e1, e2, fm.bo.div))
setMethod("/", signature(e1 = "ANY", e2 = "fm"), function(e1, e2)
		  fm.mapply2.ANY.fm(e1, e2, fm.bo.div))
setMethod("/", signature(e1 = "fmV", e2 = "ANY"), function(e1, e2)
		  fm.mapply2.fmV.ANY(e1, e2, fm.bo.div))
setMethod("/", signature(e1 = "ANY", e2 = "fmV"), function(e1, e2)
		  fm.mapply2.ANY.fmV(e1, e2, fm.bo.div))

setMethod("==", signature(e1 = "fm", e2 = "fm"), function(e1, e2)
		  fm.mapply2.fm(e1, e2, fm.bo.eq))
setMethod("==", signature(e1 = "fmV", e2 = "fmV"), function(e1, e2)
		  fm.mapply2.fmV(e1, e2, fm.bo.eq))
setMethod("==", signature(e1 = "fm", e2 = "fmV"), function(e1, e2)
		  fm.mapply.col(e1, e2, fm.bo.eq))
setMethod("==", signature(e1 = "fmV", e2 = "fm"), function(e1, e2)
		  fm.mapply.col(e2, e1, fm.bo.eq))
setMethod("==", signature(e1 = "fm", e2 = "matrix"), function(e1, e2)
		  fm.mapply2.fm.m(e1, e2, fm.bo.eq))
setMethod("==", signature(e1 = "matrix", e2 = "fm"), function(e1, e2)
		  fm.mapply2.m.fm(e1, e2, fm.bo.eq))
setMethod("==", signature(e1 = "fm", e2 = "ANY"), function(e1, e2)
		  fm.mapply2.fm.ANY(e1, e2, fm.bo.eq))
setMethod("==", signature(e1 = "ANY", e2 = "fm"), function(e1, e2)
		  fm.mapply2.fm.ANY(e2, e1, fm.bo.eq))
setMethod("==", signature(e1 = "fmV", e2 = "ANY"), function(e1, e2)
		  fm.mapply2.fmV.ANY(e1, e2, fm.bo.eq))
setMethod("==", signature(e1 = "ANY", e2 = "fmV"), function(e1, e2)
		  fm.mapply2.ANY.fmV(e1, e2, fm.bo.eq))

setMethod("!=", signature(e1 = "fm", e2 = "fm"), function(e1, e2)
		  fm.mapply2.fm(e1, e2, fm.bo.neq))
setMethod("!=", signature(e1 = "fmV", e2 = "fmV"), function(e1, e2)
		  fm.mapply2.fmV(e1, e2, fm.bo.neq))
setMethod("!=", signature(e1 = "fm", e2 = "fmV"), function(e1, e2)
		  fm.mapply.col(e1, e2, fm.bo.neq))
setMethod("!=", signature(e1 = "fmV", e2 = "fm"), function(e1, e2)
		  fm.mapply.col(e2, e1, fm.bo.neq))
setMethod("!=", signature(e1 = "fm", e2 = "matrix"), function(e1, e2)
		  fm.mapply2.fm.m(e1, e2, fm.bo.neq))
setMethod("!=", signature(e1 = "matrix", e2 = "fm"), function(e1, e2)
		  fm.mapply2.m.fm(e1, e2, fm.bo.neq))
setMethod("!=", signature(e1 = "fm", e2 = "ANY"), function(e1, e2)
		  fm.mapply2.fm.ANY(e1, e2, fm.bo.neq))
setMethod("!=", signature(e1 = "ANY", e2 = "fm"), function(e1, e2)
		  fm.mapply2.fm.ANY(e2, e1, fm.bo.neq))
setMethod("!=", signature(e1 = "fmV", e2 = "ANY"), function(e1, e2)
		  fm.mapply2.fmV.ANY(e1, e2, fm.bo.neq))
setMethod("!=", signature(e1 = "ANY", e2 = "fmV"), function(e1, e2)
		  fm.mapply2.ANY.fmV(e1, e2, fm.bo.neq))

setMethod(">", signature(e1 = "fm", e2 = "fm"), function(e1, e2)
		  fm.mapply2.fm(e1, e2, fm.bo.gt))
setMethod(">", signature(e1 = "fmV", e2 = "fmV"), function(e1, e2)
		  fm.mapply2.fmV(e1, e2, fm.bo.gt))
setMethod(">", signature(e1 = "fm", e2 = "fmV"), function(e1, e2)
		  fm.mapply.col(e1, e2, fm.bo.gt))
setMethod(">", signature(e1 = "fmV", e2 = "fm"), function(e1, e2)
		  fm.mapply.col(e2, e1, fm.bo.le))
setMethod(">", signature(e1 = "fm", e2 = "matrix"), function(e1, e2)
		  fm.mapply2.fm.m(e1, e2, fm.bo.gt))
setMethod(">", signature(e1 = "matrix", e2 = "fm"), function(e1, e2)
		  fm.mapply2.m.fm(e1, e2, fm.bo.gt))
setMethod(">", signature(e1 = "fm", e2 = "ANY"), function(e1, e2)
		  fm.mapply2.fm.ANY(e1, e2, fm.bo.gt))
setMethod(">", signature(e1 = "ANY", e2 = "fm"), function(e1, e2)
		  fm.mapply2.fm.ANY(e2, e1, fm.bo.le))
setMethod(">", signature(e1 = "fmV", e2 = "ANY"), function(e1, e2)
		  fm.mapply2.fmV.ANY(e1, e2, fm.bo.gt))
setMethod(">", signature(e1 = "ANY", e2 = "fmV"), function(e1, e2)
		  fm.mapply2.ANY.fmV(e1, e2, fm.bo.gt))

setMethod(">=", signature(e1 = "fm", e2 = "fm"), function(e1, e2)
		  fm.mapply2.fm(e1, e2, fm.bo.ge))
setMethod(">=", signature(e1 = "fmV", e2 = "fmV"), function(e1, e2)
		  fm.mapply2.fmV(e1, e2, fm.bo.ge))
setMethod(">=", signature(e1 = "fm", e2 = "fmV"), function(e1, e2)
		  fm.mapply.col(e1, e2, fm.bo.ge))
setMethod(">=", signature(e1 = "fmV", e2 = "fm"), function(e1, e2)
		  fm.mapply.col(e2, e1, fm.bo.lt))
setMethod(">=", signature(e1 = "fm", e2 = "matrix"), function(e1, e2)
		  fm.mapply2.fm.m(e1, e2, fm.bo.ge))
setMethod(">=", signature(e1 = "matrix", e2 = "fm"), function(e1, e2)
		  fm.mapply2.m.fm(e1, e2, fm.bo.ge))
setMethod(">=", signature(e1 = "fm", e2 = "ANY"), function(e1, e2)
		  fm.mapply2.fm.ANY(e1, e2, fm.bo.ge))
setMethod(">=", signature(e1 = "ANY", e2 = "fm"), function(e1, e2)
		  fm.mapply2.fm.ANY(e2, e1, fm.bo.lt))
setMethod(">=", signature(e1 = "fmV", e2 = "ANY"), function(e1, e2)
		  fm.mapply2.fmV.ANY(e1, e2, fm.bo.ge))
setMethod(">=", signature(e1 = "ANY", e2 = "fmV"), function(e1, e2)
		  fm.mapply2.ANY.fmV(e1, e2, fm.bo.ge))

setMethod("<=", signature(e1 = "fm", e2 = "fm"), function(e1, e2)
		  fm.mapply2.fm(e1, e2, fm.bo.le))
setMethod("<=", signature(e1 = "fmV", e2 = "fmV"), function(e1, e2)
		  fm.mapply2.fmV(e1, e2, fm.bo.le))
setMethod("<=", signature(e1 = "fm", e2 = "fmV"), function(e1, e2)
		  fm.mapply.col(e1, e2, fm.bo.le))
setMethod("<=", signature(e1 = "fmV", e2 = "fm"), function(e1, e2)
		  fm.mapply.col(e2, e1, fm.bo.gt))
setMethod("<=", signature(e1 = "fm", e2 = "matrix"), function(e1, e2)
		  fm.mapply2.fm.m(e1, e2, fm.bo.le))
setMethod("<=", signature(e1 = "matrix", e2 = "fm"), function(e1, e2)
		  fm.mapply2.m.fm(e1, e2, fm.bo.le))
setMethod("<=", signature(e1 = "fm", e2 = "ANY"), function(e1, e2)
		  fm.mapply2.fm.ANY(e1, e2, fm.bo.le))
setMethod("<=", signature(e1 = "ANY", e2 = "fm"), function(e1, e2)
		  fm.mapply2.fm.ANY(e2, e1, fm.bo.gt))
setMethod("<=", signature(e1 = "fmV", e2 = "ANY"), function(e1, e2)
		  fm.mapply2.fmV.ANY(e1, e2, fm.bo.le))
setMethod("<=", signature(e1 = "ANY", e2 = "fmV"), function(e1, e2)
		  fm.mapply2.ANY.fmV(e1, e2, fm.bo.le))

setMethod("<", signature(e1 = "fm", e2 = "fm"), function(e1, e2)
		  fm.mapply2.fm(e1, e2, fm.bo.lt))
setMethod("<", signature(e1 = "fmV", e2 = "fmV"), function(e1, e2)
		  fm.mapply2.fmV(e1, e2, fm.bo.lt))
setMethod("<", signature(e1 = "fm", e2 = "fmV"), function(e1, e2)
		  fm.mapply.col(e1, e2, fm.bo.lt))
setMethod("<", signature(e1 = "fmV", e2 = "fm"), function(e1, e2)
		  fm.mapply.col(e2, e1, fm.bo.ge))
setMethod("<", signature(e1 = "fm", e2 = "matrix"), function(e1, e2)
		  fm.mapply2.fm.m(e1, e2, fm.bo.lt))
setMethod("<", signature(e1 = "matrix", e2 = "fm"), function(e1, e2)
		  fm.mapply2.m.fm(e1, e2, fm.bo.lt))
setMethod("<", signature(e1 = "fm", e2 = "ANY"), function(e1, e2)
		  fm.mapply2.fm.ANY(e1, e2, fm.bo.lt))
setMethod("<", signature(e1 = "ANY", e2 = "fm"), function(e1, e2)
		  fm.mapply2.fm.ANY(e2, e1, fm.bo.ge))
setMethod("<", signature(e1 = "fmV", e2 = "ANY"), function(e1, e2)
		  fm.mapply2.fmV.ANY(e1, e2, fm.bo.lt))
setMethod("<", signature(e1 = "ANY", e2 = "fmV"), function(e1, e2)
		  fm.mapply2.ANY.fmV(e1, e2, fm.bo.lt))

setMethod("|", signature(e1 = "fm", e2 = "fm"), function(e1, e2)
		  fm.mapply2.fm(e1, e2, fm.bo.or))
setMethod("|", signature(e1 = "fmV", e2 = "fmV"), function(e1, e2)
		  fm.mapply2.fmV(e1, e2, fm.bo.or))
setMethod("|", signature(e1 = "fm", e2 = "fmV"), function(e1, e2)
		  fm.mapply.col(e1, e2, fm.bo.or))
setMethod("|", signature(e1 = "fmV", e2 = "fm"), function(e1, e2)
		  fm.mapply.col(e2, e1, fm.bo.or))
setMethod("|", signature(e1 = "fm", e2 = "matrix"), function(e1, e2)
		  fm.mapply2.fm.m(e1, e2, fm.bo.or))
setMethod("|", signature(e1 = "matrix", e2 = "fm"), function(e1, e2)
		  fm.mapply2.m.fm(e1, e2, fm.bo.or))
setMethod("|", signature(e1 = "fm", e2 = "ANY"), function(e1, e2)
		  fm.mapply2.fm.ANY(e1, e2, fm.bo.or))
setMethod("|", signature(e1 = "ANY", e2 = "fm"), function(e1, e2)
		  fm.mapply2.fm.ANY(e2, e1, fm.bo.or))
setMethod("|", signature(e1 = "fmV", e2 = "ANY"), function(e1, e2)
		  fm.mapply2.fmV.ANY(e1, e2, fm.bo.or))
setMethod("|", signature(e1 = "ANY", e2 = "fmV"), function(e1, e2)
		  fm.mapply2.ANY.fmV(e1, e2, fm.bo.or))

setMethod("&", signature(e1 = "fm", e2 = "fm"), function(e1, e2)
		  fm.mapply2.fm(e1, e2, fm.bo.and))
setMethod("&", signature(e1 = "fmV", e2 = "fmV"), function(e1, e2)
		  fm.mapply2.fmV(e1, e2, fm.bo.and))
setMethod("&", signature(e1 = "fm", e2 = "fmV"), function(e1, e2)
		  fm.mapply.col(e1, e2, fm.bo.and))
setMethod("&", signature(e1 = "fmV", e2 = "fm"), function(e1, e2)
		  fm.mapply.col(e2, e1, fm.bo.and))
setMethod("&", signature(e1 = "fm", e2 = "matrix"), function(e1, e2)
		  fm.mapply2.fm.m(e1, e2, fm.bo.and))
setMethod("&", signature(e1 = "matrix", e2 = "fm"), function(e1, e2)
		  fm.mapply2.m.fm(e1, e2, fm.bo.and))
setMethod("&", signature(e1 = "fm", e2 = "ANY"), function(e1, e2)
		  fm.mapply2.fm.ANY(e1, e2, fm.bo.and))
setMethod("&", signature(e1 = "ANY", e2 = "fm"), function(e1, e2)
		  fm.mapply2.fm.ANY(e2, e1, fm.bo.and))
setMethod("&", signature(e1 = "fmV", e2 = "ANY"), function(e1, e2)
		  fm.mapply2.fmV.ANY(e1, e2, fm.bo.and))
setMethod("&", signature(e1 = "ANY", e2 = "fmV"), function(e1, e2)
		  fm.mapply2.ANY.fmV(e1, e2, fm.bo.and))

setMethod("-", signature(e1 = "fm", e2 = "missing"), function(e1)
		  fm.sapply.fm(e1, fm.buo.neg))
setMethod("-", signature(e1 = "fmV", e2 = "missing"), function(e1)
		  fm.sapply.fmV(e1, fm.buo.neg))

pmax2 <- function(e1, e2) pmax(e1, e2)
pmin2 <- function(e1, e2) pmin(e1, e2)
setGeneric("pmax2")
setGeneric("pmin2")

setMethod("pmax2", signature(e1 = "fm", e2 = "fm"), function(e1, e2)
		  fm.mapply2.fm(e1, e2, fm.bo.max))
setMethod("pmax2", signature(e1 = "fmV", e2 = "fmV"), function(e1, e2)
		  fm.mapply2.fmV(e1, e2, fm.bo.max))
setMethod("pmax2", signature(e1 = "fm", e2 = "fmV"), function(e1, e2)
		  fm.mapply.col(e1, e2, fm.bo.max))
setMethod("pmax2", signature(e1 = "fmV", e2 = "fm"), function(e1, e2)
		  fm.mapply.col(e2, e1, fm.bo.max))
setMethod("pmax2", signature(e1 = "fm", e2 = "matrix"), function(e1, e2)
		  fm.mapply2.fm.m(e1, e2, fm.bo.max))
setMethod("pmax2", signature(e1 = "matrix", e2 = "fm"), function(e1, e2)
		  fm.mapply2.m.fm(e1, e2, fm.bo.max))
setMethod("pmax2", signature(e1 = "fm", e2 = "ANY"), function(e1, e2)
		  fm.mapply2.fm.ANY(e1, e2, fm.bo.max))
setMethod("pmax2", signature(e1 = "ANY", e2 = "fm"), function(e1, e2)
		  fm.mapply2.fm.ANY(e2, e1, fm.bo.max))
setMethod("pmax2", signature(e1 = "fmV", e2 = "ANY"), function(e1, e2)
		  fm.mapply2.fmV.ANY(e1, e2, fm.bo.max))
setMethod("pmax2", signature(e1 = "ANY", e2 = "fmV"), function(e1, e2)
		  fm.mapply2.ANY.fmV(e1, e2, fm.bo.max))

setMethod("pmin2", signature(e1 = "fm", e2 = "fm"), function(e1, e2)
		  fm.mapply2.fm(e1, e2, fm.bo.min))
setMethod("pmin2", signature(e1 = "fmV", e2 = "fmV"), function(e1, e2)
		  fm.mapply2.fmV(e1, e2, fm.bo.min))
setMethod("pmin2", signature(e1 = "fm", e2 = "fmV"), function(e1, e2)
		  fm.mapply.col(e1, e2, fm.bo.min))
setMethod("pmin2", signature(e1 = "fmV", e2 = "fm"), function(e1, e2)
		  fm.mapply.col(e2, e1, fm.bo.min))
setMethod("pmin2", signature(e1 = "fm", e2 = "matrix"), function(e1, e2)
		  fm.mapply2.fm.m(e1, e2, fm.bo.min))
setMethod("pmin2", signature(e1 = "matrix", e2 = "fm"), function(e1, e2)
		  fm.mapply2.m.fm(e1, e2, fm.bo.min))
setMethod("pmin2", signature(e1 = "fm", e2 = "ANY"), function(e1, e2)
		  fm.mapply2.fm.ANY(e1, e2, fm.bo.min))
setMethod("pmin2", signature(e1 = "ANY", e2 = "fm"), function(e1, e2)
		  fm.mapply2.fm.ANY(e2, e1, fm.bo.min))
setMethod("pmin2", signature(e1 = "fmV", e2 = "ANY"), function(e1, e2)
		  fm.mapply2.fmV.ANY(e1, e2, fm.bo.min))
setMethod("pmin2", signature(e1 = "ANY", e2 = "fmV"), function(e1, e2)
		  fm.mapply2.ANY.fmV(e1, e2, fm.bo.min))

`!.fm` <- function(o)
{
	fm.sapply.fm(o, fm.buo.not)
}

`!.fmV` <- function(o)
{
	fm.sapply.fmV(o, fm.buo.not)
}

setMethod("%*%", signature(x = "fm", y = "fm"), function(x, y) fm.multiply(x, y))
setMethod("%*%", signature(x = "fm", y = "fmV"), function(x, y) fm.multiply(x, y))

setMethod("dim", signature(x = "fm"), function(x) c(x@nrow, x@ncol))
setMethod("nrow", signature(x = "fm"), function(x) x@nrow)
setMethod("ncol", signature(x = "fm"), function(x) x@ncol)
setMethod("length", signature(x = "fmV"), function(x) x@len)
setMethod("length", signature(x = "fm"), function(x) x@nrow * x@ncol)
setMethod("typeof", signature(x = "fm"), function(x) fm.typeof(x))
setMethod("typeof", signature(x = "fmV"), function(x) fm.typeof(x))
setMethod("t", signature(x = "fm"), function(x) fm.t(x))

fm.cls <- list("fm", "fmV")

setGeneric("pmax", signature="...")
setGeneric("pmin", signature="...")

fm.mapply.list <- function(data, FUN, test.na)
{
	n <- length(data)
	res <- data[[1]]
	if (n > 1) {
		for (arg in data[2:n]) {
			if (class(res) == class(arg))
				res <- fm.mapply2(res, arg, FUN, test.na)
			else if (class(res) == "fm" && class(arg) == "fmV")
				res <- fm.mapply.col(res, arg, FUN, test.na)
			# We assume that FUN is commutative.
			# This function is used for pmin/pmax, so it's fine.
			else if (class(res) == "fmV" && class(arg) == "fm")
				res <- fm.mapply.col(arg, res, FUN, test.na)
			else
				stop("unknown arguments")
		}
	}
	res
}

get.zero <- function(type)
{
	if (type == "double")
		0
	else if (type == "integer")
		as.integer(0)
	else if (type == "logical")
		FALSE
	else
		stop("only integer and numeric value have 0")
}

get.max.val <- function(type)
{
	if (type == "double")
		.Machine$double.xmax
	else if (type == "integer")
		.Machine$integer.max
	else
		TRUE
}

get.min.val <- function(type)
{
	if (type == "double")
		-.Machine$double.xmax
	else if (type == "integer")
		-.Machine$integer.max
	else
		FALSE
}

get.na <- function(type)
{
	if (type == "double")
		as.double(NA)
	else if (type == "integer")
		as.integer(NA)
	else
		NA
}

fm.range1 <- function(x, na.rm)
{
	test.na <- TRUE
	if (na.rm) {
		x.min <- ifelse(is.na(x), get.min.val(typeof(x)), x)
		x.max <- ifelse(is.na(x), get.max.val(typeof(x)), x)
		test.na <- FALSE
		tmp1 <- fm.agg.lazy(x.max, fm.bo.min)
		tmp2 <- fm.agg.lazy(x.min, fm.bo.max)
	}
	else {
		tmp1 <- fm.agg.lazy(x, fm.bo.min)
		tmp2 <- fm.agg.lazy(x, fm.bo.max)
	}
	if (test.na) {
		x.is.na <- fm.agg.lazy(fm.is.na.only(x), fm.bo.or)
		res <- fm.materialize(tmp1, tmp2, x.is.na)
		if (fmV2scalar(res[[3]])) {
			na <- get.na(typeof(x))
			c(na, na)
		}
		else
			c(fmV2scalar(res[[1]]), fmV2scalar(res[[2]]))
	}
	else {
		res <- fm.materialize(tmp1, tmp2)
		c(fmV2scalar(res[[1]]), fmV2scalar(res[[2]]))
	}
}

# We replace both NA and NaN
replace.na <- function(obj, val)
{
	ifelse(is.na(obj), val, obj)
}
replace.na.list <- function(objs, val)
{
	if (length(objs) == 0)
		return(objs)

	for (i in 1:length(objs))
		objs[[i]] <- ifelse(is.na(objs[[i]]), val, objs[[i]])
	objs
}

fm.agg.na <- function(fm, op, test.na)
{
	if (test.na && fm.env$fm.test.na) {
		any.na <- fm.agg.lazy(fm.is.na.only(fm), fm.bo.or)
		agg.res <- fm.agg.lazy(fm, op)
		res <- fm.materialize(any.na, agg.res)
		if (fmV2scalar(res[[1]]))
			get.na(typeof(agg.res))
		else
			fmV2scalar(res[[2]])
	}
	else
		fm.agg(fm, op)
}

# Aggregation on a FlashMatrixR vector/matrix.
for (cl in fm.cls) {
	# TODO we need to handle na.rm for all of the functions here properly.
	setMethod("sum", cl, function(x, ..., na.rm) {
			  others <- list(...)
			  test.na <- TRUE
			  if (na.rm) {
				  zero <- get.zero(typeof(x))
				  x <- replace.na(x, zero)
				  others <- replace.na.list(others, zero)
				  test.na <- FALSE
			  }
			  res <- fm.agg.na(x, fm.bo.add, test.na)
			  if (length(others) >= 1) {
				  for (arg in others)
					  res <- res + fm.agg.na(arg, fm.bo.add, test.na)
			  }
			  res
		  })
	setMethod("min", cl, function(x, ..., na.rm) {
			  others <- list(...)
			  test.na <- TRUE
			  if (na.rm) {
				  max.val <- get.max.val(typeof(x))
				  x <- replace.na(x, max.val)
				  others <- replace.na.list(others, max.val)
				  test.na <- FALSE
			  }
			  res <- fm.agg.na(x, fm.bo.min, test.na)
			  if (length(others) >= 1) {
				  for (arg in others)
					  res <- min(res, fm.agg.na(arg, fm.bo.min, test.na))
			  }
			  res
		  })
	setMethod("max", cl, function(x, ..., na.rm) {
			  others <- list(...)
			  test.na <- TRUE
			  if (na.rm) {
				  min.val <- get.min.val(typeof(x))
				  x <- replace.na(x, min.val)
				  others <- replace.na.list(others, min.val)
				  test.na <- FALSE
			  }
			  res <- fm.agg.na(x, fm.bo.max, test.na)
			  if (length(others) >= 1) {
				  for (arg in others)
					  res <- max(res, fm.agg.na(arg, fm.bo.max, test.na))
			  }
			  res
		  })
	setMethod("range", cl, function(x, ..., na.rm) {
			  others <- list(...)
			  res <- fm.range1(x, na.rm)
			  if (length(others) >= 1) {
				  for (arg in others) {
					  tmp <- fm.range1(arg, na.rm)
					  # If there is NA in the data, there is no point of
					  # continuing.
					  if (is.na(tmp))
						  return(tmp)
					  res[1] <- min(res[1], tmp[1])
					  res[2] <- max(res[2], tmp[2])
				  }
			  }
			  res
		  })
	setMethod("pmin", cl, function(..., na.rm = FALSE) {
			  args <- list(...)
			  if (length(args) == 0)
				  stop("no arguments")
			  if (na.rm) {
				  args <- replace.na.list(args, get.max.val(typeof(args[[1]])))
				  fm.mapply.list(args, fm.bo.min, FALSE)
			  }
			  else
				  fm.mapply.list(args, fm.bo.min, TRUE)
		  })
	setMethod("pmax", cl, function(..., na.rm = FALSE) {
			  args <- list(...)
			  if (length(args) == 0)
				  stop("no arguments")
			  if (na.rm) {
				  args <- replace.na.list(args, get.min.val(typeof(args[[1]])))
				  fm.mapply.list(args, fm.bo.min, FALSE)
			  }
			  else
				  fm.mapply.list(args, fm.bo.max, TRUE)
		  })
	setMethod("sd", cl, function(x, na.rm) {
			  n <- length(x)
			  x2 <- x * x
			  test.na <- TRUE
			  num.na <- 0
			  if (na.rm) {
				  zero <- get.zero(typeof(x))
				  in.is.na <- is.na(x)
				  x <- ifelse(in.is.na, zero, x)
				  x2 <- ifelse(in.is.na, zero, x2)
				  test.na <- FALSE
			  }
			  sum.x <- fm.agg.lazy(x, fm.bo.add)
			  sum.x2 <- fm.agg.lazy(x2, fm.bo.add)

			  if (test.na) {
				  x.is.na <- fm.agg.lazy(fm.is.na.only(x), fm.bo.or)
				  res <- fm.materialize(sum.x, sum.x2, x.is.na)
				  if (fmV2scalar(res[[3]]))
					  return(get.na(typeof(x)))
				  sums <- res[1:2]
			  }
			  else {
				  # If we remove NA, we should calculate the number of
				  # NAs in the vector.
				  sum.na <- fm.agg.lazy(in.is.na, fm.bo.add)
				  res <- fm.materialize(sum.x, sum.x2, sum.na)
				  n <- n - fmV2scalar(res[[3]])
				  sums <- res[1:2]
			  }
			  sum.x <- fmV2scalar(sums[[1]])
			  sum.x2 <- fmV2scalar(sums[[2]])
			  avg <- sum.x / n
			  sqrt((sum.x2 - n * avg * avg) / (n - 1))
		  })
	# TODO I need to implemented trimmed mean
	setMethod("mean", cl, function(x, ...) {
			  args <- list(...)
			  na.rm <- FALSE
			  if (!is.null(args[["na.rm"]]))
				  na.rm <- args[["na.rm"]]
			  n <- length(x)
			  # TODO I need to lazily calculate it.
			  if (na.rm)
				  n <- n - sum(is.na(x))
			  sum(x, na.rm=na.rm) / n
		  })
}

# Miscellaneous Mathematical Functions
setMethod("abs", signature(x = "fm"), function(x) fm.sapply.fm(x, fm.buo.abs))
setMethod("abs", signature(x = "fmV"), function(x) fm.sapply.fmV(x, fm.buo.abs))
setMethod("sqrt", signature(x = "fm"), function(x) fm.sapply.fm(x, fm.buo.sqrt))
setMethod("sqrt", signature(x = "fmV"), function(x) fm.sapply.fmV(x, fm.buo.sqrt))
setMethod("ceiling", signature(x = "fm"), function(x) fm.sapply.fm(x, fm.buo.ceil))
setMethod("ceiling", signature(x = "fmV"), function(x) fm.sapply.fmV(x, fm.buo.ceil))
setMethod("floor", signature(x = "fm"), function(x) fm.sapply.fm(x, fm.buo.floor))
setMethod("floor", signature(x = "fmV"), function(x) fm.sapply.fmV(x, fm.buo.floor))
setMethod("round", signature(x = "fm"), function(x) fm.sapply.fm(x, fm.buo.round))
setMethod("round", signature(x = "fmV"), function(x) fm.sapply.fmV(x, fm.buo.round))
setMethod("log10", signature(x = "fm"), function(x) fm.sapply.fm(x, fm.buo.log10))
setMethod("log10", signature(x = "fmV"), function(x) fm.sapply.fmV(x, fm.buo.log10))
setMethod("log2", signature(x = "fm"), function(x) fm.sapply.fm(x, fm.buo.log2))
setMethod("log2", signature(x = "fmV"), function(x) fm.sapply.fmV(x, fm.buo.log2))
setMethod("exp", signature(x = "fm"), function(x) fm.mapply2(exp(1), x, fm.bo.pow, TRUE))
setMethod("exp", signature(x = "fmV"), function(x) fm.mapply2(exp(1), x, fm.bo.pow, TRUE))
setMethod("log", "fm", function(x, base=exp(1)) {
		  if (base == exp(1))
			  fm.sapply.fm(x, fm.buo.log)
		  else
			  fm.sapply.fm(x, fm.buo.log) / log(base)
})
setMethod("log", "fmV", function(x, base=exp(1)) {
		  if (base == exp(1))
			  fm.sapply.fmV(x, fm.buo.log)
		  else
			  fm.sapply.fmV(x, fm.buo.log) / log(base)
})

# TODO I need to handle na.rm and dims here as well.
setMethod("rowSums", signature(x = "fm", na.rm = "ANY", dims = "ANY"),
		  function(x, na.rm, dims) {
			  fm.agg.mat(x, 1, fm.bo.add)
		  })
setMethod("colSums", signature(x = "fm", na.rm = "ANY", dims = "ANY"),
		  function(x, na.rm, dims) {
			  fm.agg.mat(x, 2, fm.bo.add)
		  })
setMethod("rowMeans", signature(x = "fm", na.rm = "ANY", dims = "ANY"),
		  function(x, na.rm, dims) {
			  rowSums(x, na.rm, dims) / dim(x)[2]
		  })

setMethod("colMeans", signature(x = "fm", na.rm = "ANY", dims = "ANY"),
		  function(x, na.rm, dims) {
			  colSums(x, na.rm, dims) / dim(x)[1]
		  })

# Print FlashMatrix objects.
setMethod("print", signature(x = "fm"), function(x)
		  cat("FlashMatrixR matrix ", x@name, ": ", dim(x)[1], " rows, ", dim(x)[2],
			  " columns, is sparse: ", fm.is.sparse(x), "\n", sep=""))
setMethod("print", signature(x = "fmV"), function(x)
	cat("FlashVectorR vector ", x@name, ": length: ", length(x), "\n", sep=""))
setMethod("print", signature(x = "fm.bo"), function(x)
	cat("FLashMatrixR basic operator:", x@name, "\n"))

fm.table <- function(x)
{
	count <- fm.create.agg.op(fm.bo.count, fm.bo.add, "count")
	fm.sgroupby(x, count)
}

fm.as.integer <- function(x) fm.sapply(x, fm.buo.as.int, TRUE)
fm.as.numeric <- function(x) fm.sapply(x, fm.buo.as.numeric, TRUE)

fmV2scalar <- function(x)
{
	fm.conv.FM2R(x)[1]
}
