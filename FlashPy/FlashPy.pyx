import numpy as np
import ctypes
# "cimport" is used to import special compile-time information
# about the numpy module (this is stored in a file numpy.pxd which is
# currently part of the Cython distribution).
cimport numpy as np
from libc.stdlib cimport free, malloc
from libc.stdint cimport intptr_t
from libcpp.vector cimport vector
from libcpp.string cimport string
from libcpp cimport bool
from libc.string cimport memcpy

from cpython cimport array
import array

np.import_array()

cdef enum bulk_op_idx_t:
    OP_ADD, OP_SUB, OP_MUL, OP_DIV
    OP_MIN, OP_MAX,
    OP_POW,
    OP_EQ, OP_NEQ, OP_GT, OP_GE, OP_LT, OP_LE,
    OP_OR, OP_AND,
    OP_MOD, OP_IDIV

bop_add = OP_ADD
bop_sub = OP_SUB
bop_mul = OP_MUL
bop_div = OP_DIV
bop_min = OP_MIN
bop_max = OP_MAX
bop_pow = OP_POW
bop_eq = OP_EQ
bop_neq = OP_NEQ
bop_gt = OP_GT
bop_ge = OP_GE
bop_lt = OP_LT
bop_le = OP_LE
bop_or = OP_OR
bop_and = OP_AND
bop_mod = OP_MOD
bop_idiv = OP_IDIV

cdef enum bulk_uop_idx_t:
    UOP_NEG, UOP_SQRT, UOP_ABS, UOP_NOT, UOP_SQ,
    UOP_CEIL, UOP_FLOOR, UOP_ROUND,
    UOP_LOG, UOP_LOG2, UOP_LOG10

uop_neg = UOP_NEG
uop_sqrt = UOP_SQRT
uop_abs = UOP_ABS
uop_not = UOP_NOT
uop_sq = UOP_SQ
uop_ceil = UOP_CEIL
uop_floor = UOP_FLOOR
uop_round = UOP_ROUND
uop_log = UOP_LOG
uop_log2 = UOP_LOG2
uop_log10 = UOP_LOG10

cdef extern from "MatrixWrapper.h" namespace "flashpy":
    cdef bool init_flashpy_c(const string file)

cdef extern from "MatrixWrapper.h" namespace "flashpy":
    cdef cppclass scalar_wrapper:
        scalar_wrapper()
        const char *get_raw() const;

cdef extern from "MatrixWrapper.h" namespace "flashpy":
    cdef scalar_wrapper create_scalar_wrapper[T](T x)

cdef extern from "MatrixWrapper.h" namespace "flashpy":
    cdef cppclass matrix_wrapper:
        matrix_wrapper()
        # create a vector with data from "data_addr".
        matrix_wrapper(intptr_t data_addr, size_t length,
                const string &t)
        # create a matrix with data from "data_addr".
        matrix_wrapper(intptr_t data_addr, size_t nrow, size_t ncol,
                const string &t, const string layout)
        # create an empty vector with the specified size
        matrix_wrapper(size_t length, string t)
        # create an empty matrix with the specified size
        matrix_wrapper(size_t nrow, size_t ncol, string t, string layout)

        void init_seq[T](T start, T stride, bool byrow)
        void init_const_float(double val)
        void init_const_int(long val)

        void set_cached(bool)

        matrix_wrapper as_vector() const
        matrix_wrapper as_matrix() const

        size_t get_num_rows() const
        size_t get_num_cols() const
        size_t get_entry_size() const
        string get_type_str() const
        np.NPY_TYPES get_type_py() const
        string get_layout() const

        bool is_floating_point() const
        bool is_in_mem() const
        bool is_virtual() const
        bool is_vector() const

        bool materialize_self() const
        matrix_wrapper cast_ele_type(string dtyp) const

        matrix_wrapper get_col(long idx) const
        matrix_wrapper get_row(long idx) const
        matrix_wrapper get_cols(const vector[long] &idxs) const
        matrix_wrapper get_rows(const vector[long] &idxs) const
        matrix_wrapper get_cols(long start, long stop, long step) const
        matrix_wrapper get_rows(long start, long stop, long step) const
        matrix_wrapper get_cols(matrix_wrapper idxs) const
        matrix_wrapper get_rows(matrix_wrapper idxs) const

        matrix_wrapper set_cols(const vector[long] &idxs, matrix_wrapper cols)
        matrix_wrapper set_rows(const vector[long] &idxs, matrix_wrapper rows)
        bool copy_rows_to(char *arr, size_t len) const
        matrix_wrapper transpose() const
        matrix_wrapper conv_store(bool in_mem, int num_nodes) const
        matrix_wrapper conv_layout(const string layout) const
        matrix_wrapper inner_prod(matrix_wrapper m, bulk_op_idx_t left_op,
                bulk_op_idx_t right_op) const
        matrix_wrapper multiply(matrix_wrapper m) const
        matrix_wrapper aggregate(bulk_op_idx_t op)
        matrix_wrapper agg_row(bulk_op_idx_t op) const
        matrix_wrapper agg_col(bulk_op_idx_t op) const
        matrix_wrapper groupby_row(matrix_wrapper labels, bulk_op_idx_t op) const
        matrix_wrapper groupby_row(matrix_wrapper labels, bulk_op_idx_t op) const
        matrix_wrapper mapply_cols(matrix_wrapper vals, bulk_op_idx_t op) const
        matrix_wrapper mapply_rows(matrix_wrapper vals, bulk_op_idx_t op) const
        matrix_wrapper mapply2(matrix_wrapper m, bulk_op_idx_t op) const
        matrix_wrapper sapply(bulk_uop_idx_t op) const
        matrix_wrapper apply_scalar(scalar_wrapper var, bulk_op_idx_t op) const
        matrix_wrapper ifelse(matrix_wrapper x, matrix_wrapper y) const

class flagsobj:
    def __init__(self):
        self.c_contiguous = False
        self.f_contiguous = False
        self.owndata = True
        self.writable = False
        self.aligned = True
        self.updateifcopy = False

    def set_layout(self, layout):
        if (layout == "C"):
            self.c_contiguous = True
            self.f_contiguous = False
        elif (layout == "F"):
            self.c_contiguous = False
            self.f_contiguous = True
        else:
            raise ValueError("Invalid layout")

cdef class PyMatrix:
    cdef matrix_wrapper mat      # hold a C++ instance which we're wrapping
    cdef readonly int ndim
    cdef readonly object shape
    cdef readonly string dtype
    cdef readonly object flags
    cdef readonly long size
    cdef readonly long itemsize
    cdef readonly long nbytes
    cdef readonly PyMatrix T

    def __cinit__(self):
        self.mat = matrix_wrapper()

    def __init__(self):
        self.ndim = 0
        self.shape = (0, 0)
        self.flags = flagsobj()

    def __array__(self):
        cdef np.npy_intp shape[2]
        shape[0] = self.shape[0]
        if (self.ndim >= 2):
            shape[1] = self.shape[1]
        tmp = np.PyArray_SimpleNew(self.ndim, shape, self.mat.get_type_py())
        self.mat.copy_rows_to(np.PyArray_BYTES(tmp),
                self.mat.get_num_rows() * self.mat.get_num_cols() * self.mat.get_entry_size())
        return tmp

    # Special Methods Table
    # http://cython.readthedocs.io/en/latest/src/reference/special_methods_table.html

    def __richcmp__(PyMatrix x, y, int op):
        cdef PyMatrix ret = PyMatrix()
        # Rich comparisons:
        # http://cython.readthedocs.io/en/latest/src/userguide/special_methods.html#rich-comparisons
        # <   0
        # ==  2
        # >   4
        # <=  1
        # !=  3
        # >=  5
        if (op == 0):
            ret = x.mapply2(y, OP_LT)
        elif (op == 2):
            ret = x.mapply2(y, OP_EQ)
        elif (op == 4):
            ret = x.mapply2(y, OP_GT)
        elif (op == 1):
            ret = x.mapply2(y, OP_LE)
        elif (op == 3):
            ret = x.mapply2(y, OP_NEQ)
        elif (op == 5):
            ret = x.mapply2(y, OP_GE)
        else:
            raise ValueError("invalid argument")
        return ret

    def __add__(PyMatrix x, y):
        return x.mapply2(y, OP_ADD)

    def __sub__(PyMatrix x, y):
        return x.mapply2(y, OP_SUB)

    def __mul__(PyMatrix x, y):
        return x.mapply2(y, OP_MUL)

    def __div__(PyMatrix x, y):
        return x.mapply2(y, OP_DIV)

    def __floordiv__(PyMatrix x, y):
        return x.mapply2(y, OP_IDIV)

    def __mod__(PyMatrix x, PyMatrix y):
        return x.mapply2(y, OP_MOD)

    def __and__(PyMatrix x, PyMatrix y):
        return x.mapply2(y, OP_AND)

    def __or__(PyMatrix x, PyMatrix y):
        return x.mapply2(y, OP_OR)

    def __neg__(self):
        return self.sapply(UOP_NEG)

    def __abs__(self):
        return self.sapply(UOP_ABS)

    def __len__(self):
        return self.mat.get_num_rows()

    def __getitem__(self, key):
        cdef PyMatrix ret = PyMatrix()
        if (isinstance(key, tuple) and len(key) >= 2):
            if (len(key) > 2):
                raise IndexError("too many indices for array")
            if (self.ndim == 1 or self.shape[0] > self.shape[1]):
                ret = self.get_cols(key[1])
                ret = ret.get_rows(key[0])
            else:
                ret = self.get_rows(key[0])
                ret = self.get_cols(key[1])
        else:
            ret = self.get_rows(key)
        return ret

    def init_attr(self, T=None):
        if (self.mat.is_vector()):
            self.shape = (self.mat.get_num_rows(),)
            self.ndim = 1
        else:
            self.shape = (self.mat.get_num_rows(), self.mat.get_num_cols())
            self.ndim = 2
        self.dtype = self.mat.get_type_str()
        self.flags.set_layout(self.mat.get_layout())
        self.size = self.mat.get_num_rows() * self.mat.get_num_cols()
        self.itemsize = self.mat.get_entry_size()
        self.nbytes = self.size * self.itemsize
        if (self.ndim < 2):
            self.T = self
        elif (T is None):
            self.T = self.transpose()
        else:
            self.T = T

    # These are functions in numpy

    def astype(self, string dtype, order='K', casting='unsafe', subok=True, copy=True):
        # TODO I need to check casting
        arr = self.cast_ele_type(dtype)
        if (order == 'C' or order == 'F'):
            arr = arr.conv_layout(order)
        return arr

    def dot(self, b, out=None):
        cdef PyMatrix res = PyMatrix()
        if (self.ndim == 1 and b.ndim == 1):
            res = self.as_matrix().transpose().multiply(b)
        else:
            res = self.multiply(b)
        if (out is not None):
            out.assign(res)
        return res

    def sum(self, axis=None, dtype=None, out=None, keepdims=False):
        return self.aggregate_(OP_ADD, axis, dtype, out, keepdims)

    def prod(self, axis=None, dtype=None, out=None, keepdims=False):
        return self.aggregate_(OP_MUL, axis, dtype, out, keepdims)

    def mean(self, axis=None, dtype=None, out=None, keepdims=False):
        cdef PyMatrix a = self
        if (not a.mat.is_floating_point()):
            a = a.cast_ele_type("d")
        s = a.sum(axis, dtype, out, keepdims)
        if axis is None:
            return s/a.mat.get_num_rows()/a.mat.get_num_cols()
        else:
            return s/a.shape[axis]

    # These are specific for FlashMatrix.

    def set_cached(self, cached):
        self.mat.set_cached(cached)

    def is_in_mem(self):
        return self.mat.is_in_mem()

    def is_virtual(self):
        return self.mat.is_virtual()

    def materialize_self(self):
        return self.mat.materialize_self()

    def get_rows(self, idxs):
        cdef PyMatrix ret = PyMatrix()
        cdef long *addr
        cdef vector[long] cidxs
        cdef array.array idx_arr
        if (np.isscalar(idxs)):
            ret.mat = self.mat.get_row(idxs)
        elif (isinstance(idxs, list)):
            cidxs = idxs
            ret.mat = self.mat.get_rows(cidxs)
        elif (isinstance(idxs, slice)):
            if (idxs.start is None and idxs.stop is None and idxs.step is None):
                ret.mat = self.mat
            elif (idxs.step is None):
                ret.mat = self.mat.get_rows(idxs.start, idxs.stop, 1)
            else:
                ret.mat = self.mat.get_rows(idxs.start, idxs.stop, idxs.step)
        elif (isinstance(idxs, np.ndarray)):
            idxs = np.array(idxs, dtype='l')
            addr = <long *>np.PyArray_DATA(idxs)
            cidxs.assign(addr, addr + len(idxs))
            ret.mat = self.mat.get_rows(cidxs)
        elif (idxs is None):
            if (self.ndim >= 2):
                raise IndexError("doesn't support high dimensional array")
            # If this is a vector, we return a one-row matrix.
            return self.as_matrix().transpose()
        else:
            raise ValueError("invalid index")
        ret.init_attr()
        return ret

    def get_cols(self, idxs):
        cdef PyMatrix ret = PyMatrix()
        cdef long *addr
        cdef vector[long] cidxs
        cdef array.array idx_arr
        if (np.isscalar(idxs)):
            ret.mat = self.mat.get_col(idxs)
        elif (isinstance(idxs, list)):
            cidxs = idxs
            ret.mat = self.mat.get_cols(cidxs)
        elif (isinstance(idxs, slice)):
            if (idxs.start is None and idxs.stop is None and idxs.step is None):
                ret.mat = self.mat
            elif (idxs.step is None):
                ret.mat = self.mat.get_cols(idxs.start, idxs.stop, 1)
            else:
                ret.mat = self.mat.get_cols(idxs.start, idxs.stop, idxs.step)
        elif (isinstance(idxs, np.ndarray)):
            idxs = np.array(idxs, dtype='l')
            addr = <long *>np.PyArray_DATA(idxs)
            cidxs.assign(addr, addr + len(idxs))
            ret.mat = self.mat.get_cols(cidxs)
        elif (idxs is None):
            if (self.ndim >= 2):
                raise IndexError("doesn't support high dimensional array")
            # If this is a vector, we return a one-col matrix.
            return self.as_matrix()
        else:
            raise ValueError("invalid index")
        ret.init_attr()
        return ret

    def transpose(self):
        if (self.ndim < 2):
            return self

        cdef PyMatrix ret = PyMatrix()
        ret.mat = self.mat.transpose()
        ret.init_attr(self)
        return ret

    def conv_store(self, bool in_mem, int num_nodes):
        cdef PyMatrix ret = PyMatrix()
        ret.mat = self.mat.conv_store(in_mem, num_nodes)
        ret.init_attr()
        return ret

    def conv_layout(self, string order):
        cdef PyMatrix ret = PyMatrix()
        ret.mat = self.mat.conv_layout(order)
        ret.init_attr()
        return ret

    def cast_ele_type(self, string dtype):
        cdef PyMatrix ret = PyMatrix()
        ret.mat = self.mat.cast_ele_type(dtype)
        ret.init_attr()
        return ret

    def assign(self, PyMatrix mat):
        self.mat = mat.mat
        self.init_attr()

    def multiply(self, obj):
        cdef PyMatrix ret = PyMatrix()
        cdef PyMatrix mat
        if (isinstance(obj, PyMatrix)):
            mat = <PyMatrix>obj
        else:
            mat = array(obj)
        ret.mat = self.mat.multiply(mat.mat)
        ret.init_attr()
        return ret

    def as_vector(self):
        cdef PyMatrix ret = PyMatrix()
        ret.mat = self.mat.as_vector()
        ret.init_attr()
        return ret

    def as_matrix(self):
        cdef PyMatrix ret = PyMatrix()
        ret.mat = self.mat.as_matrix()
        ret.init_attr()
        return ret

    def aggregate_(self, op, axis=None, dtype=None, out=None, keepdims=False):
        cdef PyMatrix ret = PyMatrix()
        cdef PyMatrix a = self
        if dtype is not None:
            a = a.cast_ele_type(dtype)
        if axis is None:
            ret = a.aggregate(op).as_vector()
        elif (axis == 0):
            ret = a.agg_col(op).as_vector()
        elif (axis == 1):
            ret = a.agg_row(op).as_vector()
        else:
            raise ValueError("invalid axis")
        # TODO let's ignore keepdims for now.
        if out is not None:
            out.assign(ret)
        return ret

    # These are generalized functions.

    def inner_prod(self, PyMatrix mat, left_op, right_op):
        cdef PyMatrix ret = PyMatrix()
        ret.mat = self.mat.inner_prod(mat.mat, left_op, right_op)
        ret.init_attr()
        return ret

    def aggregate(self, op):
        cdef PyMatrix ret = PyMatrix()
        ret.mat = self.mat.aggregate(op)
        ret.init_attr()
        return ret

    def agg_row(self, op):
        cdef PyMatrix ret = PyMatrix()
        ret.mat = self.mat.agg_row(op)
        ret.init_attr()
        return ret

    def agg_col(self, op):
        cdef PyMatrix ret = PyMatrix()
        ret.mat = self.mat.agg_col(op)
        ret.init_attr()
        return ret

    def mapply_rows(self, PyMatrix vec, op):
        cdef PyMatrix ret = PyMatrix()
        ret.mat = self.mat.mapply_rows(vec.mat, op)
        ret.init_attr()
        return ret

    def mapply_cols(self, PyMatrix vec, op):
        cdef PyMatrix ret = PyMatrix()
        ret.mat = self.mat.mapply_cols(vec.mat, op)
        ret.init_attr()
        return ret

    def mapply2(self, obj, op):
        cdef PyMatrix ret = PyMatrix()
        cdef scalar_wrapper var
        cdef PyMatrix mat
        if (isinstance(obj, PyMatrix)):
            mat = <PyMatrix>obj
            ret.mat = self.mat.mapply2(mat.mat, op)
        elif (np.isscalar(obj)):
            if (isinstance(obj, float)):
                var = create_scalar_wrapper[double](obj)
            elif (isinstance(obj, long)):
                var = create_scalar_wrapper[long](obj)
            elif (isinstance(obj, int)):
                var = create_scalar_wrapper[int](obj)
            # TODO handle boolean.
            #elif (isinstance(obj, bool)):
            #    var = create_scalar_wrapper[bool](obj)
            else:
                raise ValueError("invalid scalar type")
            ret.mat = self.mat.apply_scalar(var, op)
        else:
            mat = array(obj)
            ret.mat = self.mat.mapply2(mat.mat, op)
        ret.init_attr()
        return ret

    def sapply(self, op):
        cdef PyMatrix ret = PyMatrix()
        ret.mat = self.mat.sapply(op)
        ret.init_attr()
        return ret

#        matrix_wrapper groupby_row(matrix_wrapper labels, bulk_op_idx_t op) const
#        matrix_wrapper groupby_row(matrix_wrapper labels, bulk_op_idx_t op) const

def array(arr, dtype=None, copy=True, order='K'):
    cdef np.ndarray ndarr
    cdef PyMatrix ret = PyMatrix()

    if (isinstance(arr, np.ndarray)):
        ndarr = arr
    else:
        ndarr = np.array(arr)
    if (not ndarr.flags.contiguous):
        ndarr = np.ascontiguousarray(ndarr)

    if ((order == 'K' or order == 'C') and ndarr.flags.c_contiguous):
        order = 'C'
    elif ((order == 'K' or order == 'C') and ndarr.flags.f_contiguous):
        order = 'F'

    # TODO this is a bit too hacky. Is there a better way?
    cdef intptr_t addr = ctypes.c_void_p(ndarr.ctypes.data).value
    if (ndarr.ndim == 1):
        ret.mat = matrix_wrapper(addr, ndarr.shape[0], ndarr.dtype.char)
    elif (ndarr.ndim == 2):
        ret.mat = matrix_wrapper(addr, ndarr.shape[0], ndarr.shape[1],
                ndarr.dtype.char, order)
    else:
        raise ValueError("don't support more than 2 dimensions")

    ret.init_attr()
    if dtype is None:
        return ret
    else:
        return ret.cast_ele_type(dtype)

def empty_like(a, dtype=None, order='K', subok=True):
    cdef PyMatrix ret = PyMatrix()
    shape = a.shape
    if (dtype is None):
        dtype = a.dtype

    # TODO what is the input array isn't contiguous.
    if (order == 'K' and a.flags.c_contiguous):
        order = 'C'
    elif (order == 'K' and a.flags.f_contiguous):
        order = 'F'
    if (order != 'C' and order != 'F'):
        order = 'C'

    if (len(shape) == 1):
        ret.mat = matrix_wrapper(shape[0], dtype)
    elif (len(shape) == 2):
        ret.mat = matrix_wrapper(shape[0], shape[1], dtype, order)
    else:
        raise ValueError("don't support more than 2 dimensions")
    ret.init_attr()
    return ret

def empty(shape, dtype='f', order='C'):
    cdef PyMatrix ret = PyMatrix()
    if (np.isscalar(shape)):
        ret.mat = matrix_wrapper(shape, dtype)
    elif (len(shape) == 1):
        ret.mat = matrix_wrapper(shape[0], dtype)
    elif (len(shape) == 2):
        ret.mat = matrix_wrapper(shape[0], shape[1], dtype, order)
    else:
        raise ValueError("don't support more than 2 dimensions")
    ret.init_attr()
    return ret

def init_val(PyMatrix data, dtype, val):
    if (dtype == 'f' or dtype == 'd' or dtype == 'g'):
        data.mat.init_const_float(val)
    else:
        data.mat.init_const_int(val)

def ones(shape, dtype='f', order='C'):
    cdef PyMatrix ret = empty(shape, dtype, order)
    init_val(ret, dtype, 1)
    ret.init_attr()
    return ret

def zeros(shape, dtype='f', order='C'):
    cdef PyMatrix ret = empty(shape, dtype, order)
    init_val(ret, dtype, 0)
    ret.init_attr()
    return ret

def arange(start, stop, step=1, dtype=None):
    cdef size_t l = (stop - start) / step
    cdef PyMatrix ret
    if (isinstance(start, float)):
        ret = empty(l, 'd', order='F')
        ret.mat.init_seq[float](start, step, 0)
    elif (isinstance(start, long)):
        ret = empty(l, 'l', order='F')
        ret.mat.init_seq[long](start, step, 0)
    elif (isinstance(start, int)):
        ret = empty(l, 'i', order='F')
        ret.mat.init_seq[int](start, step, 0)
    else:
        raise ValueError("invalid scalar type")
    ret.init_attr()
    return ret

def sum(PyMatrix a, axis=None, dtype=None, out=None, keepdims=False):
    return a.sum(axis, dtype, out, keepdims)

def prod(PyMatrix a, axis=None, dtype=None, out=None, keepdims=False):
    return a.prod(axis, dtype, out, keepdims)

def mean(PyMatrix a, axis=None, dtype=None, out=None, keepdims=False):
    return a.mean(axis, dtype, out, keepdims)

def average(PyMatrix a, axis=None, weights=None, returned=False):
    if weights is not None and axis is None:
        if (a.shape != weights.shape):
            raise ValueError("weights need to have the same shape as a")
        else:
            a = a * weights
            wsum = sum(weights)
    elif weights is not None:
        if (weights.ndim > 1):
            raise ValueError("weights need to be a 1D array")
        elif (axis == 0):
            a = a.mapply_cols(weights, OP_MUL)
        elif (axis == 1):
            a = a.mapply_rows(weights, OP_MUL)
        else:
            raise ValueError("invalid axis")
        wsum = sum(weights)
    elif axis is None:
        wsum = a.mat.get_num_rows() * a.mat.get_num_cols()
    else:
        wsum = a.shape[axis]
    if (returned):
        return (sum(a, axis)/wsum, wsum)
    else:
        return sum(a, axis)/wsum

def dot(PyMatrix a, b, out=None):
    return a.dot(b, out)

def sqrt(PyMatrix x):
    return x.sapply(UOP_SQRT)

def absolute(PyMatrix x):
    return x.sapply(UOP_ABS)

def where(PyMatrix condition, PyMatrix x, PyMatrix y):
    cdef PyMatrix ret = PyMatrix()
    ret.mat = condition.mat.ifelse(x.mat, y.mat)
    ret.init_attr()
    return ret

def init_flashpy(conf_file=""):
    return init_flashpy_c(conf_file)
