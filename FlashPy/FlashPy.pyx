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

cdef enum bulk_uop_idx_t:
    UOP_NEG, UOP_SQRT, UOP_ABS, UOP_NOT, UOP_SQ,
    UOP_CEIL, UOP_FLOOR, UOP_ROUND,
    UOP_LOG, UOP_LOG2, UOP_LOG10

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

#        void init_seq[T](T start, T stride, bool byrow)
#        void init_const[T](T val)

        size_t get_num_rows() const
        size_t get_num_cols() const
        size_t get_entry_size() const
        string get_type_str() const
        np.NPY_TYPES get_type_py() const
        bool is_in_mem() const
        bool is_virtual() const
        bool is_vector() const
        bool materialize_self() const
        matrix_wrapper get_cols(const vector[long] &idxs) const
        matrix_wrapper get_rows(const vector[long] &idxs) const
        matrix_wrapper get_cols(matrix_wrapper idxs) const
        matrix_wrapper get_rows(matrix_wrapper idxs) const
        matrix_wrapper get_cols(size_t start, size_t end) const
        matrix_wrapper get_rows(size_t start, size_t end) const
        matrix_wrapper set_cols(const vector[long] &idxs, matrix_wrapper cols)
        matrix_wrapper set_rows(const vector[long] &idxs, matrix_wrapper rows)
        const char *get_raw_arr() const
        matrix_wrapper transpose() const
        matrix_wrapper conv_store(bool in_mem, int num_nodes) const
        matrix_wrapper inner_prod(matrix_wrapper m, bulk_op_idx_t left_op,
                bulk_op_idx_t right_op) const
        matrix_wrapper multiply(matrix_wrapper m) const
        matrix_wrapper agg_row(bulk_op_idx_t op) const
        matrix_wrapper agg_col(bulk_op_idx_t op) const
        matrix_wrapper groupby_row(matrix_wrapper labels, bulk_op_idx_t op) const
        matrix_wrapper groupby_row(matrix_wrapper labels, bulk_op_idx_t op) const
        matrix_wrapper mapply_cols(matrix_wrapper vals, bulk_op_idx_t op) const
        matrix_wrapper mapply_rows(matrix_wrapper vals, bulk_op_idx_t op) const
        matrix_wrapper mapply2(matrix_wrapper m, bulk_op_idx_t op) const
        matrix_wrapper sapply(bulk_uop_idx_t op) const

cdef class PyMatrix:
    cdef matrix_wrapper mat      # hold a C++ instance which we're wrapping
    cdef readonly int ndim
    cdef readonly object shape

    def __cinit__(self):
        self.mat = matrix_wrapper()

    def __init__(self):
        self.ndim = 0
        self.shape = (0, 0)

    def __array__(self):
        cdef char *src = self.mat.get_raw_arr()
        if (src == NULL):
            return None
        cdef np.npy_intp shape[2]
        shape[0] = self.shape[0]
        shape[1] = self.shape[1]
        return np.PyArray_SimpleNewFromData(self.ndim, shape,
                self.mat.get_type_py(), src)

    # Special Methods Table
    # http://cython.readthedocs.io/en/latest/src/reference/special_methods_table.html

    def __richcmp__(PyMatrix x, PyMatrix y, int op):
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
            ret.mat = x.mat.mapply2(y.mat, OP_LT)
        elif (op == 2):
            ret.mat = x.mat.mapply2(y.mat, OP_EQ)
        elif (op == 4):
            ret.mat = x.mat.mapply2(y.mat, OP_GT)
        elif (op == 1):
            ret.mat = x.mat.mapply2(y.mat, OP_LE)
        elif (op == 3):
            ret.mat = x.mat.mapply2(y.mat, OP_NEQ)
        elif (op == 5):
            ret.mat = x.mat.mapply2(y.mat, OP_GE)
        else:
            print("invalid argument")
        return ret

    def __add__(PyMatrix x, PyMatrix y):
        cdef PyMatrix ret = PyMatrix()
        ret.mat = x.mat.mapply2(y.mat, OP_ADD)
        ret.init_shape()
        return ret

    def __sub__(PyMatrix x, PyMatrix y):
        cdef PyMatrix ret = PyMatrix()
        ret.mat = x.mat.mapply2(y.mat, OP_SUB)
        ret.init_shape()
        return ret

    def __mul__(PyMatrix x, PyMatrix y):
        cdef PyMatrix ret = PyMatrix()
        ret.mat = x.mat.mapply2(y.mat, OP_MUL)
        ret.init_shape()
        return ret

    def __div__(PyMatrix x, PyMatrix y):
        cdef PyMatrix ret = PyMatrix()
        ret.mat = x.mat.mapply2(y.mat, OP_DIV)
        ret.init_shape()
        return ret

    def __floordiv__(PyMatrix x, PyMatrix y):
        cdef PyMatrix ret = PyMatrix()
        ret.mat = x.mat.mapply2(y.mat, OP_IDIV)
        ret.init_shape()
        return ret

    def __mod__(PyMatrix x, PyMatrix y):
        cdef PyMatrix ret = PyMatrix()
        ret.mat = x.mat.mapply2(y.mat, OP_MOD)
        ret.init_shape()
        return ret

    def __and__(PyMatrix x, PyMatrix y):
        cdef PyMatrix ret = PyMatrix()
        ret.mat = x.mat.mapply2(y.mat, OP_AND)
        ret.init_shape()
        return ret

    def __or__(PyMatrix x, PyMatrix y):
        cdef PyMatrix ret = PyMatrix()
        ret.mat = x.mat.mapply2(y.mat, OP_OR)
        ret.init_shape()
        return ret

    def __neg__(self):
        cdef PyMatrix ret = PyMatrix()
        ret.mat = self.mat.sapply(UOP_NEG)
        ret.init_shape()
        return ret

    def __abs__(self):
        cdef PyMatrix ret = PyMatrix()
        ret.mat = self.mat.sapply(UOP_ABS)
        ret.init_shape()
        return ret

    def __len__(self):
        return self.mat.get_num_rows()

    def init_shape(self):
        self.shape = (self.mat.get_num_rows(), self.mat.get_num_cols())
        if (self.mat.is_vector()):
            self.ndim = 1
        else:
            self.ndim = 2

    def is_in_mem(self):
        return self.mat.is_in_mem()

    def is_virtual(self):
        return self.mat.is_virtual()

    def materialize_self(self):
        return self.mat.materialize_self()

    def get_cols(self, array.array idxs):
        cdef vector[long] cidxs
        cdef long *p = idxs.data.as_longs
        cidxs.assign(p, p + len(idxs))

        cdef PyMatrix ret = PyMatrix()
        ret.mat = self.mat.get_cols(cidxs)
        ret.init_shape()
        return ret

    def transpose(self):
        cdef PyMatrix ret = PyMatrix()
        ret.mat = self.mat.transpose()
        ret.init_shape()
        return ret

    def conv_store(self, bool in_mem, int num_nodes):
        cdef PyMatrix ret = PyMatrix()
        ret.mat = self.mat.conv_store(in_mem, num_nodes)
        ret.init_shape()
        return ret
"""
# TODO this function should have the same interface as numpy.array.
def array(np.ndarray arr, string dtype):
    cdef PyMatrix ret = PyMatrix()
    # TODO this is a bit too hacky. Is there a better way?
    cdef intptr_t addr = ctypes.c_void_p(arr.ctypes.data).value
    if (arr.ndim == 1):
        ret.mat = matrix_wrapper(addr, arr.shape[0], dtype)
    elif (arr.ndim == 2):
        ret.mat = matrix_wrapper(addr, arr.shape[0], arr.shape[1], dtype, "c")
    else:
        raise ValueError("don't support more than 2 dimensions")
    ret.init_shape()
    return ret

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

    if (len(shape) == 1):
        ret.mat = matrix_wrapper(shape[0], dtype)
    elif (len(shape) == 2):
        ret.mat = matrix_wrapper(shape[0], shape[1], dtype, order)
    else:
        raise ValueError("don't support more than 2 dimensions")
    ret.init_shape()
    return ret

def empty(shape, dtype='f', order='C'):
    cdef PyMatrix ret = PyMatrix()
    if (len(shape) == 1):
        ret.mat = matrix_wrapper(shape[0], dtype)
    elif (len(shape) == 2):
        ret.mat = matrix_wrapper(shape[0], shape[1], dtype, order)
    else:
        raise ValueError("don't support more than 2 dimensions")
    ret.init_shape()
    return ret

def init_val(PyMatrix data, dtype, val):
    if (dtype == 'f'):
        data.mat.init_const[float](val)

def ones(shape, dtype='f', order='C'):
    cdef PyMatrix ret = PyMatrix()
    if (len(shape) == 1):
        ret.mat = matrix_wrapper(1, shape[0], dtype)
    elif (len(shape) == 2):
        ret.mat = matrix_wrapper(1, shape[0], shape[1], dtype, order)
    else:
        raise ValueError("don't support more than 2 dimensions")
    init_val(ret, dtype, 1)
    ret.init_shape()
    return ret

def zeros(shape, dtype='f', order='C'):
    cdef PyMatrix ret = PyMatrix()
    if (len(shape) == 1):
        ret.mat = matrix_wrapper(0, shape[0], dtype)
    elif (len(shape) == 2):
        ret.mat = matrix_wrapper(0, shape[0], shape[1], dtype, order)
    else:
        raise ValueError("don't support more than 2 dimensions")
    ret.init_shape()
    return ret
"""

#        matrix_wrapper inner_prod(matrix_wrapper m, bulk_op_idx_t left_op,
#                bulk_op_idx_t right_op) const
#        matrix_wrapper agg_row(bulk_op_idx_t op) const
#        matrix_wrapper agg_col(bulk_op_idx_t op) const
#        matrix_wrapper groupby_row(matrix_wrapper labels, bulk_op_idx_t op) const
#        matrix_wrapper groupby_row(matrix_wrapper labels, bulk_op_idx_t op) const
#        matrix_wrapper mapply_cols(matrix_wrapper vals, bulk_op_idx_t op) const
#        matrix_wrapper mapply_rows(matrix_wrapper vals, bulk_op_idx_t op) const
#        matrix_wrapper mapply2(matrix_wrapper m, bulk_op_idx_t op) const
#        matrix_wrapper sapply(bulk_uop_idx_t op) const
