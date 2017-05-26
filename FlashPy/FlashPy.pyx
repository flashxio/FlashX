import numpy as np
# "cimport" is used to import special compile-time information
# about the numpy module (this is stored in a file numpy.pxd which is
# currently part of the Cython distribution).
cimport numpy as np
from libc.stdlib cimport free, malloc
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
        matrix_wrapper(size_t nrow, size_t ncol, string t, string layout)
        void init_seq[T](T start, T stride, size_t nrow, size_t ncol,
                string layout, bool byrow, int num_nodes, bool in_mem)
        size_t get_num_rows() const
        size_t get_num_cols() const
        size_t get_entry_size() const
        string get_type_str() const
        np.NPY_TYPES get_type_py() const
        bool is_in_mem() const
        bool is_virtual() const
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
    def __cinit__(self):
        self.mat = matrix_wrapper()

    def __array__(self):
        cdef char *src = self.mat.get_raw_arr()
        if (src == NULL):
            return None

        cdef np.npy_intp shape[2]
        shape[0] = <np.npy_intp> self.mat.get_num_rows()
        shape[1] = <np.npy_intp> self.mat.get_num_cols()
        return np.PyArray_SimpleNewFromData(2, shape, self.mat.get_type_py(), src)

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
        return ret

    def __sub__(PyMatrix x, PyMatrix y):
        cdef PyMatrix ret = PyMatrix()
        ret.mat = x.mat.mapply2(y.mat, OP_SUB)
        return ret

    def __mul__(PyMatrix x, PyMatrix y):
        cdef PyMatrix ret = PyMatrix()
        ret.mat = x.mat.mapply2(y.mat, OP_MUL)
        return ret

    def __div__(PyMatrix x, PyMatrix y):
        cdef PyMatrix ret = PyMatrix()
        ret.mat = x.mat.mapply2(y.mat, OP_DIV)
        return ret

    def __floordiv__(PyMatrix x, PyMatrix y):
        cdef PyMatrix ret = PyMatrix()
        ret.mat = x.mat.mapply2(y.mat, OP_IDIV)
        return ret

    def __mod__(PyMatrix x, PyMatrix y):
        cdef PyMatrix ret = PyMatrix()
        ret.mat = x.mat.mapply2(y.mat, OP_MOD)
        return ret

    def __and__(PyMatrix x, PyMatrix y):
        cdef PyMatrix ret = PyMatrix()
        ret.mat = x.mat.mapply2(y.mat, OP_AND)
        return ret

    def __or__(PyMatrix x, PyMatrix y):
        cdef PyMatrix ret = PyMatrix()
        ret.mat = x.mat.mapply2(y.mat, OP_OR)
        return ret

    def __neg__(self):
        cdef PyMatrix ret = PyMatrix()
        ret.mat = self.mat.sapply(UOP_NEG)
        return ret

    def __abs__(self):
        cdef PyMatrix ret = PyMatrix()
        ret.mat = self.mat.sapply(UOP_ABS)
        return ret

    def __len__(self):
        return self.mat.get_num_rows() * self.mat.get_num_cols()

    @staticmethod
    def create_seq(long start, long stride, unsigned long nrow, unsigned long ncol,
            string layout, bool byrow, int num_node, bool in_mem):
        cdef PyMatrix mat = PyMatrix()
        mat.mat.init_seq[long](start, stride, nrow, ncol, layout, byrow,
                num_node, in_mem)
        return mat

    def get_num_rows(self):
        return self.mat.get_num_rows()

    def get_num_cols(self):
        return self.mat.get_num_cols()

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
        return ret

    def transpose(self):
        cdef PyMatrix ret = PyMatrix()
        ret.mat = self.mat.transpose()
        return ret

    def conv_store(self, bool in_mem, int num_nodes):
        cdef PyMatrix ret = PyMatrix()
        ret.mat = self.mat.conv_store(in_mem, num_nodes)
        return ret

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
