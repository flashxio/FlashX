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

ctypedef int bulk_op_idx_t
ctypedef int bulk_uop_idx_t

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
