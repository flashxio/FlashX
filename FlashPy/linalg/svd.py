from numpy import array as np_array
from numpy import where as np_where
from numpy import sqrt as np_sqrt
from numpy import sum as np_sum
from scipy import linalg
from scipy.sparse.linalg import eigsh
from scipy.sparse.linalg import LinearOperator

from ..sparse import issparse
from .. import sqrt as fp_sqrt
from .. import array as fp_array
from ..mat import bop_div as fp_bop_div

def svds(a, k=6, tol=0):
    if a.ndim != 2:
        raise ValueError("expect a matrix")
    n, p = a.shape
    comp_right = False
    if n > p:
        comp_right = True

    x_prod = None
    if (issparse(a) and comp_right):
        size = p
        multiply = LinearOperator(matvec=lambda v: a.T.dot(a.dot(v)),
                shape=(p, p))
    elif (issparse(a)):
        size = n
        multiply = LinearOperator(matvec=lambda v: a.dot(a.T.dot(v)),
                shape=(n, n))
    elif (comp_right):
        size = p
        x_prod = np_array(a.T.dot(a))
        multiply = LinearOperator(matvec=lambda v: x_prod.dot(v),
                shape=(p, p))
    else:
        size = n
        x_prod = np_array(a.dot(a.T))
        multiply = LinearOperator(matvec=lambda v: x_prod.dot(v),
                shape=(n, n))

    if (x_prod is not None and (size < 100 or k >= size / 2)):
        x_prod = np_array(x_prod)
        vals, vecs = linalg.eigh(x_prod)
        vals = vals[::-1][0:k]
        vecs = vecs[:,::-1][:,0:k]
    else:
        vals, vecs = eigsh(multiply, k=k, tol=tol)

    def rescale(x):
        x.set_cached(True)
        scal = fp_sqrt(np_sum(x * x, axis=0))
        return x.mapply_rows(scal, fp_bop_div)

    if (comp_right):
        v = fp_array(vecs)
        u = rescale(a.dot(vecs))
    else:
        u = fp_array(vecs)
        v = rescale(a.T.dot(vecs))
    s = np_where(vals > 0, np_sqrt(vals), 0)
    return u, s, v.T
