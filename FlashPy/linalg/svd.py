import numpy as np
import FlashPy as fp
from scipy import linalg
from scipy.sparse.linalg import eigsh
from scipy.sparse.linalg import LinearOperator

#from ..sparse import issparse

def issparse(a):
    return False

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
        x_prod = a.T.dot(a)
        multiply = LinearOperator(matvec=lambda v: x_prod.dot(v),
                shape=(p, p))
    else:
        size = n
        x_prod = a.dot(a.T)
        multiply = LinearOperator(matvec=lambda v: x_prod.dot(v),
                shape=(n, n))

    if (x_prod is not None and (size < 100 or k >= size / 2)):
        x_prod = np.array(x_prod)
        vals, vecs = linalg.eigh(x_prod)
        vals = vals[::-1][0:k]
        vecs = vecs[:,::-1][:,0:k]
    else:
        vals, vecs = eigsh(multiply, k=k, tol=tol)

    def rescale(x):
        x.set_cached(True)
        scal = fp.sqrt(fp.sum(x * x, axis=0))
        return x.mapply_rows(scal, fp.bop_div)

    if (comp_right):
        v = fp.array(vecs)
        u = rescale(a.dot(vecs))
    else:
        u = fp.array(vecs)
        v = rescale(a.T.dot(vecs))
    s = np.where(vals > 0, np.sqrt(vals), 0)
    return u, s, v.T
