import numpy as np
import FlashPy as fp
from scipy import linalg
from scipy.sparse.linalg import LinearOperator

#from ..sparse import issparse

def issparse(a):
    return False

def svd(a, compute_uv=True, nu=None, nv=None, tol=1e-8):
    if nu is None:
        nu = min(a.shape)
    if nv is None:
        nv = min(a.shape)

    if a.ndim != 2:
        raise ValueError("expect a matrix")
    n, p = a.shape
    comp_right = False
    if n > p:
        comp_right = True

    nev = max(nu, nv)
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

    if (x_prod is not None and (size < 100 or nev >= size / 2)):
        x_prod = np.array(x_prod)
        vals, vecs = linalg.eigh(x_prod)
        vals = vals[::-1][0:nev]
        vecs = vecs[:,::-1][:,0:nev]
    else:
        vals, vecs = linalg.eigsh(multiply, k=nev, tol=tol)

    def rescale(x):
        x.set_cached(True)
        scal = fp.sqrt(fp.sum(x * x, axis=0))
        return x.mapply_rows(scal, fp.bop_div)

    if (comp_right):
        v = None
        if (nv > 0):
            if (nv < nev):
                v = fp.array(vecs[:,0:nv])
            else:
                v = fp.array(vecs)
        u = None
        if (nu > 0):
            u = rescale(a.dot(vecs[:,0:nu]))
    else:
        u = None
        if (nu > 0):
            if (nu < nev):
                u = fp.array(vecs[:,0:nu])
            else:
                u = fp.array(vecs)
        v = None
        if (nv > 0):
            v = rescale(a.T.dot(vecs[:,0:nv]))
    s = np.where(vals > 0, np.sqrt(vals), 0)
    return u, s, v.T
