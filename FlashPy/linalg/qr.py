import numpy as np

from .norm import norm
from .. import dot
from .. import concatenate
from .. import atleast_2d

def gram_schmidt(mat):
    mat = mat.copy(order='F')
    mat.materialize_self()
    basis = []
    ncol = mat.shape[1]
    r = np.zeros(shape=[ncol, ncol])
    for i in xrange(0, mat.shape[1]):
        v = mat[:,i]
        bv = None
        if (len(basis) > 0):
            B = concatenate(basis, axis=1)
            B = atleast_2d(B)
            bv = np.array(dot(B.T, v))
            w = v - dot(B, bv)
            w.materialize_self()
            # TODO somehow caching the data has much worse performance.
            # I need to fix it.
            #w.set_cached(True)
        else:
            w = v
        if (np.array(abs(w).max())[0] > 1e-10):
            e = w/norm(w)
            basis.append(e)
            n = len(basis) - 1
            r[n,n] = np.array(dot(e.T, v))[0]
            if (n > 0):
                r[0:n,n] = bv.ravel()
    return (concatenate(basis, axis=1), r)

def qr(a, mode='reduced'):
    return gram_schmidt(a)
