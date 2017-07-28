import numpy as np

from .norm import norm
from .. import dot
from .. import concatenate

def gram_schmidt(mat):
    basis = []
    ncol = mat.shape[1]
    r = np.zeros(shape=[ncol, ncol])
    for i in xrange(0, mat.shape[1]):
        v = mat[:,i]
        bv = None
        if (len(basis) > 0):
            B = concatenate(basis, axis=1)
            bv = np.array(dot(B.T, v))
            w = v - dot(B, bv)
            w.set_cached(True)
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
