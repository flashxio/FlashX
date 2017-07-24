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
        if (len(basis) > 0):
            B = concatenate(basis, axis=1)
            bv = np.array(dot(B.T, v))
            r[0:i,i] = bv.ravel()
            w = v - dot(B, bv)
            w.set_cached(True)
        else:
            w = v
        if (w > 1e-10).any():  
            e = w/norm(w)
            basis.append(e)
            r[i,i] = np.array(dot(e.T, v))[0]
    return (concatenate(basis, axis=1), r)

def qr(a, mode='reduced'):
    return gram_schmidt(a)
