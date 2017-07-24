import numpy as np

from .. import asarray
from .. import sqrt
from .svd import svds

def norm(x, ord=None, axis=None, keepdims=False):
    x = asarray(x)
    if (x.ndim > 2):
        raise ValueError("norm doesn't work on high-dim array")
    if (ord is 'fro'):
        if (x.ndim == 1):
            raise ValueError("Frobenius norm doesn't run on a vector")
        ord = None
    if (ord is None):
        ret = sqrt((x * x).sum(axis))
    elif (ord == np.inf):
        if (x.ndim == 2 and axis is None):
            ret = abs(x).sum(axis=1).max()
        elif (x.ndim == 2):
            ret = abs(x).max(axis=axis)
        else:
            ret = abs(x).max()
    elif (ord == -np.inf):
        if (x.ndim == 2 and axis is None):
            ret = abs(x).sum(axis=1).min()
        elif (x.ndim == 2):
            ret = abs(x).min(axis=axis)
        else:
            ret = abs(x).min()
    elif (ord == 0):
        if (x.ndim == 2):
            raise ValueError("0-norm doesn't work on a matrix")
        return (x != 0).sum()
    elif (type(ord) == np.int and x.ndim == 2 and axis is None):
        if (ord == 1):
            ret = abs(x).sum(axis=0).max()
        elif (ord == -1):
            ret = abs(x).sum(axis=0).min()
        elif (ord == 2):
            u, s, v = svds(x, k=1)
            ret = s[0]
        elif (ord == -2):
            # TODO
            ret = None
        else:
            raise ValueError("unknown order of norm")
    elif (type(ord) == np.int and x.ndim == 2):
        ret = (abs(x) ** ord).sum(axis=axis) ** (1./ord)
    elif (type(ord) == np.int):
        ret = (abs(x) ** ord).sum() ** (1./ord)
    else:
        raise ValueError("we don't support other norms yet")
    return ret
