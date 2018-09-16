from sklearn.datasets import make_blobs
import numpy as np

nsamples = 3145728
nfeatures = 16
ncenters = 100
data = make_blobs(nsamples, nfeatures, ncenters,
        center_box=(0,20), shuffle=False)[0]

fn = "r"+str(nfeatures)+"_c"+str(nsamples)+"_k"+str(ncenters)
"""
# Text save
# rowwise
np.savetxt(fn + "_rw.txt", data, fmt="%.4f", delimiter=" ")
# colwise
np.savetxt(fn + "_cw.txt", data.T, fmt="%.4f", delimiter=" ")

# Binary save
# rowwise
f = open(fn + "_rw.dat", "wb")
data.tofile(f, format="%.4f")
f.close()
"""

# colwise
f = open(fn + "_cw.dat", "wb")
data.T.tofile(f, format="%.4f")
f.close()