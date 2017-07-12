import flashpy as fp
import numpy as np
from scipy import linalg
import scipy.sparse.linalg as sp_linalg
import sparse as fp_sparse

fp.init_flashpy()

def verify(fp_arr, np_arr):
    assert fp_arr.ndim == np_arr.ndim
    assert fp_arr.shape == np_arr.shape
    assert fp_arr.size == np_arr.size
    assert fp_arr.itemsize == np_arr.itemsize
    assert fp_arr.nbytes == np_arr.nbytes
    tmp = np.array(fp_arr, copy=True)
    assert np.all(np.absolute(tmp - np_arr) < 1e-12)

def verify_init(dtype):
    print("verify " + dtype)
    fp_mat1 = fp.empty([25, 10], dtype=dtype)
    assert fp_mat1.shape[0] == 25
    assert fp_mat1.shape[1] == 10
    assert fp_mat1.dtype == dtype

    fp_mat2 = fp.empty_like(fp_mat1)
    assert fp_mat1.shape == fp_mat2.shape
    assert fp_mat1.dtype == fp_mat2.dtype

    fp_mat1 = fp.ones([25, 10], dtype=dtype)
    assert fp_mat1.shape[0] == 25
    assert fp_mat1.shape[1] == 10
    assert fp_mat1.dtype == dtype

    fp_mat1 = fp.zeros([25, 10], dtype=dtype)
    assert fp_mat1.shape[0] == 25
    assert fp_mat1.shape[1] == 10
    assert fp_mat1.dtype == dtype

verify_init("f")
verify_init("d")
verify_init("g")
verify_init("b")
verify_init("h")
verify_init("i")
verify_init("l")
verify_init("H")
verify_init("I")
verify_init("L")

print("test transpose a vector")
np_mat1 = np.arange(100)
fp_mat1 = fp.array(np_mat1)
verify(fp_mat1, np_mat1)
verify(fp_mat1.T, np_mat1.T)
verify(fp_mat1.transpose(), np_mat1.transpose())
print("test transpose a one-col matrix")
np_mat1 = np.random.normal(scale=100, size=[10, 1])
fp_mat1 = fp.array(np_mat1)
verify(fp_mat1, np_mat1)
verify(fp_mat1.T, np_mat1.T)
verify(fp_mat1.transpose(), np_mat1.transpose())
print("test transpose a matrix")
np_mat1 = np.random.normal(scale=100, size=[10, 5])
fp_mat1 = fp.array(np_mat1)
verify(fp_mat1, np_mat1)
verify(fp_mat1.T, np_mat1.T)
verify(fp_mat1.transpose(), np_mat1.transpose())

print("test C-contiguous array")
np_mat1 = np.random.normal(scale=100, size=[10, 10])
fp_mat1 = fp.array(np_mat1)
verify(fp_mat1, np_mat1)
np_mat1 = np.array(fp_mat1, copy=True)
verify(fp_mat1, np_mat1)

print("test F-contiguous array")
np_mat1 = np.array(np_mat1, order='F', copy=True)
fp_mat1 = fp.array(np_mat1)
verify(fp_mat1, np_mat1)
np_mat1 = np.array(fp_mat1, copy=True)
verify(fp_mat1, np_mat1)

print("test indexing (get)")
np_mat1 = np.random.normal(scale=100, size=[25, 10])
fp_mat1 = fp.array(np_mat1)
np_mat2 = np.random.normal(scale=100, size=[25, 10])
fp_mat2 = fp.array(np_mat2)
np_vec2 = np.random.normal(scale=100, size=25)
fp_vec2 = fp.array(np_vec2)
np_vec3 = np.random.normal(scale=100, size=10)
fp_vec3 = fp.array(np_vec3)

verify(fp_mat1[1], np_mat1[1])
verify(fp_mat1[1:3], np_mat1[1:3])
verify(fp_mat1[1:7:2], np_mat1[1:7:2])
verify(fp_mat1[1:3, 1:3], np_mat1[1:3, 1:3])
verify(fp_mat1[1:7:2, 1:7:2], np_mat1[1:7:2, 1:7:2])
idx = [1, 3, 5]
verify(fp_mat1[idx], np_mat1[idx])
#verify(fp_mat1[idx, idx], np_mat1[idx, idx])

print("test indexing (set)")
np_mat1[:] = np_mat2
fp_mat1[:] = fp_mat2
verify(fp_mat1, np_mat1)

print("test1")
tmp = np.random.normal(scale=100, size=25)
np_mat1[:,1] = tmp
fp_mat1[:,1] = fp.array(tmp)
verify(fp_mat1, np_mat1)

print("test2")
tmp = np.random.normal(scale=100, size=[25, 4])
np_mat1[:,1:9:2] = tmp
fp_mat1[:,1:9:2] = fp.array(tmp)
verify(fp_mat1, np_mat1)

np_mat1 = np_mat1.T
fp_mat1 = fp_mat1.T

print("test3")
tmp = np.array(np.random.normal(scale=100, size=[1, 25]))
np_mat1[1,:] = tmp
tmp2 = fp_mat1
fp_mat1[1,:] = fp.array(tmp)
verify(fp_mat1, np_mat1)

print("test4")
tmp = np.random.normal(scale=100, size=[4, 25])
np_mat1[1:9:2,:] = tmp
fp_mat1[1:9:2,:] = fp.array(tmp)
verify(fp_mat1, np_mat1)

np_mat1 = np_mat1.T
fp_mat1 = fp_mat1.T

np_mat1 = np.random.normal(scale=100, size=[25, 10])
fp_mat1 = fp.array(np_mat1)
print("test +")
fp_res = fp_mat1 + fp_mat2
np_res = np_mat1 + np_mat2
verify(fp_res, np_res)
fp_res = fp_mat1 + 1
np_res = np_mat1 + 1
verify(fp_res, np_res)
fp_res = 1 + fp_mat1
np_res = 1 + np_mat1
verify(fp_res, np_res)
fp_res = fp_mat1 + fp_vec3
np_res = np_mat1 + np_vec3
verify(fp_res, np_res)
fp_res = fp_vec3 + fp_mat1
np_res = np_vec3 + np_mat1
verify(fp_res, np_res)
fp_res = fp_mat1 + fp_vec3[np.newaxis,:]
np_res = np_mat1 + np_vec3[np.newaxis,:]
verify(fp_res, np_res)
fp_res = fp_mat1 + fp_vec2[:,np.newaxis]
np_res = np_mat1 + np_vec2[:,np.newaxis]
verify(fp_res, np_res)
fp_res = fp_vec3[np.newaxis,:] + fp_mat1
np_res = np_vec3[np.newaxis,:] + np_mat1
verify(fp_res, np_res)
fp_res = fp_vec2[:,np.newaxis] + fp_mat1
np_res = np_vec2[:,np.newaxis] + np_mat1
verify(fp_res, np_res)

fp_tmp = fp_mat1
np_tmp = np_mat1.copy()
fp_tmp += fp_mat2
np_tmp += np_mat2
verify(fp_tmp, np_tmp)

print("test -")
fp_res = fp_mat1 - fp_mat2
np_res = np_mat1 - np_mat2
verify(fp_res, np_res)

fp_res = fp_mat1 - 1
np_res = np_mat1 - 1
verify(fp_res, np_res)

fp_res = 1 - fp_mat1
np_res = 1 - np_mat1
verify(fp_res, np_res)

fp_tmp = fp_mat1
np_tmp = np_mat1.copy()
fp_tmp -= fp_mat2
np_tmp -= np_mat2
verify(fp_tmp, np_tmp)

print("test *")
fp_res = fp.array(np_mat1, "l") * fp.array(np_mat2, "l")
np_res = np.array(np_mat1, "l") * np.array(np_mat2, "l")
verify(fp_res, np_res)

fp_res = fp.array(np_mat1, "l") * 2
np_res = np.array(np_mat1, "l") * 2
verify(fp_res, np_res)

fp_res = 2 * fp.array(np_mat1, "l")
np_res = 2 * np.array(np_mat1, "l")
verify(fp_res, np_res)

fp_tmp = fp.array(np_mat1, "l")
np_tmp = np.array(np_mat1, "l", copy=True)
fp_tmp *= fp.array(fp_mat2, "l")
np_tmp *= np.array(np_mat2, "l")
verify(fp_tmp, np_tmp)

print("test /")
fp_res = fp_mat1 / fp_mat2
np_res = np_mat1 / np_mat2
verify(fp_res, np_res)

fp_res = fp_mat1 / 2
np_res = np_mat1 / 2
verify(fp_res, np_res)

fp_res = 2 / fp_mat1
np_res = 2 / np_mat1
verify(fp_res, np_res)

fp_tmp = fp_mat1
np_tmp = np_mat1.copy()
fp_tmp /= fp_mat2
np_tmp /= np_mat2
verify(fp_tmp, np_tmp)

print("test //")
fp_res = fp_mat1 // fp_mat2
np_res = np_mat1 // np_mat2
verify(fp_res, np_res)

fp_res = fp_mat1 // 2
np_res = np_mat1 // 2
verify(fp_res, np_res)

fp_res = 100 // fp_mat1
np_res = 100 // np_mat1
verify(fp_res, np_res)

fp_tmp = fp_mat1
np_tmp = np_mat1.copy()
fp_tmp //= fp_mat2
np_tmp //= np_mat2
verify(fp_tmp, np_tmp)

#print("test %")
#fp_res = fp_mat1 % fp_mat2
#np_res = np_mat1 % np_mat2
#verify(fp_res, np_res)
#
#fp_res = fp_mat1 % 2
#np_res = np_mat1 % 2
#verify(fp_res, np_res)
#
#fp_res = 100 % fp_mat1
#np_res = 100 % np_mat1
#verify(fp_res, np_res)
#
#fp_tmp = fp_mat1
#np_tmp = np_mat1.copy()
#fp_tmp %= fp_mat2
#np_tmp %= np_mat2
#verify(fp_tmp, np_tmp)

print("test >")
fp_res = fp_mat1 > fp_mat2
np_res = np_mat1 > np_mat2
verify(fp_res, np_res)

fp_res = fp_mat1 > 2
np_res = np_mat1 > 2
verify(fp_res, np_res)

fp_res = 100 > fp_mat1
np_res = 100 > np_mat1
verify(fp_res, np_res)

print("test <")
fp_res = fp_mat1 < fp_mat2
np_res = np_mat1 < np_mat2
verify(fp_res, np_res)

fp_res = fp_mat1 < 2
np_res = np_mat1 < 2
verify(fp_res, np_res)

fp_res = 100 < fp_mat1
np_res = 100 < np_mat1
verify(fp_res, np_res)

print("test >=")
fp_res = fp_mat1 >= fp_mat2
np_res = np_mat1 >= np_mat2
verify(fp_res, np_res)

fp_res = fp_mat1 >= 2
np_res = np_mat1 >= 2
verify(fp_res, np_res)

fp_res = 100 >= fp_mat1
np_res = 100 >= np_mat1
verify(fp_res, np_res)

print("test <=")
fp_res = fp_mat1 <= fp_mat2
np_res = np_mat1 <= np_mat2
verify(fp_res, np_res)

fp_res = fp_mat1 <= 2
np_res = np_mat1 <= 2
verify(fp_res, np_res)

fp_res = 100 <= fp_mat1
np_res = 100 <= np_mat1
verify(fp_res, np_res)

print("test ==")
fp_res = fp_mat1 == fp_mat2
np_res = np_mat1 == np_mat2
verify(fp_res, np_res)

fp_res = fp_mat1 == 2
np_res = np_mat1 == 2
verify(fp_res, np_res)

fp_res = 100 == fp_mat1
np_res = 100 == np_mat1
verify(fp_res, np_res)

print("test ==")
fp_res = fp_mat1 != fp_mat2
np_res = np_mat1 != np_mat2
verify(fp_res, np_res)

fp_res = fp_mat1 != 2
np_res = np_mat1 != 2
verify(fp_res, np_res)

fp_res = 100 != fp_mat1
np_res = 100 != np_mat1
verify(fp_res, np_res)

print("test abs")
fp_res = abs(fp_mat1)
np_res = abs(np_mat1)
verify(fp_res, np_res)

print("test neg")
fp_res = -fp_mat1
np_res = -np_mat1
verify(fp_res, np_res)

def verify_cast(dtype1, dtype2):
    print("cast " + dtype1 + " to " + dtype2)
    fp_mat1 = fp.empty([25, 10], dtype=dtype1)
    assert fp_mat1.dtype == dtype1
    fp_mat2 = fp_mat1.cast_ele_type(dtype2)
    assert fp_mat2.dtype == dtype2

verify_cast("f", "d")
verify_cast("d", "f")
verify_cast("b", "h")
verify_cast("i", "l")

print("test sum")
np_mat1 = np.random.normal(scale=100, size=[25, 10])
fp_mat1 = fp.array(np_mat1, dtype="l")
np_mat1 = np.array(np_mat1, dtype="l")
fp_res = fp.sum(fp_mat1)
np_res = np.sum(np_mat1)
tmp = np.array(fp_res, copy=True)
assert abs(tmp - np_res) < 1e-13

fp_res = fp.sum(fp_mat1, axis=0)
np_res = np.sum(np_mat1, axis=0)
verify(fp_res, np_res)

fp_res = fp.sum(fp_mat1, axis=0, dtype="d")
np_res = np.sum(np_mat1, axis=0, dtype="d")
verify(fp_res, np_res)

fp_res = fp.sum(fp_mat1, axis=1)
np_res = np.sum(np_mat1, axis=1)
verify(fp_res, np_res)

print("test mean")
fp_res = fp.mean(fp_mat1)
np_res = np.mean(np_mat1)
tmp = np.array(fp_res, copy=True)
assert abs(tmp[0] - np_res) < 1e-14

fp_res = fp.mean(fp_mat1, axis=0)
np_res = np.mean(np_mat1, axis=0)
verify(fp_res, np_res)

fp_res = fp.mean(fp_mat1, axis=1)
np_res = np.mean(np_mat1, axis=1)
verify(fp_res, np_res)

print("test min/max")
fp_res = fp_mat1.min()
np_res = np_mat1.min()
tmp = np.array(fp_res, copy=True)
assert tmp[0] == np_res

fp_res = fp_mat1.min(axis=1)
np_res = np_mat1.min(axis=1)
verify(fp_res, np_res)

fp_res = fp_mat1.min(axis=0)
np_res = np_mat1.min(axis=0)
verify(fp_res, np_res)

print("test argmin/argmax")
fp_res = fp_mat1.argmin(axis=1)
np_res = np_mat1.argmin(axis=1)
verify(fp_res, np_res)

fp_res = fp_mat1.argmax(axis=1)
np_res = np_mat1.argmax(axis=1)
verify(fp_res, np_res)

print("test minimum/maximum")
fp_res = fp.minimum(fp_mat1, fp_mat2)
np_res = np.minimum(np_mat1, np_mat2)
verify(fp_res, np_res)

fp_res = fp.maximum(fp_mat1, fp_mat2)
np_res = np.maximum(np_mat1, np_mat2)
verify(fp_res, np_res)

print("test dot on matrix")
np_mat1 = np.random.normal(scale=100, size=[25, 10])
fp_mat1 = fp.array(np_mat1, dtype="l")
np_mat1 = np.array(np_mat1, dtype="l")
fp_res = fp.dot(fp_mat1, fp_mat1.transpose())
np_res = np.dot(np_mat1, np_mat1.transpose())
verify(fp_res, np_res)

fp_res = fp.dot(fp_mat1, np_mat1.transpose())
verify(fp_res, np_res)

fp_res = fp.dot(np_mat1, fp_mat1.transpose())
verify(fp_res, np_res)

fp_res = fp.dot(np_mat1, np_mat1.transpose())
verify(fp_res, np_res)


print("test dot on vector")
np_mat1 = np.random.normal(scale=100, size=25)
fp_mat1 = fp.array(np_mat1, dtype="l")
np_mat1 = np.array(np_mat1, dtype="l")
fp_res = fp.dot(fp_mat1, fp_mat1)
np_res = np.dot(np_mat1, np_mat1)
assert np.array(fp_res, copy=True) == np_res

print("test casting float to int")
np_mat1 = np.random.normal(scale=100, size=[25, 10])
fp_mat1 = fp.array(np_mat1)
fp_res = np_mat1.astype('i')
np_res = fp_mat1.astype('i')
verify(fp_res, np_res)
np_res = np_mat1.astype('i', order='F')
fp_res = fp_mat1.astype('i', order='F')
verify(fp_res, np_res)
fp_res = np_mat1.astype('i', order='C')
np_res = fp_mat1.astype('i', order='C')
verify(fp_res, np_res)

print("test casting int to float")
np_mat1 = np.arange(1, 100)
fp_mat1 = fp.arange(1, 100)
fp_res = np_mat1.astype('d')
np_res = fp_mat1.astype('d')
verify(fp_res, np_res)

print("test svd")
np_mat1 = np.random.normal(scale=100, size=[1000, 10])
fp_mat1 = fp.array(np_mat1)
np_u, np_s, np_v = linalg.svd(np_mat1, full_matrices=False)
fp_u, fp_s, fp_v = fp.linalg.svds(fp_mat1, k=min(fp_mat1.shape))
verify(fp.absolute(fp_u), np.absolute(np_u))
verify(fp.absolute(fp_v), np.absolute(np_v))
assert all(abs(fp_s - np_s) / np_s < 1e-14)

np_mat1 = np.random.normal(scale=100, size=[10, 1000])
fp_mat1 = fp.array(np_mat1)
np_u, np_s, np_v = linalg.svd(np_mat1, full_matrices=False)
fp_u, fp_s, fp_v = fp.linalg.svds(fp_mat1, k=min(fp_mat1.shape))
verify(fp.absolute(fp_v), np.absolute(np_v))
verify(fp.absolute(fp_v), np.absolute(np_v))
assert all(abs(fp_s - np_s) / np_s < 1e-14)

np_mat1 = np.random.normal(scale=100, size=[1000, 200])
fp_mat1 = fp.array(np_mat1)
np_u, np_s, np_v = sp_linalg.svds(np_mat1, k=10)
fp_u, fp_s, fp_v = fp.linalg.svds(fp_mat1, k=10)
verify(fp.absolute(fp_u), np.absolute(np_u))
verify(fp.absolute(fp_v), np.absolute(np_v))
assert all(abs(fp_s - np_s) / np_s < 1e-14)

print("test where")
np_mat1 = np.random.normal(scale=100, size=[1000, 10])
np_mat2 = np.random.normal(scale=100, size=[1000, 10])
fp_mat1 = fp.array(np_mat1)
fp_mat2 = fp.array(np_mat2)
np_res = np.where(np_mat1 > np_mat2, np_mat1, np_mat2)
fp_res = fp.where(fp_mat1 > fp_mat2, fp_mat1, fp_mat2)
verify(fp_res, np_res)

print("test average")
fp_res = fp.average(fp_mat1)
np_res = np.average(np_mat1)
tmp = np.array(fp_res, copy=True)
assert abs(tmp[0] - np_res) < 1e-12

fp_res = fp.average(fp_mat1, axis=0)
np_res = np.average(np_mat1, axis=0)
verify(fp_res, np_res)

fp_res = fp.average(fp_mat1, axis=1)
np_res = np.average(np_mat1, axis=1)
verify(fp_res, np_res)

np_weights = np.random.normal(scale=100, size=fp_mat1.shape)
fp_weights = fp.array(np_weights, dtype="l")
np_weights = np.array(fp_weights, copy=True)
fp_res = fp.average(fp_mat1, weights=fp_weights)
np_res = np.average(np_mat1, weights=np_weights)
tmp = np.array(fp_res, copy=True)
assert abs(tmp[0] - np_res) < 1e-11

np_weights = np.random.normal(scale=100, size=fp_mat1.shape[1])
fp_weights = fp.array(np_weights, dtype="l")
np_weights = np.array(fp_weights, copy=True)
fp_res = fp.average(fp_mat1, axis=1, weights=fp_weights, returned=True)
np_res = np.average(np_mat1, axis=1, weights=np_weights, returned=True)
verify(fp_res[0], np_res[0])

np_weights = np.random.normal(scale=100, size=fp_mat1.shape[0])
fp_weights = fp.array(np_weights, dtype="l")
np_weights = np.array(fp_weights, copy=True)
fp_res = fp.average(fp_mat1, axis=0, weights=fp_weights)
np_res = np.average(np_mat1, axis=0, weights=np_weights)
verify(fp_res, np_res)
