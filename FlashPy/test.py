import FlashPy
import numpy as np

def verify(fp_arr, np_arr):
    assert fp_arr.ndim == np_arr.ndim
    assert fp_arr.shape == np_arr.shape
    assert fp_arr.size == np_arr.size
    assert fp_arr.itemsize == np_arr.itemsize
    assert fp_arr.nbytes == np_arr.nbytes
    tmp = np.array(fp_arr, copy=True)
    assert np.all(tmp == np_arr)

def verify_init(dtype):
    print("verify " + dtype)
    fp_mat1 = FlashPy.empty([25, 10], dtype=dtype)
    assert fp_mat1.shape[0] == 25
    assert fp_mat1.shape[1] == 10
    assert fp_mat1.dtype == dtype

    fp_mat2 = FlashPy.empty_like(fp_mat1)
    assert fp_mat1.shape == fp_mat2.shape
    assert fp_mat1.dtype == fp_mat2.dtype

    fp_mat1 = FlashPy.ones([25, 10], dtype=dtype)
    assert fp_mat1.shape[0] == 25
    assert fp_mat1.shape[1] == 10
    assert fp_mat1.dtype == dtype

    fp_mat1 = FlashPy.zeros([25, 10], dtype=dtype)
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

print("test C-contiguous array")
np_mat1 = np.random.normal(scale=100, size=[10, 10])
fp_mat1 = FlashPy.array(np_mat1)
verify(fp_mat1, np_mat1)
np_mat1 = np.array(fp_mat1, copy=True)
verify(fp_mat1, np_mat1)

print("test F-contiguous array")
np_mat1 = np.array(np_mat1, order='F', copy=True)
fp_mat1 = FlashPy.array(np_mat1)
verify(fp_mat1, np_mat1)
np_mat1 = np.array(fp_mat1, copy=True)
verify(fp_mat1, np_mat1)

np_mat1 = np.random.normal(scale=100, size=[25, 10])
fp_mat1 = FlashPy.array(np_mat1)
np_mat2 = np.random.normal(scale=100, size=[25, 10])
fp_mat2 = FlashPy.array(np_mat2)

print("test +")
fp_res = fp_mat1 + fp_mat2
np_res = np_mat1 + np_mat2
verify(fp_res, np_res)

print("test -")
fp_res = fp_mat1 - fp_mat2
np_res = np_mat1 - np_mat2
verify(fp_res, np_res)

print("test *")
fp_res = FlashPy.array(np_mat1, "l") * FlashPy.array(np_mat2, "l")
np_res = np.array(np_mat1, "l") * np.array(np_mat2, "l")
verify(fp_res, np_res)

print("test /")
fp_res = fp_mat1 / fp_mat2
np_res = np_mat1 / np_mat2
verify(fp_res, np_res)

print("test //")
fp_res = fp_mat1 // fp_mat2
np_res = np_mat1 // np_mat2
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
    fp_mat1 = FlashPy.empty([25, 10], dtype=dtype1)
    assert fp_mat1.dtype == dtype1
    fp_mat2 = fp_mat1.cast_ele_type(dtype2)
    assert fp_mat2.dtype == dtype2

verify_cast("f", "d")
verify_cast("d", "f")
verify_cast("b", "h")
verify_cast("i", "l")

print("test sum")
np_mat1 = np.random.normal(scale=100, size=[25, 10])
fp_mat1 = FlashPy.array(np_mat1, dtype="l")
np_mat1 = np.array(np_mat1, dtype="l")
fp_res = FlashPy.sum(fp_mat1)
np_res = np.sum(np_mat1)
tmp = np.array(fp_res, copy=True)
assert abs(tmp - np_res) < 1e-13

fp_res = FlashPy.sum(fp_mat1, axis=0)
np_res = np.sum(np_mat1, axis=0)
verify(fp_res, np_res)

fp_res = FlashPy.sum(fp_mat1, axis=0, dtype="d")
np_res = np.sum(np_mat1, axis=0, dtype="d")
verify(fp_res, np_res)

fp_res = FlashPy.sum(fp_mat1, axis=1)
np_res = np.sum(np_mat1, axis=1)
verify(fp_res, np_res)

print("test mean")
fp_res = FlashPy.mean(fp_mat1)
np_res = np.mean(np_mat1)
tmp = np.array(fp_res, copy=True)
assert tmp[0] == np_res

fp_res = FlashPy.mean(fp_mat1, axis=0)
np_res = np.mean(np_mat1, axis=0)
verify(fp_res, np_res)

fp_res = FlashPy.mean(fp_mat1, axis=1)
np_res = np.mean(np_mat1, axis=1)
verify(fp_res, np_res)

print("test dot on matrix")
np_mat1 = np.random.normal(scale=100, size=[25, 10])
fp_mat1 = FlashPy.array(np_mat1, dtype="l")
np_mat1 = np.array(np_mat1, dtype="l")
fp_res = FlashPy.dot(fp_mat1, fp_mat1.transpose())
np_res = np.dot(np_mat1, np_mat1.transpose())
verify(fp_res, np_res)

print("test dot on vector")
np_mat1 = np.random.normal(scale=100, size=25)
fp_mat1 = FlashPy.array(np_mat1, dtype="l")
np_mat1 = np.array(np_mat1, dtype="l")
fp_res = FlashPy.dot(fp_mat1, fp_mat1)
np_res = np.dot(np_mat1, np_mat1)
assert np.array(fp_res, copy=True) == np_res

print("test average")
fp_res = FlashPy.average(fp_mat1)
np_res = np.average(np_mat1)
tmp = np.array(fp_res, copy=True)
assert tmp[0] == np_res

fp_res = FlashPy.average(fp_mat1, axis=0)
np_res = np.average(np_mat1, axis=0)
verify(fp_res, np_res)

fp_res = FlashPy.average(fp_mat1, axis=1)
np_res = np.average(np_mat1, axis=1)
verify(fp_res, np_res)

np_weights = np.random.normal(scale=100, size=[25, 10])
fp_weights = FlashPy.array(np_weights, dtype="l")
np_weights = np.array(fp_weights, copy=True)
fp_res = FlashPy.average(fp_mat1, weights=fp_mat1)
np_res = np.average(np_mat1, weights=np_mat1)
tmp = np.array(fp_res, copy=True)
assert tmp[0] == np_res

np_weights = np.random.normal(scale=100, size=10)
fp_weights = FlashPy.array(np_weights, dtype="l")
np_weights = np.array(fp_weights, copy=True)
fp_res = FlashPy.average(fp_mat1, axis=0, weights=fp_weights)
np_res = np.average(np_mat1, axis=0, weights=np_weights)
verify(fp_res, np_res)

np_weights = np.random.normal(scale=100, size=25)
fp_weights = FlashPy.array(np_weights, dtype="l")
np_weights = np.array(fp_weights, copy=True)
fp_res = FlashPy.average(fp_mat1, axis=1, weights=fp_weights)
np_res = np.average(np_mat1, axis=1, weights=np_weights)
verify(fp_res, np_res)
