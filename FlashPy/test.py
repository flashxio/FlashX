import FlashPy
import numpy as np

def verify(fp_arr, np_arr):
    assert fp_arr.ndim == np_arr.ndim
    assert fp_arr.shape == np_arr.shape
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

print("test array")
np_mat1 = np.random.normal(size=[25, 10])
fp_mat1 = FlashPy.array(np_mat1, "d")
np_mat2 = np.random.normal(size=[25, 10])
fp_mat2 = FlashPy.array(np_mat2, "d")
verify(fp_mat1, np_mat1)
verify(fp_mat2, np_mat2)

print("test +")
fp_res = fp_mat1 + fp_mat2
np_res = np_mat1 + np_mat2
verify(fp_res, np_res)

print("test -")
fp_res = fp_mat1 - fp_mat2
np_res = np_mat1 - np_mat2
verify(fp_res, np_res)

print("test *")
fp_res = fp_mat1 * fp_mat2
np_res = np_mat1 * np_mat2
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
