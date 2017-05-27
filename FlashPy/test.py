import FlashPy
import numpy as np

def verify(fp_arr, np_arr):
    assert fp_arr.ndim == np_arr.ndim
    assert fp_arr.shape == np_arr.shape
    tmp = np.array(fp_arr, copy=True)
    assert np.all(tmp == np_arr)


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
