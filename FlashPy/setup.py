from distutils.core import setup, Extension
from Cython.Build import cythonize
import os 
dir_path = os.path.dirname(os.path.realpath(__file__))
print(dir_path)

setup(ext_modules = cythonize(Extension(
           "FlashPy",                             # the extension name
           sources=["FlashPy.pyx"],               # the Cython source and
                                                  # additional C++ source files
           include_dirs = ["../matrix", "../libsafs"],
           libraries = ["hwloc", "cblas"],
           language="c++",                        # generate and compile C++ code
           extra_compile_args=['-fopenmp', '-std=c++11'],
           extra_link_args=['-fopenmp'],
           extra_objects = ["../build/matrix/libFMatrix.a", "../build/libsafs/libsafs.a"],
      )))
