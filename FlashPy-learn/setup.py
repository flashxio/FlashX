from distutils.core import setup, Extension
from Cython.Build import cythonize
import os 
dir_path = os.path.dirname(os.path.realpath(__file__))

setup(
    name="fplearn",
    version="0.0.1",
    description="A parallel and scalable machine learning library",
    long_description="This parallelizes and scales the machine learning algorithms, " +\
            "in Scikit-learn. It extends " +\
            "memory capacity with SSDs and is optimized for NUMA machines",
    url="https://github.com/flashxio/FlashX",
    author="Da Zheng",
    author_email="dzheng5@jhu.edu",
    license="Apache License, Version 2.0",
    keywords="parallel scalable machine-learning Scikit-learn",
    install_requires=[
        "numpy",
        "Cython==0.23.5",
        "cython==0.23.5",
        ],
    package_dir = {"fplearn": os.path.join("fplearn")},
    packages=["fplearn", "fplearn.cluster",
        "fplearn.externals", "fplearn.linear_model",
        "fplearn.metrics", "fplearn.utils"],
#    libraries =libraries,
#    cmdclass = {'build_clib': flashpy_clib, 'build_ext': build_ext},
)
