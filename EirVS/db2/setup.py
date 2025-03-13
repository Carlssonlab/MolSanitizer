from setuptools import setup
from Cython.Build import cythonize

setup(name="mol2db2", ext_modules=cythonize(["buckets2.py", "clash.py"]))
