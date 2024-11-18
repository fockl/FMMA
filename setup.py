from glob import glob
from setuptools import setup
from pybind11.setup_helpers import Pybind11Extension

ext_modules = [
    Pybind11Extension("pyfmma", ["cpython/pyfmma.cpp"] + sorted(glob("src/*.cpp")))
    ]

setup(
    name="pyfmma",
    ext_modules=ext_modules,
    install_requires=[])
