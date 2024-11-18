from glob import glob
from setuptools import setup
from pybind11.setup_helpers import Pybind11Extension

ext_modules = [
    Pybind11Extension("pyfmma", ["cpython/pyfmma.cpp"] + sorted(glob("src/*.cpp")))
    ]

URL = 'https://github.com/fockl/FMMA'
AUTHOR = 'fockl'
AUTHOR_EMAIL = 'j4xn3gvbxqtfdnxg@gmail.com'

with open('README.md', encoding='utf-8') as f:
  readme = f.read()
LONG_DESCRIPTION = readme

setup(
    name="pyfmma",
    version="0.1.1",
    description="pyfmma: Fast Multipole Method for arbitrary functions",
    long_description=LONG_DESCRIPTION,
    long_description_content_type='text/markdown',
    author=AUTHOR,
    author_email=AUTHOR_EMAIL,
    url=URL,
    download_url=URL,
    license='MIT',
    ext_modules=ext_modules,
    classifiers=[
      'License :: OSI Approved :: MIT License',
      'Programming Language :: C++',
      'Programming Language :: Python',
      ],
    install_requires=['setuptools<72.2.0', 'pybind11'],
    )
