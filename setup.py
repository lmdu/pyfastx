#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os
import sys
import glob
#from distutils.core import setup, Extension
from setuptools import setup, Extension

sources = glob.glob('src/*.c')

extensions = [
    Extension('pyfastx', sources, extra_link_args=['-lz', '-lsqlite3'])
]

description = (
    "pyfastx is a python module for fast random "
    "access to sequences from plain and gzipped "
    "FASTA file"
)

with open('README.rst') as fh:
    long_description = fh.read()

with open(os.path.join('src', 'version.h')) as fh:
    version = fh.read().split()[2].strip('"')

setup(
    name = 'pyfastx',
    version = version,
    ext_modules = extensions,
    description = description,
    long_description = long_description,
    author = 'Lianming Du',
    author_email = 'adullb@qq.com',
    url = 'https://github.com/lmdu/pyfastx',
    license = 'MIT',
    keywords = 'fasta sequence bioinformatics',
    classifiers = [
            "Development Status :: 3 - Alpha",
            "Intended Audience :: Developers",
            "Intended Audience :: Education",
            "Intended Audience :: Science/Research",
            "Natural Language :: English",
            "License :: OSI Approved :: MIT License",
            "Programming Language :: C",
            "Programming Language :: Python :: 3.5",
            "Programming Language :: Python :: 3.6",
            "Programming Language :: Python :: 3.7",
            "Operating System :: Microsoft :: Windows",
            "Operating System :: POSIX :: Linux",
            "Operating System :: Unix",
            "Topic :: Scientific/Engineering :: Bio-Informatics"
    ],
    test_suite="tests"
)
