#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os
import sys
import glob
import platform
#import unittest##
#from distutils.core import setup, Extension
from setuptools import setup, Extension

#link_args = ['-lz', '-lsqlite3']
comp_args = []
include_dirs = []
libs = []
lib_dirs = []

#if os.name == 'nt' and '64' in platform.architecture()[0]:
    #link_args.append('-DMS_WIN64')
    #comp_args.append('-DMS_WIN64')
    #comp_args.append('-D_FILE_OFFSET_BITS=64')
    #comp_args.append('-pedantic')
    #comp_args.append('-Wno-unused-function')

if os.name == 'nt':
	zlib_home = os.environ.get('ZLIB_HOME')
	include_dirs.append(os.path.join(zlib_home, 'include'))
	libs.append('zlib')
	libs.append('sqlite3')
	lib_dirs.append(os.path.join(zlib_home, 'lib'))

	if sys.version_info[0] == 2:
		include_dirs.append('compat')
		comp_args.append('-DNO_C99')

extension = Extension('pyfastx',
    sources = glob.glob('src/*.c'),
    extra_compile_args = comp_args,
    #extra_link_args = link_args,
    include_dirs = include_dirs,
    libraries = libs,
    library_dirs = library_dirs
    #define_macros = [("_FILE_OFFSET_BITS", 64)]
)

description = (
    "pyfastx is a python module for fast random "
    "access to sequences from plain and gzipped "
    "FASTA/Q file"
)

with open('README.rst') as fh:
    long_description = fh.read()

with open(os.path.join('src', 'version.h')) as fh:
    version = fh.read().split()[2].strip('"')

#def my_test_suite():
#    test_loader = unittest.TestLoader()
#    test_suite = test_loader.discover('tests', pattern='test_*.py')
#    return test_suite

setup(
    name = 'pyfastx',
    version = version,
    ext_modules = [extension],
    description = description,
    long_description = long_description,
    author = 'Lianming Du',
    author_email = 'adullb@qq.com',
    url = 'https://github.com/lmdu/pyfastx',
    license = 'MIT',
    keywords = 'fasta sequence bioinformatics',
    classifiers = [
            "Development Status :: 5 - Production/Stable",
            "Intended Audience :: Developers",
            "Intended Audience :: Education",
            "Intended Audience :: Science/Research",
            "Natural Language :: English",
            "License :: OSI Approved :: MIT License",
            "Programming Language :: C",
            "Programming Language :: Python :: 2.7",
            "Programming Language :: Python :: 3.4",
            "Programming Language :: Python :: 3.5",
            "Programming Language :: Python :: 3.6",
            "Programming Language :: Python :: 3.7",
            "Operating System :: Microsoft :: Windows",
            "Operating System :: POSIX :: Linux",
            "Operating System :: Unix",
            "Topic :: Scientific/Engineering :: Bio-Informatics"
    ],
    entry_points = {
        'console_scripts': ['pyfastx = pyfastxcli:main']
    },
    py_modules = ["pyfastxcli"],
    test_suite = "tests"
)
