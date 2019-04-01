#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os
import sys
from distutils.core import setup, Extension
#from Cython.Build import cythonize

extensions = [
	Extension('pyfastx', ['pyfastx.c', 'zran.c'], extra_link_args=['-lz', '-lsqlite3']),
]

setup(
	name = 'pyfastx',
	version = '1.0.0',
	ext_modules = extensions
)
