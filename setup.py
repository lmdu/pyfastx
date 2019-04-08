#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os
import sys
from distutils.core import setup, Extension
#from Cython.Build import cythonize

sources = ['src/module.c', 'src/fastx.c', 'src/zran.c', 'src/util.c']

extensions = [
	Extension('pyfastx', sources, extra_link_args=['-lz', '-lsqlite3']),
]

setup(
	name = 'pyfastx',
	version = '1.0.0',
	ext_modules = extensions
)
