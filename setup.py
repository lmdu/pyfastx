import os
import sys
import glob
import tarfile
import zipfile
import urllib.request
from setuptools import setup, Extension

root_dir = os.path.dirname(os.path.abspath(__file__))

sources = glob.glob(os.path.join(root_dir, 'src', '*.c'))
link_args = []
comp_args = []
include_dirs = []

def prepare_zlib():
    global include_dirs
    global sources

    zlib_dir = os.path.join(root_dir, "zlib-1.2.13")
    zlib_file = os.path.join(root_dir, "zlib-1.2.13.zip")
    url = "https://github.com/madler/zlib/releases/download/v1.2.13/zlib1213.zip"

    if not os.path.exists(zlib_dir):
        urllib.request.urlretrieve(url, zlib_file)
        with zipfile.ZipFile(zlib_file) as _zip:
            _zip.extractall()

    include_dirs.append(zlib_dir)
    sources.extend(glob.glob(os.path.join(zlib_dir, '*.c')))

def prepare_sqlite3():
    global include_dirs
    global sources

    sqlite_dir = os.path.join(root_dir, "sqlite-amalgamation-3400100")
    sqlite_file = os.path.join(root_dir, "sqlite-amalgamation-3400100.zip")
    url = "https://www.sqlite.org/2022/sqlite-amalgamation-3400100.zip"

    if not os.path.exists(sqlite_dir):
        urllib.request.urlretrieve(url, sqlite_file)
        with zipfile.ZipFile(sqlite_file) as _zip:
            _zip.extractall()

    include_dirs.append(sqlite_dir)
    sources.append(os.path.join(sqlite_dir, 'sqlite3.c'))

def prepare_indexed_gzip():
    global include_dirs
    global sources

    igzip_dir = os.path.join(root_dir, "indexed_gzip-1.7.0", "indexed_gzip")
    igzip_file = os.path.join(root_dir, "indexed_gzip-1.7.0.zip")
    url = "https://github.com/pauldmccarthy/indexed_gzip/archive/refs/tags/v1.7.0.zip"

    if not os.path.exists(igzip_dir):
        urllib.request.urlretrieve(url, igzip_file)
        with zipfile.ZipFile(igzip_file) as _zip:
            _zip.extractall()

    include_dirs.append(igzip_dir)
    sources.extend(glob.glob(os.path.join(igzip_dir, 'zran*.c')))


if sys.platform.startswith('win'):
    comp_args.extend([
        '/D_LFS64_LARGEFILE',
        '/D_LARGEFILE64_SOURCE',
        '/D_FILE_OFFSET_BITS=64'
    ])
else:
    comp_args.extend([
        '-Wno-unused-result',
        '-D_FILE_OFFSET_BITS=64'
    ])

    if sys.platform.startswith('linux'):
        link_args.extend(['-lz', '-lsqlite3'])
        comp_args.extend([
            '-D_LFS64_LARGEFILE',
            '-D_LARGEFILE64_SOURCE',
        ])

    elif sys.platform.startswith('darwin'):
        comp_args.append('-DHAVE_UNISTD_H')

if not sys.platform.startswith('linux'):
    prepare_zlib()
    prepare_sqlite3()

prepare_indexed_gzip()

extension = Extension('pyfastx',
    sources = sources,
    include_dirs = include_dirs,
    extra_compile_args = comp_args,
    extra_link_args = link_args
)

description = (
    "pyfastx is a python module for fast random "
    "access to sequences from plain and gzipped "
    "FASTA/Q file"
)

with open(os.path.join(root_dir, 'README.rst')) as fh:
    long_description = fh.read()

with open(os.path.join(root_dir, 'src', 'version.h')) as fh:
    version = fh.read().split()[2].strip('"')

setup(
    name = 'pyfastx',
    version = version,
    ext_modules = [extension],
    description = description,
    long_description = long_description,
    #long_description_content_type = 'text/x-rst',
    author = 'Lianming Du',
    author_email = 'adullb@qq.com',
    url = 'https://github.com/lmdu/pyfastx',
    license = 'MIT',
    keywords = 'fasta fastq sequence bioinformatics',
    classifiers = [
            "Development Status :: 5 - Production/Stable",
            "Intended Audience :: Developers",
            "Intended Audience :: Education",
            "Intended Audience :: Science/Research",
            "Natural Language :: English",
            "License :: OSI Approved :: MIT License",
            "Programming Language :: C",
            "Programming Language :: Python :: 3.6",
            "Programming Language :: Python :: 3.7",
            "Programming Language :: Python :: 3.8",
            "Programming Language :: Python :: 3.9",
            "Programming Language :: Python :: 3.10",
            "Programming Language :: Python :: 3.11",
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
