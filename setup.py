import os
import sys
import glob
import tarfile
import zipfile
import urllib.request
from setuptools import setup, Extension
from setuptools.command import build_ext

root_dir = os.path.dirname(os.path.abspath(__file__))

sources = glob.glob(os.path.join(root_dir, 'src', '*.c'))
link_args = []
comp_args = []
include_dirs = []

def prepare_zlib():
    global include_dirs
    global sources

    zlib_dir = os.path.join(root_dir, "zlib-1.3.1")
    zlib_file = os.path.join(root_dir, "zlib-1.3.1.zip")
    url = "https://github.com/madler/zlib/releases/download/v1.3.1/zlib131.zip"

    if not os.path.exists(zlib_dir):
        if not os.path.isfile(zlib_file):
            urllib.request.urlretrieve(url, zlib_file)

        with zipfile.ZipFile(zlib_file) as _zip:
            _zip.extractall()

    include_dirs.append(zlib_dir)
    sources.extend(glob.glob(os.path.join(zlib_dir, '*.c')))

def prepare_sqlite3():
    global include_dirs
    global sources

    sqlite_dir = os.path.join(root_dir, "sqlite-amalgamation-3470200")
    sqlite_file = os.path.join(root_dir, "sqlite-amalgamation-3470200.zip")
    url = "https://www.sqlite.org/2024/sqlite-amalgamation-3470200.zip"

    if not os.path.exists(sqlite_dir):
        if not os.path.isfile(sqlite_file):
            urllib.request.urlretrieve(url, sqlite_file)

        with zipfile.ZipFile(sqlite_file) as _zip:
            _zip.extractall()

    include_dirs.append(sqlite_dir)
    sources.append(os.path.join(sqlite_dir, 'sqlite3.c'))

def prepare_indexed_gzip():
    global include_dirs
    global sources

    igzip_dir = os.path.join(root_dir, "indexed_gzip-1.9.4", "indexed_gzip")
    igzip_file = os.path.join(root_dir, "indexed_gzip-1.9.4.zip")
    url = "https://github.com/pauldmccarthy/indexed_gzip/archive/refs/tags/v1.9.4.zip"

    if not os.path.exists(igzip_dir):
        if not os.path.isfile(igzip_file):
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

    if sys.platform.startswith('darwin'):
        comp_args.append('-DHAVE_UNISTD_H')
    else:
        comp_args.extend([
            '-D_LFS64_LARGEFILE',
            '-D_LARGEFILE64_SOURCE'
        ])

class CustomBuildExt(build_ext.build_ext):
    def build_extensions(self):
        if not sys.platform.startswith(('win', 'darwin')):
            if '--debug' in sys.argv:
                self.compiler.compiler_so = [
                    opt
                    for opt in self.compiler.compiler_so
                    if opt != '-O2'
                ]
                self.compiler.compiler_so.append('-O0')
                self.compiler.compiler_so.append('--coverage')

                self.compiler.linker_so = [
                    opt
                    for opt in self.compiler.linker_so
                    if opt != '-O2'
                ]
                self.compiler.linker_so.append('-O0')
                self.compiler.linker_so.append('--coverage')

            else:
                self.compiler.compiler_so = [
                    opt
                    for opt in self.compiler.compiler_so
                    if opt != '-g'
                ]

                self.compiler.linker_so = [
                    opt
                    for opt in self.compiler.linker_so
                    if opt != '-g'
                ]

        super().build_extensions()

#if not sys.platform.startswith('linux'):
prepare_sqlite3()
prepare_zlib()
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
            "Programming Language :: Python :: 3.12",
            "Programming Language :: Python :: 3.13",
            "Operating System :: Microsoft :: Windows",
            "Operating System :: POSIX :: Linux",
            "Operating System :: Unix",
            "Topic :: Scientific/Engineering :: Bio-Informatics"
    ],
    entry_points = {
        'console_scripts': ['pyfastx = pyfastxcli:main']
    },
    py_modules = ["pyfastxcli"],
    test_suite = "tests",
    cmdclass = {
        'build_ext': CustomBuildExt
    }
)
