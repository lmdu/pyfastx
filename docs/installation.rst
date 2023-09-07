Installation
============

You can install pyfastx via the Python Package Index (PyPI) (**recommended**) or from source.

Make sure you have installed both `pip <https://pip.pypa.io/en/stable/installing/>`_ and Python before starting.

Currently, ``pyfastx`` supports Python 3.6, 3.7, 3.8, 3.9, 3.10, 3.11 and can work on Windows, Linux, MacOS.

Install from PyPI
-----------------

::

	pip install pyfastx

Update pyfastx using ``pip``

::

	pip install -U pyfastx

Install from source
-------------------

``pyfastx`` depends on `zlib <https://zlib.net/>`_ and `sqlite3 <https://www.sqlite.org/index.html>`_. If you want to compile and install pyfastx from source code. First, you should install zlib and sqlite3.

On Centos

::

	yum install zlib-devel
	yum install sqlite-devel

On Ubuntu

::

	apt install zlib1g-dev
	apt install libsqlite3-dev

On Windows and MacOS, it will automatically download zlib and sqlite3 library to build.


Second, clone pyfastx using ``git`` or download latest `release <https://github.com/lmdu/pyfastx/releases>`_:

::

	git clone https://github.com/lmdu/pyfastx.git

Then, ``cd`` to the pyfastx folder and run install command:

::

	cd pyfastx
	python setup.py install
