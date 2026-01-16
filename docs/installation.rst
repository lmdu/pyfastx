Installation
============

You can install pyfastx via the Python Package Index (PyPI) (**recommended**) or from source.

Make sure you have installed both `pip <https://pip.pypa.io/en/stable/installing/>`_ and Python before starting.

Currently, ``pyfastx`` supports Python 3.8, 3.9, 3.10, 3.11, 3.12, 3.13, 3.14 and can work on Windows, Linux, MacOS.

Install from PyPI
-----------------

::

	pip install pyfastx

Update pyfastx using ``pip``

::

	pip install -U pyfastx

Install from source
-------------------

``pyfastx`` depends on `zlib <https://zlib.net/>`_, `sqlite3 <https://www.sqlite.org/index.html>`_ and `indexed_gzip <https://github.com/pauldmccarthy/indexed_gzip>`_. In latest version, pyfastx will automatically download these libraries to build.


First, clone pyfastx using ``git`` or download latest `release <https://github.com/lmdu/pyfastx/releases>`_:

::

	git clone https://github.com/lmdu/pyfastx.git

Then, ``cd`` to the pyfastx folder and run install command:

::

	cd pyfastx
	python setup.py install

Or just build:

::

	cd pyfastx
	python setup.py build_ext -i
