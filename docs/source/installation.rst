Installation
============

Dependencies
------------

Make sure you have both `pip <https://pip.pypa.io/en/stable/installing/>`_ and at least version 3.4 of Python before starting.

Currently, pyfastx only support Python 3.4, 3.5, 3.6, 3.7 and can work on Windows, Linux, MacOS.

.. note::
	
	Python 2.7 will not be maintained past 2020.

pyfastx depends on `zlib <https://zlib.net/>`_ and `sqlite3 <https://www.sqlite.org/index.html>`_. If you want to compile and install pyfastx from source code, you should install zlib and sqlite3 first.

On Centos 

::

	yum install zlib-devel
	yum install sqlite-devel

On Ubuntu

::

	apt install zlib1g-dev
	apt install libsqlite3-dev

On MacOS

::

	brew install zlib
	brew install sqlite3


Installation
------------

You can install pyfastx via the Python Package Index (PyPI) (**recommended**) or from source.

Install from PyPI
-----------------

::

	pip install pyfastx

Update pyfastx using ``pip``

::

	pip install -U pyfastx

Install from source
-------------------

First, clone pyfastx using ``git`` or download latest `release <https://github.com/lmdu/pyfastx/releases>`_:

::

	git clone https://github.com/lmdu/pyfastx.git

Then, ``cd`` to the pyfastx folder and run install command:

::

	cd pyfastx
	python setup.py install
