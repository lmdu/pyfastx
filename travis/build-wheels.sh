#!/bin/bash
set -e -x

yum install -y zlib-devel

for PYBIN in /opt/python/*/bin; do
	"${PYBIN}/pip" install --upgrade pip wheel
	"${PYBIN}/pip" install --upgrade setuptools
	"${PYBIN}/pip" wheel /io/ -w wheelhouse/
done

for whl in wheelhouse/*.whl; do
	auditwheel repair "$whl" --plat $PLAT -w /io/wheelhouse/
done

#for PYBIN in /opt/python/*/bin; do
#	"${PYBIN}/pip" install pyfastx --no-index -f /io/wheelhouse
#	(cd "$HOME";)
#done
