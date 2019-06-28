#!/bin/bash
set -e -x

yum install -y atlas-devel zlib-devel

for PYBIN in /opt/python/*/bin; do
	"${PYBIN}/pip" install nose
	"${PYBIN}/pip" wheel /io/ -w wheelhouse/
done

for whl in wheelhouse/*.whl; do
	auditwheel repair "$whl" --plat $PLAT -w /io/wheelhouse/
done

#chown -R --reference=/io/setup.py /io/wheelhouse

#for PYBIN in /opt/python/*/bin; do
#	"${PYBIN}/pip" install pyfastx --no-index -f /io/wheelhouse
#	(cd "$HOME";)
#done
