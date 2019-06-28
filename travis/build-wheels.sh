#!/usr/bin/env bash
set -e -x

# Install a system package required by our library
yum install -y zlib-devel

# Compile wheels
for PYBIN in /opt/python/*/bin; do
	if [ `${PYBIN}/python -c "import sys; print(1 if sys.version_info[:2] >= (3, 5) else 0)"` -eq 1 ]; then
        "${PYBIN}/pip" install -U setuptools
        #"${PYBIN}/pip" install -U wheel
        "${PYBIN}/pip" wheel /io/ -w wheelhouse/
	fi
done

# Bundle external shared libraries into the wheels
for whl in wheelhouse/*.whl; do
    auditwheel repair "$whl" --plat $PLAT -w /io/wheelhouse/
done

# Install packages and test
#for PYBIN in /opt/python/*/bin/; do
#    "${PYBIN}/pip" install python-manylinux-demo --no-index -f /io/wheelhouse
#    (cd "$HOME"; "${PYBIN}/nosetests" pymanylinuxdemo)
#done
