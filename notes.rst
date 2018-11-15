A few notes:
************

Building wheels
===============

Binary wheels for linux need to be built with manylinux1 image, i.e. a Centos5
base image. Easiest is for now to build this with docker (see `https://github.com/pypa/manylinux`_):

.. sh:

    docker pull quay.io/pypa/manylinux1_x86_64

    docker run --rm -v `pwd`:/io quay.io/pypa/manylinux1_x86_64 /io/build-linux-wheels.sh

This will compile wheels and store them under wheelhouse. upload them with twine to pypi.

on commit id 47dda7b this has been used like this on travis. checkout this release 
for reference.


Building wheels for windows and mac:
------------------------------------

This has not yet been fully setteled. MacPython seems to provide a good way how to do
it for macos, e.g. `https://github.com/MacPython/pytables-wheels`_. However, it seems
that `https://github.com/matthew-brett/multibuild`_ provides even a better template
how to build wheels also for windows. So far this has not been investigated further.

