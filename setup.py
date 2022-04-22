from setuptools import setup, find_packages
from setuptools.extension import Extension
from Cython.Build import cythonize
import numpy
import os

package_name = 'pyopa'

c_dir = 'c_source/'
c_sources = ['Python_extension.c', 'DynProgr_sse_short.c', 'DynProgr_sse_byte.c',
             'DynProgr_scalar.c', 'EstimatePam.c', 'Page_size.c', 'DynProgr_sse_double.c']

__version__ = "Undefined"
for line in open('{}/__init__.py'.format(package_name)):
    if line.startswith('__version__'):
        exec(line.strip())

# Get the long description from the README file
with open('README.rst') as f:
    long_description = f.read()

setup(
    ext_modules=cythonize([
        Extension(
            package_name + ".backend.pyopa",
            sources=[os.path.join(package_name, 'backend', 'pyopa.pyx')] +
                    list(map(lambda c: os.path.join(c_dir, c), c_sources)),
            extra_compile_args=['-O2', '-DSIMDE_ENABLE_NATIVE_ALIASES'], #, '-DPY_DEBUG'],
            include_dirs=[c_dir, numpy.get_include(), '.']
        )]
    ),
    name=package_name,
    version=__version__,
    author='OMA Browser',
    author_email='contact@omabrowser.org',
    install_requires=['numpy'],
    packages=find_packages(), 
    include_package_data=True,
    license='MPL 2.0',
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Developers',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'License :: OSI Approved :: Mozilla Public License 2.0 (MPL 2.0)',
        'Programming Language :: Python :: 2',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
    ],
    keywords='sequence alignments Smith-Waterman Needleman-Wunsch dynamic programming bioinformatics',
    description='PyOPA - optimal pairwise sequence alignments',
    long_description=long_description,
    url='http://omabrowser.org/',
)
