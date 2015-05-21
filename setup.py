from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize
import numpy
import sys
import os


package_name = 'pyopa'

c_dir = 'c_source/'
c_sources = ['Python_extension.c', 'DynProgr_sse_short.c', 'DynProgr_sse_byte.c',
             'DynProgr_scalar.c', 'EstimatePam.c', 'Page_size.c', 'DynProgr_sse_double.c']

data_dir = os.path.join(sys.prefix, package_name + '_test')


setup(
    ext_modules = cythonize([
	    Extension(
		package_name,
		sources=[package_name+'.pyx'] + map(lambda c: os.path.join(c_dir, c), c_sources),
		include_dirs=[c_dir, numpy.get_include()]
	)]
    ),
    name=package_name,
    version='0.6',
    author='OMA Browser',
    author_email='contact@omabrowser.org',
    install_requires=['numpy', 'cython'],
    data_files=[
        (os.path.join(data_dir, 'data'), ['test/data/reference_test_results.dat', 'test/data/testseqs.txt']),
        (os.path.join(data_dir, 'data/matrices/json'),
         ['test/data/matrices/json/all_matrices.json', 'test/data/matrices/json/logPAM1.json']),
        (data_dir, ['test/runall.py', 'test/test_darwin.py', 'test/test_input.py'])
    ],
    license='MPL 2.0',
    classifiers = [
         'Development Status :: 3 - Alpha',
         'Intended Audience :: Developers',
         'Intended Audience :: Science/Research',
         'Topic :: Scientific/Engineering :: Bio-Informatics',
         'License :: OSI Approved :: Mozilla Public License 2.0 (MPL 2.0)',
         'Programming Language :: Python :: 2',
         'Programming Language :: Python :: 2.7',
         'Programming Language :: Python :: 3',
         'Programming Language :: Python :: 3.3',
         'Programming Language :: Python :: 3.4',
         ],
    keywords='sequence alignments Smith-Waterman dynamic programming',
    description='pyopa computes optimal pairwise sequence alignments using a vectorized implementation of the Smith-Waterman algorithm',
    url='http://omabrowser.org/',
)
