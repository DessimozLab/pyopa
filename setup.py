from setuptools import setup, Extension
from Cython.Build import cythonize
import numpy
import sys
import os


package_name = 'pyopa'

c_dir = 'c_source/'
c_sources = ['Python_extension.c', 'DynProgr_sse_short.c', 'DynProgr_sse_byte.c',
             'DynProgr_scalar.c', 'EstimatePam.c', 'Page_size.c', 'DynProgr_sse_double.c']

data_dir = os.path.join(sys.prefix, package_name + '_test')

here = os.path.abspath(os.path.dirname(__file__))
# Get the long description from the README file
with open(os.path.join(here, 'README.rst')) as f:
    long_description = f.read()

setup(
    ext_modules = cythonize([
	    Extension(
		package_name,
		sources=[package_name+'.pyx'] + list(map(lambda c: os.path.join(c_dir, c), c_sources)),
		#libraries=["swps3"],
		#extra_link_args=["-Wl,-stack_size", '-Wl,0x10000000'],
                #extra_compile_args = ['-DPY_DEBUG'],
		include_dirs=[c_dir, numpy.get_include()]
	)]
    ),
    name=package_name,
    version='0.7.0',
    author='OMA Browser',
    author_email='contact@omabrowser.org',
    install_requires=['numpy', 'cython'],
    data_files=[
        (os.path.join(data_dir, 'data'), ['test/data/reference_test_results.dat', 'test/data/testseqs.txt']),
        (os.path.join(data_dir, 'data/matrices/json'),
         ['test/data/matrices/json/all_matrices.json', 'test/data/matrices/json/logPAM1.json']),
        (data_dir, ['test/test_darwin.py', 'test/test_input.py', 'test/test_env.py'])
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
    keywords='sequence alignments Smith-Waterman Needleman-Wunsch dynamic programming bioinformatics',
    description='PyOPA - optimal pairwise sequence alignments',
    long_description=long_description,
    
    url='http://omabrowser.org/',
)
