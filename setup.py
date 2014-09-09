from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize
import numpy
import sys
import os


package_name = 'pyopa'

c_dir = 'c_source/'
c_sources = ['Python_extension.c', 'DynProgr_sse_short.c', 'DynProgr_sse_byte.c',
             'DynProgr_scalar.c', 'EstimatePam.c', 'Page_size.c']

data_dir = os.path.join(sys.prefix, package_name + '_test')


setup(
    ext_modules = cythonize([
	    Extension(
		package_name,
		sources=[package_name+'.pyx'] + map(lambda c: os.path.join(c_dir, c), c_sources),
		#libraries=["swps3"],
		#extra_link_args=["-L../"],
        #extra_compile_args = ['-DPY_DEBUG'],
		include_dirs=[c_dir, numpy.get_include()]
	)]
    ),
    name=package_name,
    version='0.6',
    author='Your Name',
    author_email='your@email',
    #cmdclass = {'build_ext': build_ext},
    #install_requires=['numpy', 'cython'],
    data_files=[
        (os.path.join(data_dir, 'data'), ['test/data/reference_test_results.dat', 'test/data/testseqs.txt']),
        (os.path.join(data_dir, 'data/matrices/json'),
         ['test/data/matrices/json/all_matrices.json', 'test/data/matrices/json/logPAM1.json']),
        (data_dir, ['test/runall.py', 'test/test_darwin.py', 'test/test_input.py'])
    ],
    license='MIT',
    description='Example package that says hello'
)
