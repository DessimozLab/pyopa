from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize
import numpy
import subprocess


#subprocess.Popen("rm -rf build", shell=True, executable="/bin/bash")
#subprocess.Popen("rm -rf *.c", shell=True, executable="/bin/bash")
#subprocess.Popen("rm -rf *.so", shell=True, executable="/bin/bash")

c_dir = 'c_source/'
c_sources = ['Python_extension.c', 'DynProgr_sse_short.c', 'DynProgr_sse_byte.c', 'DynProgr_scalar.c']

setup(
    ext_modules = cythonize([
	    Extension(
		"cython_swps3",
		sources=["cython_swps3.pyx"] + map(lambda c: c_dir + c, c_sources),
		#libraries=["swps3"],
		#extra_link_args=["-L../"], 
		include_dirs=[c_dir]
	)]
    ),
    name='cython_swps3',
    version='0.3',
    author='Your Name',
    author_email='your@email',
    #cmdclass = {'build_ext': build_ext},
    license='MIT',
    description='Example package that says hello'
)