from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize
import numpy
import subprocess


subprocess.Popen("rm -rf build", shell=True, executable="/bin/bash")
subprocess.Popen("rm -rf *.c", shell=True, executable="/bin/bash")
subprocess.Popen("rm -rf *.so", shell=True, executable="/bin/bash")


setup(
    ext_modules = cythonize([
	    Extension(
		"cython_swps3",
		["cython_swps3.pyx"],
		libraries=["swps3"],
		extra_link_args=["-L../"], 
		include_dirs=[numpy.get_include(), "../"]
	)]
    )
)
