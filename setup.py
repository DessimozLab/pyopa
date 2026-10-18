from setuptools import setup
from setuptools.extension import Extension
from Cython.Build import cythonize
import numpy
import os

# project metadata lives in pyproject.toml; this file only defines the C extension
package_name = 'pyopa'

c_dir = 'c_source/'
c_sources = ['Python_extension.c', 'DynProgr_sse_short.c', 'DynProgr_sse_byte.c',
             'DynProgr_scalar.c', 'EstimatePam.c', 'Page_size.c', 'DynProgr_sse_double.c']

setup(
    ext_modules=cythonize([
        Extension(
            package_name + ".backend.pyopa",
            sources=[os.path.join(package_name, 'backend', 'pyopa.pyx')] +
                    list(map(lambda c: os.path.join(c_dir, c), c_sources)),
            extra_compile_args=['-O2', '-DSIMDE_ENABLE_NATIVE_ALIASES', '-fopenmp', '-march=native'], #, '-DPY_DEBUG'],
            extra_link_args=['-fopenmp'],
            include_dirs=[c_dir, os.path.join(package_name, "backend"), numpy.get_include(), '.']
        )],
        compiler_directives={"language_level": "3"},
    ),
    include_package_data=True,
)
