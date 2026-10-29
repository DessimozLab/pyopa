from setuptools import setup
from setuptools.extension import Extension
from setuptools.command.build_ext import build_ext
from setuptools.errors import CompileError, LinkError
from Cython.Build import cythonize
import numpy
import os
import shutil
import subprocess
import sys
import tempfile

# project metadata lives in pyproject.toml; this file only defines the C extension
package_name = 'pyopa'

c_dir = 'c_source/'
c_sources = ['Python_extension.c', 'DynProgr_sse_short.c', 'DynProgr_sse_byte.c',
             'DynProgr_scalar.c', 'EstimatePam.c', 'EstimatePamEigen.c', 'Page_size.c', 'DynProgr_sse_double.c']


def openmp_flags(compiler_type):
    """returns the (compile, link) flags to enable OpenMP for the given compiler"""
    if compiler_type == 'msvc':
        return ['/openmp'], []
    if sys.platform == 'darwin':
        # Apple clang has no OpenMP runtime; use libomp from OPENMP_PREFIX
        # (set for the wheel builds) or from Homebrew
        prefix = os.environ.get('OPENMP_PREFIX')
        if not prefix and shutil.which('brew'):
            prefix = subprocess.run(['brew', '--prefix', 'libomp'], capture_output=True,
                                    text=True).stdout.strip()
        if not prefix:
            return None
        return (['-Xpreprocessor', '-fopenmp', '-I' + os.path.join(prefix, 'include')],
                ['-L' + os.path.join(prefix, 'lib'), '-lomp', '-Wl,-rpath,' + os.path.join(prefix, 'lib')])
    return ['-fopenmp'], ['-fopenmp']


class build_ext_openmp(build_ext):
    """Enables OpenMP if the compiler supports it. Otherwise the extension is built
    without it and the parallel loops (e.g. estimate_many_pam) run serially, unless
    PYOPA_REQUIRE_OPENMP is set, in which case the build fails."""

    def build_extensions(self):
        flags = openmp_flags(self.compiler.compiler_type)
        if flags is not None and self._compiles_with_openmp(*flags):
            for ext in self.extensions:
                ext.extra_compile_args += flags[0]
                ext.extra_link_args += flags[1]
        elif os.environ.get('PYOPA_REQUIRE_OPENMP'):
            raise RuntimeError('OpenMP is required (PYOPA_REQUIRE_OPENMP is set), but not available')
        else:
            print('WARNING: OpenMP not available, building pyopa without parallel support')
        super().build_extensions()

    def _compiles_with_openmp(self, compile_args, link_args):
        with tempfile.TemporaryDirectory() as tmpdir:
            src = os.path.join(tmpdir, 'test_openmp.c')
            with open(src, 'w') as f:
                f.write('#include <omp.h>\nint main(void) { return omp_get_max_threads() > 0 ? 0 : 1; }\n')
            try:
                objs = self.compiler.compile([src], output_dir=tmpdir, extra_postargs=compile_args)
                self.compiler.link_executable(objs, os.path.join(tmpdir, 'test_openmp'),
                                              extra_postargs=link_args)
            except (CompileError, LinkError):
                return False
        return True


setup(
    ext_modules=cythonize([
        Extension(
            package_name + ".backend.pyopa",
            sources=[os.path.join(package_name, 'backend', 'pyopa.pyx')] +
                    list(map(lambda c: os.path.join(c_dir, c), c_sources)),
            extra_compile_args=['-O2', '-DSIMDE_ENABLE_NATIVE_ALIASES'], #, '-DPY_DEBUG'],
            include_dirs=[c_dir, os.path.join(package_name, "backend"), numpy.get_include(), '.']
        )],
        compiler_directives={"language_level": "3"},
    ),
    cmdclass={'build_ext': build_ext_openmp},
    include_package_data=True,
)
