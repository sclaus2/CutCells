import os
import shlex
import subprocess
import sys
import sysconfig

from setuptools import Extension, setup
from setuptools.command.build_ext import build_ext

if sys.version_info < (3, 8):
    print("Python 3.8 or higher required, please upgrade.")
    sys.exit(1)

VERSION = "0.1"

REQUIREMENTS = [ "numpy>=1.21" ]


class CMakeExtension(Extension):
    def __init__(self, name, sourcedir=''):
        Extension.__init__(self, name, sources=[])
        self.sourcedir = os.path.abspath(sourcedir)


class CMakeBuild(build_ext):
    def run(self):
        try:
            _ = subprocess.check_output(['cmake', '--version'])
        except OSError:
            raise RuntimeError("CMake must be installed to build the following extensions: "
                               + ", ".join(e.name for e in self.extensions))

        for ext in self.extensions:
            self.build_extension(ext)

    def build_extension(self, ext):
        extdir = os.path.abspath(os.path.dirname(self.get_ext_fullpath(ext.name)))
        cmake_args = shlex.split(os.environ.get("CMAKE_ARGS", ""))
        cmake_args += ['-DCMAKE_LIBRARY_OUTPUT_DIRECTORY=' + extdir,
                       '-DPython3_EXECUTABLE=' + sys.executable,
                       f'-DPython3_LIBRARIES={sysconfig.get_config_var("LIBDEST")}',
                       f'-DPython3_INCLUDE_DIRS={sysconfig.get_config_var("INCLUDEPY")}']

        cfg = 'Debug' if self.debug else 'Release'
        build_args = ['--config', cfg]
        cmake_args += ['-DCMAKE_BUILD_TYPE=' + cfg]

        env = os.environ.copy()
        # default to 3 build threads
        if "CMAKE_BUILD_PARALLEL_LEVEL" not in env:
            env["CMAKE_BUILD_PARALLEL_LEVEL"] = "3"

        import pybind11
        env['pybind11_DIR'] = pybind11.get_cmake_dir()
        env['CXXFLAGS'] = '{} -DVERSION_INFO=\\"{}\\"'.format(env.get('CXXFLAGS', ''),
                                                              self.distribution.get_version())

        if not os.path.exists(self.build_temp):
            os.makedirs(self.build_temp)
        subprocess.check_call(['cmake', ext.sourcedir] + cmake_args, cwd=self.build_temp, env=env)
        subprocess.check_call(['cmake', '--build', '.'] + build_args, cwd=self.build_temp, env=env)


setup(name='cutcells',
      version=VERSION,
      author='Susanne Claus (ONERA)',
      description='CutCells Python interface',
      long_description='',
      packages=["cutcells"],
      package_data={'cutcells': ['*.h']},
      ext_modules=[CMakeExtension('cutcells._cutcellscpp')],
      cmdclass=dict(build_ext=CMakeBuild),
      install_requires=REQUIREMENTS,
      setup_requires=["pybind11"],
      zip_safe=False)
