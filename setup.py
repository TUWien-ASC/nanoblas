from skbuild import setup

import sys
_cmake_args = []

# if 'win32' in sys.platform:
#    _cmake_args += ['-G', 'MinGW Makefiles']

    
setup(
    name="nanoblas",
    version="1.0.0",
    author="J. Schoeberl",
    license="MIT",
    packages=["nanoblas"],
    cmake_args=_cmake_args
)
