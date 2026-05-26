"""Build the trviz Cython extension.

All package metadata lives in pyproject.toml. This file exists only to
declare the C extension because setuptools' declarative pyproject support
for compiled extensions is still rough.

Build from .pyx when Cython is available, fall back to the committed .c
file otherwise.
"""

import numpy
from setuptools import Extension, setup

try:
    from Cython.Build import cythonize

    HAVE_CYTHON = True
except ImportError:
    HAVE_CYTHON = False

source = "trviz/cy/decompose.pyx" if HAVE_CYTHON else "trviz/cy/decompose.c"

extensions = [
    Extension(
        "trviz.cy.decompose",
        [source],
        extra_compile_args=["-O3"],
        include_dirs=[numpy.get_include()],
    )
]

if HAVE_CYTHON:
    extensions = cythonize(extensions)

setup(ext_modules=extensions)
