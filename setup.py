from setuptools import Extension, setup, find_packages
import numpy
import sys
import platform

try:
    from Cython.Build import cythonize
except ImportError:
    USE_CYTHON = False
else:
    USE_CYTHON = True

ext = 'pyx' if USE_CYTHON else 'c'

with open('README.md', 'r', encoding='utf-8') as f:
    long_description = f.read()

use_cython = '--cython' in sys.argv
if '--cython' in sys.argv:
    print("Use cython")
    sys.argv.remove('--cython')
    USE_CYTHON = True
else:
    USE_CYTHON = False

extensions = [
    Extension(
        "trviz.cy.decompose",
        ["trviz/cy/decompose.{}".format(ext)],
        extra_compile_args=["-O3"],
        include_dirs=[numpy.get_include()],
    )
]

if USE_CYTHON:
    try:
        extensions = cythonize(extensions)
    except Exception as e:
        print(e)
        print("Cythonizing failed. Continue without cythonizing")

setup(
    name='trviz',
    version="1.3.0",
    author='Jonghun Park',
    author_email='jop002@ucsd.edu',
    description='A python library for decomposing and visualizing tandem repeat sequences',
    url="https://github.com/Jong-hun-Park/trviz",
    long_description=long_description,
    long_description_content_type='text/markdown',
    classifiers=[
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3.10',
        'Programming Language :: Python :: 3.11',
        'Programming Language :: Python :: 3.12',
        'License :: OSI Approved :: BSD License',
        'Operating System :: OS Independent',
    ],
    packages=find_packages(),
    include_package_data=True,
    install_requires=[
        'matplotlib',
        'numpy',
        'biopython',
        'scipy',
        'distinctipy',
    ],
    python_requires='>=3.8',
    ext_modules=extensions if not use_cython else [],
)
