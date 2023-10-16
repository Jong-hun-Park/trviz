from setuptools import Extension, setup, find_packages
from Cython.Build import cythonize
import numpy

with open('README.md', 'r', encoding='utf-8') as f:
    long_description = f.read()

extensions = [
    Extension(
        "trviz.cy.decompose",
        ["trviz/cy/decompose.pyx"],
        extra_compile_args=["-O3", '-march=native'],
        include_dirs=[numpy.get_include()],
    )
]

setup(
    name='trviz',
    version="1.0.1",
    author='Jonghun Park',
    author_email='jop002@ucsd.edu',
    description='A python library for decomposing and visualizing tandem repeat sequences',
    url="https://github.com/Jong-hun-Park/trviz",
    long_description=long_description,
    long_description_content_type='text/markdown',
    classifiers=[
        'Programming Language :: Python :: 3.8',
        'License :: OSI Approved :: BSD License',
        'Operating System :: OS Independent',
    ],
    packages=find_packages(),
    include_package_data=True,
    install_requires=[
        'matplotlib',
        'numpy',
        'biopython',
    ],
    python_requires='>=3.8',
    ext_modules=cythonize(extensions),
)
