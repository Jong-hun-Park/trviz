import setuptools

with open('README.md', 'r', encoding='utf-8') as f:
    long_description = f.read()

setuptools.setup(
    name='trviz',
    use_scm_version=True,
    setup_requires=['setuptools_scm'],
    author='Jonghun Park',
    author_email='jop002@eng.ucsd.edu',
    description='A python library for decomposing and visualizing tandem repeat sequences',
    long_description=long_description,
    long_description_content_typ='text/markdown',
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: BSD-3-Clause',
        'Operating System :: OS Independent',
    ],
    packages=setuptools.find_packages(),
    include_package_data=True,
    install_requires=[
        'matplotlib',
        'numpy',
        'biopython',
    ],
    python_requires='>=3.8',
)
