import setuptools

with open('README.md', 'r', encoding='utf-8') as f:
    long_description = f.read()

setuptools.setup(
    name='trviz',
    use_scm_version=True,
    setup_requires=['setuptools_scm'],
    author='Jonghun Park',
    author_email='jop002@eng.ucsd.edu',
    description='A tool to visualize the tandem repeat polymorphisms',
    long_description=long_description,
    long_description_content_typ='text/markdown',
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: MIT License',
        'Operating System :: OS Independent',
    ],
    packages=setuptools.find_packages(),
    include_package_data=True,
    install_requires=[
        'matplotlib',
        'numpy',
        'biopython'
        'mafft'
    ],
    python_requires='>=3.8',
)