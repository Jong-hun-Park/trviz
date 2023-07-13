import setuptools

with open('README.md', 'r', encoding='utf-8') as f:
    long_description = f.read()

setuptools.setup(
    name='trviz',
    version="1.0.1",
    author='Jonghun Park',
    author_email='jop002@eng.ucsd.edu',
    description='A python library for decomposing and visualizing tandem repeat sequences',
    url="https://github.com/Jong-hun-Park/trviz",
    long_description=long_description,
    long_description_content_type='text/markdown',
    classifiers=[
        'Programming Language :: Python :: 3.8',
        'License :: OSI Approved :: BSD License',
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
