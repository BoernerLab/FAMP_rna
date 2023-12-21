from setuptools import setup, find_packages

setup(
    name='famp',
    version='0.0.1',
    packages=find_packages("src/famp"),
    package_dir={'': 'src'},
    author="Felix Erichson",
    license="MIT License",
    classifiers=["Programming Language :: Python ::"],
)