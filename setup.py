from setuptools import setup, find_packages

setup(
    name='famp',
    version='0.0.1',
    packages=find_packages(where="src"),
    package_dir={'': 'src'},
    author="Felix Erichson",
    license="GPL-3.0 license",
    classifiers=["Programming Language :: Python :: 3"],
    scripts=["src/famp/scripts/mac_os/rosetta/extract_pdb.sh"],
    package_data={'': ['src/famp/scripts/gromacs']},
    include_package_data=True,
)