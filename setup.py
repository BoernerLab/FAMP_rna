from setuptools import setup, find_packages

setup(
    name='famp',
    version='0.0.1',
    packages=find_packages(where="src"),
    package_dir={'': 'src'},
    author="Felix Erichson",
    license="GPL-3.0 license",
    classifiers=["Programming Language :: Python :: 3"],
    scripts=["src/famp/scripts/mac_os/rosetta/mac_extract_pdb.sh","src/famp/scripts/linux/rosetta/linux_extract_pdb.sh", "src/famp/scripts/dye_properties.csv"],
    package_data={'': ['src/famp/scripts/gromacs/mdp']},
    include_package_data=True,
)