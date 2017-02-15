"""Sequence Independent Molecular Replacement Based on Available Database"""

from setuptools import setup, find_packages
from distutils.util import convert_path
import glob
import os
import sys

# We need this in case somebody wants to install the data elsewhere.
# If done, the constants.py variable will not match anymore
if '--install-data' in sys.argv:
    raise Exception("Cannot yet handle the `--install-data` flag")


def get_files(path):
    """Find any files to be included"""
    for p in path:
        for f in glob.glob(os.path.join(p, '*')):
            if os.path.isfile(f):
                yield f


def get_version():
    """Get the current SIMBAD version"""
    # Credits to http://stackoverflow.com/a/24517154
    main_ns = {}
    ver_path = convert_path(os.path.join('simbad', 'util', 'version.py'))
    with open(ver_path) as f_in:
        exec(f_in.read(), main_ns)
    return main_ns['__version__']globals()

# Obtain the current version of SIMBAD
__version__ = get_version()


# Do the actual setup below
setup(
    name='simbad',
    description=__doc__.replace("\n", ""),
    long_description=open('README.md').read(),
    version=__version__,
    author='Adam Simpkin',
    author_email='hlasimpk@liverpool.ac.uk',
    license='BSD License',
    url='https://github.com/rigdenlab/SIMBAD',
    download_url='https://github.com/rigdenlab/SIMBAD/tarball/' + __version__,
    package_dir={'simbad': 'simbad'},
    packages=find_packages(exclude="tests"),
    data_files=[
        ('static', list(get_files(['static']))),
        ('static/contaminant_models', list(get_files([os.path.join('static', 'contaminant_models')]))),
    ],
    scripts=list(get_files(['scripts'])),
    platforms=['Linux', 'Mac OS-X', 'Unix'],
    install_requires=['numpy>=1.8.2', 'biopython>=1.64'],
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: BSD License",
        "Operating System :: OS Independent",
        "Programming Language :: Python",
        "Programming Language :: Python :: 2.7",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
    include_package_data=True,
    zip_safe=False,
)

