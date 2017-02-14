"""Sequence Independent Molecular replacement Based on Available Database"""

from setuptools import setup, find_packages
from distutils.util import convert_path
import glob

def get_version():
    # Credits to http://stackoverflow.com/a/24517154
    main_ns = {}
    ver_path = convert_path('simbad/util/version.py')
    with open(ver_path) as f_in:
        exec(f_in.read(), main_ns)
    return main_ns['__version__']

# Obtain the current version of ConKit
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
    scripts=[script for script in glob.glob("scripts/*")],
    platforms=['Linux', 'Mac OS-X', 'Unix'],
    # Any non-standard dependencies here
    install_requires=['numpy', 'biopython'],
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








