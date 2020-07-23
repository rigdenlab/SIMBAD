import os
import sys

from distutils.command.build import build
from distutils.util import convert_path
from setuptools import setup

import numpy as np


class BuildCommand(build):
    user_options = build.user_options + [
        ("script-python-path=", None, "Path to Python interpreter to be included in the scripts")
    ]

    def initialize_options(self):
        build.initialize_options(self)
        self.script_python_path = None

    def finalize_options(self):
        build.finalize_options(self)

    def run(self):
        global script_python_path
        script_python_path = self.script_python_path
        build.run(self)


def dependencies():
    with open("requirements.txt", "r") as f_in:
        deps = f_in.read().splitlines()
    return deps


def readme():
    with open("README.rst", "r") as f_in:
        return f_in.read()

def version():
    # Credits to http://stackoverflow.com/a/24517154
    main_ns = {}
    ver_path = convert_path("simbad/version.py")
    with open(ver_path) as f_in:
        exec(f_in.read(), main_ns)
    return main_ns["__version__"]


# ==============================================================
# Determine the Python executable
# ==============================================================
PYTHON_EXE = None
for arg in sys.argv:
    if arg[0:20] == "--script-python-path" and len(arg) == 20:
        option, value = arg, sys.argv[sys.argv.index(arg) + 1]
        PYTHON_EXE = value
    elif arg[0:20] == "--script-python-path" and arg[20] == "=":
        option, value = arg[:20], arg[21:]
        PYTHON_EXE = value

if not PYTHON_EXE:
    PYTHON_EXE = sys.executable

AUTHOR = "Adam Simpkin"
AUTHOR_EMAIL = "hlasimpk@liv.ac.uk"
DESCRIPTION = "Sequence independent MR pipeline"
DEPENDENCIES = dependencies()
LICENSE = "BSD License"
LONG_DESCRIPTION = readme()
PACKAGE_DIR = "simbad"
PACKAGE_NAME = "simbad"
PLATFORMS = ["Mac OS", "Windows", "Unix"]
URL = "http://www.simbad.rtfd.io/en/latest/"
VERSION = version()

PACKAGES = [
    "simbad",
    "simbad/command_line",
    "simbad/core",
    "simbad/db",
    "simbad/lattice",
    "simbad/mr",
    "simbad/parsers",
    "simbad/rotsearch",
    "simbad/util",
]

CLASSIFIERS = [
    "Development Status :: 4 - Beta",
    "Intended Audience :: Science/Research",
    "License :: OSI Approved :: BSD License",
    "Operating System :: OS Independent",
    "Programming Language :: Python",
    "Programming Language :: Python :: 3.6",
    "Programming Language :: Python :: 3.7",
    "Programming Language :: Python :: 3.8",
    "Topic :: Scientific/Engineering :: Bio-Informatics",
]

setup(
    cmdclass={"build": BuildCommand},
    author=AUTHOR,
    author_email=AUTHOR_EMAIL,
    name=PACKAGE_NAME,
    description=DESCRIPTION,
    include_dirs=[np.get_include()],
    long_description=LONG_DESCRIPTION,
    license=LICENSE,
    version=VERSION,
    url=URL,
    packages=PACKAGES,
    package_dir={PACKAGE_NAME: PACKAGE_DIR},
    platforms=PLATFORMS,
    classifiers=CLASSIFIERS,
    install_requires=DEPENDENCIES,
    include_package_data=True,
    zip_safe=False,
)

