"""Sequence Independent Molecular Replacement Based on Available Database"""

from distutils.command.build import build
from distutils.util import convert_path
from setuptools import setup

import distutils 
print(distutils.__file__)

import glob
import os
import sys

# ==============================================================
# Setup.py command extensions
# ==============================================================

# Credits to http://stackoverflow.com/a/33181352
class BuildCommand(build):
    user_options = build.user_options + [ 
        ('script-python-path=', None, 'Path to Python interpreter to be included in the scripts')
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

# ==============================================================
# Functions, functions, functions ... 
# ==============================================================

def dependencies():
    with open('requirements.txt', 'r') as f_in:
        return [l for l in f_in.read().rsplit(os.linesep) 
                if l and not l.startswith("#")]

def files(path):
    """Find any files to be included"""
    fs = []
    for p in path:
        for f in glob.glob(os.path.join(p, '*')):
            if os.path.isfile(f):
                fs.append(f)
    return fs

def readme():
    with open('README.md', 'r') as f_in:
        return f_in.read()


def scripts():
    extension = '.bat' if sys.platform.startswith('win') else ''
    header = '' if sys.platform.startswith('win') else '#!/bin/sh'
    bin_dir = 'bin'
    command_dir = convert_path('simbad/command_line')
    scripts = []
    for file in os.listdir(command_dir):
        if not file.startswith('_') and file.endswith('.py'):
            # Make sure we have a workable name
            f_name = os.path.basename(file).rsplit('.', 1)[0]
            for c in ['.', '_']:
                new_f_name = f_name.replace(c, '-')
            # Special case for simbad-main
            if new_f_name == "simbad-main":
                new_f_name = "simbad"
            # Write the content of the script
            script = os.path.join(bin_dir, new_f_name + extension)
            with open(script, "w") as f_out:
                f_out.write(header + os.linesep)
                # BATCH file
                if sys.platform.startswith('win'):
                    string = "@{0} -m simbad.command_line.{1} %*"
                # BASH file
                else:
                    string = "{0} -m simbad.command_line.{1} \"$@\""
                f_out.write(string.format(PYTHON_EXE, f_name) + os.linesep)
            os.chmod(script, 0o777)
            scripts.append(script)
    return scripts


def version():
    """Get the current SIMBAD version"""
    # Credits to http://stackoverflow.com/a/24517154
    main_ns = {}
    ver_path = convert_path(os.path.join('simbad', 'util', 'version.py'))
    with open(ver_path) as f_in:
        exec(f_in.read(), main_ns)
    return main_ns['__version__']


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

# ==============================================================
# Define all the relevant options
# ==============================================================
AUTHOR = "Adam Simpkin"
AUTHOR_EMAIL = "hlasimpk@liverpool.ac.uk"
DESCRIPTION = __doc__.replace("\n", "")
DEPENDENCIES = dependencies()
LICENSE = "BSD License"
LONG_DESCRIPTION = readme()
PACKAGE_DIR = "simbad"
PACKAGE_NAME = "simbad"
SCRIPTS = scripts()
URL = "https://github.com/rigdenlab/SIMBAD"
VERSION = version()

PACKAGES = [
    'simbad',
    'simbad/command_line',
    'simbad/lattice',
    'simbad/parsers',
    'simbad/util',
]

DATA_FILES = [
    ('static', files(['static'])),
    ('static/contaminant_models', files([os.path.join('static', 'contaminant_models')])),
]

CLASSIFIERS = [
    "Development Status :: 4 - Beta",
    "Intended Audience :: Science/Research",
    "License :: OSI Approved :: BSD License",
    "Programming Language :: Python",
    "Programming Language :: Python :: 2.7",
    "Topic :: Scientific/Engineering :: Bio-Informatics",
]

# Do the actual setup below
setup(
    cmdclass={
        'build': BuildCommand,
    },
    author=AUTHOR,
    author_email=AUTHOR_EMAIL,
    name=PACKAGE_NAME,
    description=DESCRIPTION,
    long_description=LONG_DESCRIPTION,
    license=LICENSE,
    version=VERSION,
    url=URL,
    packages=PACKAGES,
    package_dir={PACKAGE_NAME: PACKAGE_DIR},
    scripts=SCRIPTS,
    install_requires=DEPENDENCIES,
    data_files=DATA_FILES,
    classifiers=CLASSIFIERS,
    test_suite='nose.collector',
    tests_require=['nose >=1.3.7'],
    include_package_data=True,
    zip_safe=False,
)

