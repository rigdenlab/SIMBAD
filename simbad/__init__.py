"""This is SIMBAD

Sequence Independent Molecular replacement Based on Available Database
"""

__authors__ = "Adam Simpkin & Felix Simkovic"
__contributing_authors__ = "Jens Thomas & Ronan Keegan"
__credits__ = "Daniel Rigden, William Shepard, Charles Ballard, Villi Uski, Andrey Lebedev"
__email__ = "hlasimpk@liv.ac.uk"

from distutils.version import StrictVersion
from simbad import version

__version__ = version.__version__

import os
import pyjob

if "CCP4" not in os.environ:
    raise RuntimeError("Cannot find CCP4 root directory")

if StrictVersion(pyjob.__version__) < StrictVersion("0.1.3"):
    raise RuntimeError("Please upgrade PyJOB")

SIMBAD_SHARE_STATIC_DIR = os.path.join(os.environ["CCP4"], "share", "simbad", "static")
LATTICE_DB = os.path.join(SIMBAD_SHARE_STATIC_DIR, "niggli_database.npz")
CONTAMINANT_MODELS = os.path.join(SIMBAD_SHARE_STATIC_DIR, "contaminants")
MORDA_MODELS = os.path.join(SIMBAD_SHARE_STATIC_DIR, "morda")
