"""This is SIMBAD

Sequence Independent Molecular replacement Based on Available Database
"""

__authors__ = "Adam Simpkin & Felix Simkovic"
__contributing_authors__ = "Jens Thomas & Ronan Keegan"
__credits__ = "Daniel Rigden, William Shepard, Charles Ballard, Villi Uski, Andrey Lebedev"
__email__ = "hlasimpk@liv.ac.uk"
from simbad import version
__version__ = version.__version__

import os

if "CCP4" not in os.environ:
    msg = "Cannot find CCP4 root directory"
    raise RuntimeError(msg)

# SIMBAD egg directory
EGG_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))

# SIMBAD database for lattice search
LATTICE_DB = os.path.join(EGG_ROOT, 'static', 'niggli_database.npz')

# SIMBAD database of contaminant models
CONTAMINANT_MODELS = os.path.join(EGG_ROOT, 'static', 'contaminants')
