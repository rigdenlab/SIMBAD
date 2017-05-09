import os

if "CCP4" not in os.environ.keys():
    msg = "Cannot find CCP4 root directory"
    raise RuntimeError(msg)

# SIMBAD egg directory
SIMBAD_EGG_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))

# SIMBAD database for lattice search
SIMBAD_LATTICE_DB = os.path.join(SIMBAD_EGG_ROOT, 'static', 'niggli_database.npz')

# SIMBAD database of contaminant models
CONTAMINANT_MODELS = os.path.join(SIMBAD_EGG_ROOT, 'static', 'contaminants')
