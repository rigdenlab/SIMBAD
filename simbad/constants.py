import os

if "CCP4" not in os.environ.keys():
    msg = "Cannot find CCP4 root directory"
    raise RuntimeError(msg)

__all__ = ["SIMBAD_DIR", "SHARE_DIR", "SIMBAD_CONFIG_FILE", "SIMBAD_LATTICE_DB"]

# SIMBAD Egg root
SIMBAD_EGG_ROOT = os.path.dirname(__file__)

# SIMBAD configuration file
SIMBAD_CONFIG_FILE = os.path.join(SIMBAD_EGG_ROOT, 'static', 'simbad.ini')

# SIMBAD database for lattice search
SIMBAD_LATTICE_DB = os.path.join(SIMBAD_EGG_ROOT, 'static', 'niggli_database.cpk')