import os

if "CCP4" not in os.environ.keys():
    msg = "Cannot find CCP4 root directory"
    raise RuntimeError(msg)

__all__ = ["SIMBAD_DIR", "SHARE_DIR", "SIMBAD_CONFIG_FILE"]

# SIMBAD source code directory
SIMBAD_DIR = os.path.join(os.environ["CCP4"], "lib", "py2", "simbad")

# SIMBAD share directory
SHARE_DIR = os.path.join(os.environ["CCP4"], "share", "simbad")

# SIMBAD configuration file
SIMBAD_CONFIG_FILE = os.path.join(SHARE_DIR, "include", "simbad.ini")