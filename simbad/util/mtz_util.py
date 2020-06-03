"""Module for MTZ file I/O and manipulation"""

__author__ = "Adam Simpkin & Jens Thomas & Felix Simkovic"
__date__ = "21 Jul 2019"
__version__ = "0.3"

import shutil

from mrbump.ccp4 import MRBUMP_ctruncate
from pyjob import cexec
from pyjob.script import EXE_EXT

from simbad.parsers.mtz_parser import MtzParser


def ctruncate(hklin, hklout):
    """Function to run Ctruncate on input MTZ to generate any missing columns"""

    ctr_colin = None
    ctr_colin_sig = None
    plus_minus = None

    mp = MtzParser(hklin)
    mp.parse()

    ctr = MRBUMP_ctruncate.Ctruncate()

    log_file = hklout.rsplit(".", 1)[0] + ".log"
    ctr.setlogfile(log_file)

    input_f = bool(mp.f)

    if mp.f and mp.sigf or mp.i and mp.sigi:
        plus_minus = False
        if mp.i:
            ctr_colin = mp.i
            ctr_colin_sig = mp.sigi
        else:
            ctr_colin = mp.f
            ctr_colin_sig = mp.sigf

    elif mp.iplus:
        plus_minus = True
        ctr_colin = []
        ctr_colin_sig = []
        ctr_colin.append(mp.iplus)
        ctr_colin.append(mp.iminus)
        ctr_colin_sig.append(mp.sigiplus)
        ctr_colin_sig.append(mp.sigiminus)

    elif mp.fplus:
        plus_minus = True
        ctr_colin = []
        ctr_colin_sig = []
        ctr_colin.append(mp.fplus)
        ctr_colin.append(mp.fminus)
        ctr_colin_sig.append(mp.sigfplus)
        ctr_colin_sig.append(mp.sigfminus)

    if mp.i and mp.sigi and mp.f and mp.sigf and mp.free:
        shutil.copyfile(hklin, hklout)
    elif mp.i and mp.free:
        ctr.ctruncate(
            hklin, hklout, ctr_colin, ctr_colin_sig, colout="from_SIMBAD", colinFREE=mp.free, USEINTEN=True,
            INPUTF=input_f, PLUSMINUS=plus_minus
        )
    elif mp.i and not mp.free:
        ctr.ctruncate(hklin, hklout, ctr_colin, ctr_colin_sig, colout="from_SIMBAD", USEINTEN=True, INPUTF=input_f,
                      PLUSMINUS=plus_minus)
    elif mp.free:
        ctr.ctruncate(hklin, hklout, ctr_colin, ctr_colin_sig, colout="from_SIMBAD", colinFREE=mp.free,
                      USEINTEN=False, PLUSMINUS=plus_minus)
    else:
        ctr.ctruncate(hklin, hklout, ctr_colin, ctr_colin_sig, colout="from_SIMBAD", USEINTEN=False,
                      PLUSMINUS=plus_minus)


def reindex(hklin, hklout, sg):
    """Function to reindex input hkl using pointless"""

    cmd = ["pointless" + EXE_EXT, "hklin", hklin, "hklout", hklout]
    stdin = """
spacegroup {0}
    """

    stdin = stdin.format(sg)
    cexec(cmd, stdin=stdin)
