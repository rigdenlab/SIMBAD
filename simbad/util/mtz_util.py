"""Module for MTZ file I/O and manipulation"""

__author__ = "Adam Simpkin & Jens Thomas & Felix Simkovic"
__date__ = "21 Jul 2019"
__version__ = "0.3"

import shutil
import warnings

from pyjob import cexec
from pyjob.script import EXE_EXT

from simbad.parsers.mtz_parser import MtzParser


def deprecate(version, msg=None):
    """Decorator to deprecate Python classes and functions

    Parameters
    ----------
    version : str
       A string containing the version with which the callable is removed
    msg : str, optional
       An additional message that will be displayed alongside the default message
    """

    def deprecate_decorator(callable_):
        def warn(*args, **kwargs):
            message = "%s has been deprecated and will be removed in version %s!" % (callable_.__name__, version)
            if msg:
                message += " - %s" % msg
            warnings.warn(message, DeprecationWarning)
            return callable_(*args, **kwargs)

        return warn

    return deprecate_decorator


def ctruncate(hklin, hklout):
    """Function to run Ctruncate on input MTZ to generate any missing columns"""
    from mrbump.ccp4 import MRBUMP_ctruncate

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

    elif mp.i_plus:
        plus_minus = True
        ctr_colin = []
        ctr_colin_sig = []
        ctr_colin.append(mp.i_plus)
        ctr_colin.append(mp.i_minus)
        ctr_colin_sig.append(mp.sigi_plus)
        ctr_colin_sig.append(mp.sigi_minus)

    elif mp.f_plus:
        plus_minus = True
        ctr_colin = []
        ctr_colin_sig = []
        ctr_colin.append(mp.f_plus)
        ctr_colin.append(mp.f_minus)
        ctr_colin_sig.append(mp.sigf_plus)
        ctr_colin_sig.append(mp.sigf_minus)

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


@deprecate('0.2.1', msg="Use simbad.parsers.mtz_parser.MtzParser instead")
class GetLabels(object):
    def __init__(self, hklin):
        self.f = None
        self.sigf = None
        self.i = None
        self.sigi = None
        self.fplus = None
        self.sigfplus = None
        self.fminus = None
        self.sigfminus = None
        self.iplus = None
        self.sigiplus = None
        self.iminus = None
        self.sigiminus = None
        self.dano = None
        self.sigdano = None
        self.free = None

        self.run(hklin)

    def run(self, hklin):
        mp = MtzParser(hklin)
        mp.parse()

        self.f = mp.f
        self.sigf = mp.sigf
        self.i = mp.i
        self.sigi = mp.sigi
        self.fplus = mp.f_plus
        self.sigfplus = mp.sigf_plus
        self.fminus = mp.f_minus
        self.sigfminus = mp.sigf_plus
        self.iplus = mp.i_plus
        self.sigiplus = mp.sigi_plus
        self.iminus = mp.i_minus
        self.sigiminus = mp.sigi_minus
        self.dano = mp.dp
        self.sigdano = mp.sigdp
        self.free = mp.free


@deprecate('0.2.1', msg="Use simbad.parsers.mtz_parser.MtzParser instead")
def crystal_data(hklin):
    mp = MtzParser(hklin)
    space_group = "".join(mp.spacegroup_symbol.encode("ascii").split())
    resolution = mp.resolution
    cell_parameters = (mp.cell.a,
                       mp.cell.b,
                       mp.cell.c,
                       mp.cell.alpha,
                       mp.cell.beta,
                       mp.cell.gamma)
    return space_group, resolution, cell_parameters





