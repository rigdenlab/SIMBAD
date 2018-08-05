"""Module for MTZ file I/O and manipulation"""

__author__ = "Adam Simpkin & Jens Thomas"
__date__ = "17 May 2017"
__version__ = "0.2"

import logging
import numpy
import os
import sys

sys.path.append(os.path.join(os.environ["CCP4"], "share", "mrbump", "include", "ccp4"))
import MRBUMP_ctruncate

from iotbx import reflection_file_reader
from pyjob import cexec

logger = logging.getLogger(__name__)


def ctruncate(hklin, hklout):
    """Function to run Ctruncate on input MTZ to generate any missing columns"""

    ctr_colin = None
    ctr_colin_sig = None
    plus_minus = None

    mtz_obj = GetLabels(hklin)

    ctr = MRBUMP_ctruncate.Ctruncate()

    log_file = hklout.rsplit(".", 1)[0] + '.log'
    ctr.setlogfile(log_file)

    if mtz_obj.f:
        input_f = True
    else:
        input_f = False

    if mtz_obj.f or mtz_obj.i:
        plus_minus = False
        if mtz_obj.i:
            ctr_colin = mtz_obj.i
            ctr_colin_sig = mtz_obj.sigi
        else:
            ctr_colin = mtz_obj.f
            ctr_colin_sig = mtz_obj.sigf

    elif mtz_obj.iplus:
        plus_minus = True
        ctr_colin = []
        ctr_colin_sig = []
        ctr_colin.append(mtz_obj.fplus)
        ctr_colin.append(mtz_obj.fminus)
        ctr_colin_sig.append(mtz_obj.sigfplus)
        ctr_colin_sig.append(mtz_obj.sigfminus)

    elif mtz_obj.fplus:
        plus_minus = True
        ctr_colin = []
        ctr_colin_sig = []
        ctr_colin.append(mtz_obj.fplus)
        ctr_colin.append(mtz_obj.fminus)
        ctr_colin_sig.append(mtz_obj.sigfplus)
        ctr_colin_sig.append(mtz_obj.sigfminus)

    if mtz_obj.i and mtz_obj.free:
        ctr.ctruncate(hklin, hklout, ctr_colin, ctr_colin_sig, colout="from_SIMBAD", colinFREE=mtz_obj.free,
                      USEINTEN=True, INPUTF=input_f, PLUSMINUS=plus_minus)
    elif mtz_obj.i and not mtz_obj.free:
        ctr.ctruncate(hklin, hklout, ctr_colin, ctr_colin_sig, colout="from_SIMBAD", USEINTEN=True, INPUTF=input_f,
                      PLUSMINUS=plus_minus)
    elif mtz_obj.free:
        ctr.ctruncate(hklin, hklout, ctr_colin, ctr_colin_sig, colout="from_SIMBAD", colinFREE=mtz_obj.free,
                      USEINTEN=False, PLUSMINUS=plus_minus)
    else:
        ctr.ctruncate(hklin, hklout, ctr_colin, ctr_colin_sig, colout="from_SIMBAD", USEINTEN=False,
                      PLUSMINUS=plus_minus)


def reindex(hklin, hklout, sg):
    """Function to reindex input hkl using pointless"""

    cmd = ["pointless", "hklin", hklin, "hklout", hklout]

    stdin = """
reindex k,h,-l
spacegroup {0}
    """

    stdin = stdin.format(sg)
    cexec(cmd, stdin=stdin)


def crystal_data(mtz_file):
    """Set crystallographic parameters from mtz file

    Parameters
    ----------
    mtz_file : str
       The path to the mtz file

    Returns
    -------
    space_group : str
       The space group
    resolution : str
       The resolution
    cell_parameters : tuple
       The cell parameters

    """

    reflection_file = reflection_file_reader.any_reflection_file(file_name=mtz_file)
    content = reflection_file.file_content()
    space_group = content.space_group_name().replace(" ", "")
    resolution = content.max_min_resolution()[1]
    cell_parameters = content.crystals()[0].unit_cell_parameters()

    return space_group, resolution, cell_parameters


class GetLabels(object):
    """Class to get the column labels for input mtz file

    Attributes
    ----------
    mtz_file : str
       The path to the mtz file

    Returns
    -------
    f : str
        f column label
    sigf : str
        fp column label
    i : str
        i column label
    sigi : str
        sigi column label
    fplus : str
        f(+) column label
    sigfplus : str
        sigf(+) column label
    fminus : str
        f(-) column label
    sigfminus : str
        sigf(-) column label
    iplus : str
        i(+) column label
    sigiplus : str
        sigi{+} column label
    iminus : str
        i(-) column label
    sigiminus : str
        sigi(-) column label
    dano : str
        dano column label
    sigdano : str
        sigdano column label
    free : str
        free column label
    """

    def __init__(self, mtz_file):
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

        self.run(mtz_file)

    def run(self, mtz_file):
        reflection_file = reflection_file_reader.any_reflection_file(file_name=mtz_file)
        if not reflection_file.file_type() == "ccp4_mtz":
            msg="File is not of type ccp4_mtz: {0}".format(mtz_file)
            logging.critical(msg)
            raise RuntimeError(msg)

        content = reflection_file.file_content()
        ctypes = content.column_types()
        npctypes = numpy.array(ctypes)
        clabels = content.column_labels()

        dtype = 'D'
        ftype = 'F'
        ktype = 'K'
        gtype = 'G'
        jtype = 'J'
        ltype = 'L'
        mtype = 'M'

        if ftype in ctypes:
            self.f = clabels[ctypes.index(ftype)]
            self.sigf = 'SIG' + self.f
            fp_alt = 'PHI' + self.f
            if self.sigf in clabels:
                pass
            elif fp_alt in clabels:
                self.sigf = fp_alt
                pass
            else:
                msg = "Cannot find label {0} or {1} in file: {2}".format(self.sigf, fp_alt, mtz_file)
                logging.warning(msg)

        if jtype in ctypes:
            self.i = clabels[ctypes.index(jtype)]
            self.sigi = 'SIG' + self.i
            if self.sigi not in clabels:
                msg = "Cannot find label {0} in file: {1}".format(self.sigi, mtz_file)
                raise RuntimeError(msg)

        if ktype in ctypes:
            indices = numpy.where(npctypes == ktype)[0]
            for index in indices:
                if '+' in clabels[index]:
                    self.iplus = clabels[index]
                elif '-' in clabels[index]:
                    self.iminus = clabels[index]

        if mtype in ctypes:
            indices = numpy.where(npctypes == mtype)[0]
            for index in indices:
                if '+' in clabels[index]:
                    self.sigiplus = clabels[index]
                elif '-' in clabels[index]:
                    self.sigiminus = clabels[index]

        if gtype in ctypes:
            indices = numpy.where(npctypes == gtype)[0]
            for index in indices:
                if '+' in clabels[index]:
                    self.fplus = clabels[index]
                elif '-' in clabels[index]:
                    self.fminus = clabels[index]

        if ltype in ctypes:
            indices = numpy.where(npctypes == ltype)[0]
            for index in indices:
                if '+' in clabels[index]:
                    self.sigfplus = clabels[index]
                elif '-' in clabels[index]:
                    self.sigfminus = clabels[index]

        if dtype in ctypes:
            self.dano = clabels[ctypes.index(dtype)]
            self.sigdano = 'SIG' + self.dano
            sigdano_alt = 'SD' + self.dano
            if self.sigdano in clabels:
                pass
            elif sigdano_alt in clabels:
                self.sigdano = sigdano_alt
                pass
            else:
                msg = "Cannot find label {0} in file: {1}".format(self.sigdano, mtz_file)
                raise RuntimeError(msg)

        for label in clabels:
            for word in ["free", "test", "cross", "status", "flag"]:
                if label.lower().find(word) >= 0:
                    if self.free:
                        logger.warning("FOUND >1 R FREE label in file!")
                    self.free = label
                    break


