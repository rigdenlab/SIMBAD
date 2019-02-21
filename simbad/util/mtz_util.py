"""Module for MTZ file I/O and manipulation"""

__author__ = "Adam Simpkin & Jens Thomas"
__date__ = "17 May 2017"
__version__ = "0.2"

import logging
import os
import shutil
import sys

sys.path.append(os.path.join(os.environ["CCP4"], "share", "mrbump", "include", "ccp4"))
import MRBUMP_ctruncate

from iotbx import reflection_file_reader
from iotbx.reflection_file_utils import looks_like_r_free_flags_info
from pyjob import cexec
from pyjob.script import EXE_EXT

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

    if mtz_obj.i and mtz_obj.f and mtz_obj.free:
        shutil.copyfile(hklin, hklout)
    elif mtz_obj.i and mtz_obj.free:
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

    cmd = ["pointless" + EXE_EXT, "hklin", hklin, "hklout", hklout]
    stdin = """
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

        miller_arrays = reflection_file.as_miller_arrays()

        for m_a in miller_arrays:
            if looks_like_r_free_flags_info(m_a.info()) and not self.free:
                self.free = m_a.info().labels[0]
            elif self.check_anomalous(m_a):
                if self.check_for_dano_labels(m_a):
                    if len(m_a.info().labels) == 5:
                        self.f, self.sigf, self.dano, self.sigdano, isym = m_a.info().labels
                    elif len(m_a.info().labels) == 4:
                        self.f, self.sigf, self.dano, self.sigdano = m_a.info().labels
                    elif len(m_a.info().labels) == 2:
                        self.dano, self.sigdano = m_a.info().labels
                    else:
                        msg = "Unexpected number of columns found in anomalous miller array"
                        logging.critical(msg)
                elif self.check_for_plus_minus_labels(m_a):
                    if m_a.is_xray_amplitude_array():
                        self.fplus, self.sigfplus, self.fminus, self.sigfminus = m_a.info().labels
                    elif m_a.is_xray_intensity_array():
                        self.iplus, self.sigiplus, self.iminus, self.sigiminus = m_a.info().labels
                    else:
                        msg = "Type of anomalous miller array unknown"
                        logging.critical(msg)
                else:
                    msg = "Type of anomalous miller array unknown"
                    logging.critical(msg)
            elif m_a.is_xray_intensity_array() and len(m_a.info().labels) == 2 and not self.i:
                self.i, self.sigi = m_a.info().labels
            elif m_a.is_xray_amplitude_array() and len(m_a.info().labels) == 2 and not self.f:
                self.f, self.sigf = m_a.info().labels
            else:
                pass

    def check_anomalous(self, miller_array):
        if miller_array.anomalous_flag():
            return True
        elif miller_array.info().type_hints_from_file == 'anomalous_difference':
            return True
        # Check for anomalous miller arrays which aren't properly labeled
        elif self.check_for_dano_labels(miller_array):
            return True
        elif self.check_for_plus_minus_labels(miller_array):
            return True
        return False

    @staticmethod
    def check_for_dano_labels(miller_array):
        if any(['DANO' in i.upper() for i in miller_array.info().labels]):
            return True
        return False

    @staticmethod
    def check_for_plus_minus_labels(miller_array):
        if any(['(+)' in i for i in miller_array.info().labels]):
            return True
        return False



