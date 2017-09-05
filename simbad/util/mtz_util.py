"""Module for MTZ file I/O and manipulation"""

__author__ = "Adam Simpkin & Jens Thomas"
__date__ = "17 May 2017"
__version__ = "0.2"

from iotbx import reflection_file_reader

import logging

logger = logging.getLogger(__name__)


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


def get_labels(mtz_file):
    """Function to get the column labels for input mtz file

    Parameters
    ----------
    mtz_file : str
       The path to the mtz file

    Returns
    -------
    f : str
        f column label
    fp : str
        fp column label
    dano : str
        dano column label
    sigdano : str
        sigdano column label
    free : str
        free column label
    """

    reflection_file = reflection_file_reader.any_reflection_file(file_name=mtz_file)
    if not reflection_file.file_type() == "ccp4_mtz":
        msg="File is not of type ccp4_mtz: {0}".format(mtz_file)
        logging.critical(msg)
        raise RuntimeError(msg)

    content = reflection_file.file_content()
    ctypes = content.column_types()
    clabels = content.column_labels()
    ftype = 'F'
    dtype = 'D'

    if ftype not in ctypes:
        msg = "Cannot find any structure amplitudes in: {0}".format(mtz_file)
        raise RuntimeError(msg)
    f = clabels[ctypes.index(ftype)]

    # FP derived from F
    fp = 'SIG' + f
    if fp not in clabels:
        msg = "Cannot find label {0} in file: {1}".format(fp, mtz_file)
        raise RuntimeError(msg)

    try:
        if dtype not in ctypes:
            msg = "Cannot find any structure amplitudes in: {0}".format(mtz_file)
            raise RuntimeError(msg)
        dano = clabels[ctypes.index(dtype)]

        # SIGDANO derived from DANO
        sigdano = 'SIG' + dano
        if sigdano not in clabels:
            msg = "Cannot find label {0} in file: {1}".format(sigdano, mtz_file)
            raise RuntimeError(msg)
    except RuntimeError:
        dano, sigdano = None, None

    free = None
    for label in content.column_labels():
        if 'free' in label.lower():
            column = content.get_column(label=label)
            selection_valid = column.selection_valid()
            flags = column.extract_values()
            sel_0 = (flags == 0)
            # extract number of work/test reflections
            n0 = (sel_0 & selection_valid).count(True)
            n1 = (~sel_0 & selection_valid).count(True)
            if n0 > 0 and n1 > 0:
                if free:
                    logger.warning("FOUND >1 R FREE label in file!")
                free = label

    return f, fp, dano, sigdano, free

