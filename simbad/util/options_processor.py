"""
Code to check options added to SIMBAD

@author: hlasimpk
"""

import logging
import os

import exit_util
import mtz_util

logger = logging.getLogger(__name__)

def check_mandatory_options(optd):
    """
    Check that the mandatory options have been submitted and check for correctness
    :param optd:
    :return:
    """

    def _exit(msg, wdir):
        exit_util.exit_error(msg)

    if not 'mtz' in optd.keys() or not optd['mtz'] or not os.path.exists(optd['mtz']):
        msg = "MTZ not defined"
        _exit(msg)

    mtz_util.processReflectionFile(optd)
    if not optd['FREE']:
        msg = "FREE flag missing from MTZ"
        _exit(msg)

    if (optd['mtz'] and optd['sf_cif']):
        msg = "Please supply a single crystallographic data file."
        _exit(msg, optd['work_dir'])

    ####################################################################################################################
    # NEED TO ADD A SECTION HERE TO CHECK THAT DATABASES ARE INSTALLED
    ####################################################################################################################

    return optd

def process_options(optd):
    """
    Process input options

    :param optd:
    :return:
    """

    # Check if there is anomalous signal in the MTZ file and all options have been provided
    if optd['DANO'] != None and optd['SIGDANO'] != None and optd['hatom_num'] and optd['hatom_type']:
        optd['anomalous_signal'] = True

    # Get crystal data
    mtz_util.set_crystal_data(optd)

    return optd



