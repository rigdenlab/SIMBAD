"""Exit utility for catching errors and printing unified error messages"""

__author__ = "Jens Thomas & Felix Simkovic"
__date__ = "08 May 2017"
__version__ = "1.1"

import logging
import os
import sys
import traceback

try:
    import pyrvapi
except ImportError:
    pyrvapi = None


def _debug_logfile(logger):
    """Get the debug logfile"""
    if logger.handlers:
        for d in logger.handlers:
            if hasattr(d, "baseFilename") and d.level == logging.DEBUG:
                return getattr(d, "baseFilename")
    return None


def exit_error(exc_type, exc_value, exc_traceback):
    """Exit on error collecting as much information as we can.

    Parameters
    ----------
    exc_type : str
       The exception type
    exc_value : str
       The exception value
    exc_traceback
       The exception traceback
    
    Warnings
    --------
    This function terminates the program after printing appropriate
    error messages.
    
    """
    # Get the root logger
    logger = logging.getLogger(__name__)

    # Traceback info
    traceback_value_msg = exc_value
    traceback_full_msg = traceback.format_exception(exc_type, exc_value, exc_traceback)

    # Find debug log file
    debug_log = _debug_logfile(logger)

    # Construct the message
    main_msg = (
        "%(sep)s%(hashish)s%(sep)s"
        + "%(short_hash)s%(msg)s%(short_hash)s%(sep)s"
        + "%(hashish)s%(sep)s%(sep)s"
        + "SIMBAD exited with message: %(tb_value)s"
        + "%(sep)s%(sep)s%(hashish)s%(sep)s%(sep)s"
    )
    if debug_log:
        main_msg += "More information may be found in the debug log file: %(logfile)s%(sep)s"
    main_msg += "%(sep)sIf you believe that this is an error with SIMBAD, please email: %(email)s%(sep)s"
    main_msg += "providing as much information as you can about how you ran the program.%(sep)s"
    if debug_log:
        main_msg += "%(sep)sPlease static the debug logfile with your email: %(logfile)s%(sep)s"

    nhashes = 70
    main_msg_kwargs = {
        "sep": os.linesep,
        "hashish": "*" * nhashes,
        "short_hash": "*" * 19,
        "msg": "SIMBAD_ERROR".center(32, " "),
        "tb_value": traceback_value_msg,
        "logfile": debug_log,
        "email": "ccp4@stfc.ac.uk",
    }

    # String it all together
    logger.critical(main_msg, main_msg_kwargs)

    logger.critical("SIMBAD EXITING AT...")
    logger.critical("".join(traceback_full_msg))

    # Make sure the error widget is updated
    if pyrvapi:
        pyrvapi.rvapi_flush()

    sys.exit(1)
