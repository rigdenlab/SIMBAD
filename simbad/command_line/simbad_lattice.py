#!/usr/bin/env python

__author__ = "Felix Simkovic & Adam Simpkin"
__date__ = "06 Mar 2017"
__version__ = "0.1"

import argparse
import os
import platform
import sys
import time

import simbad.command_line
import simbad.util.exit_util
import simbad.util.simbad_util
import simbad.version

__version__ = simbad.version.__version__


def lattice_argparse():
    """Create the argparse options"""
    p = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    simbad.command_line._argparse_core_options(p)
    simbad.command_line._argparse_lattice_options(p)
    simbad.command_line._argparse_mtz_options(p)
    simbad.command_line._argparse_mr_options(p)
    p.add_argument('mtz', help="The path to the input mtz file")
    return p.parse_args()
    

def main():
    """Main function to run SIMBAD's lattice search"""
    args = lattice_argparse()

    if args.work_dir and os.path.isdir(args.work_dir):
        raise ValueError("Named working directory exists, please rename or remove")
    elif args.work_dir:
        os.mkdir(args.work_dir)
        args.work_dir = args.work_dir
    elif args.run_dir and os.path.isdir(args.run_dir):
        args.work_dir = simbad.command_line.make_workdir(args.run_dir, ccp4_jobid=args.ccp4_jobid)
    elif args.run_dir:
        os.mkdir(args.run_dir)
        args.work_dir = simbad.command_line.make_workdir(args.run_dir, ccp4_jobid=args.ccp4_jobid)
    else:
        raise RuntimeError("Not entirely sure what has happened here but I should never get to here")

    # Logger setup
    debug_log = os.path.join(args.work_dir, 'debug.log')
    logger = simbad.command_line.setup_logging(logfile=debug_log)

    # Check the CCP4 installation
    ccp4_root = simbad.command_line.setup_ccp4()
    ccp4_version = simbad.util.simbad_util.ccp4_version()
                                
    # Print some fancy info
    logger.info(simbad.command_line.header)
    logger.info("SIMBAD version: %s", __version__)
    logger.info("Running with CCP4 version: %s from directory: %s", ccp4_version, ccp4_root)
    logger.info("Running on host: %s", platform.node())
    logger.info("Running on platform: %s", platform.platform())
    logger.info("Job started at: %s", time.strftime("%a, %d %b %Y %H:%M:%S", time.gmtime()))
    logger.info("Invoked with command-line:\n%s\n", " ".join(map(str, sys.argv)))
    logger.info("Running in directory: %s\n", args.work_dir)

    # Take the start time
    time_start = time.time()

    # Perform the contaminante search
    solution_found = simbad.command_line._simbad_lattice_search(args)
    if solution_found:
        logger.info("Lucky you! SIMBAD worked its charm and found a lattice match for you.")
    else:
        logger.info("No results found - lattice search was unsuccessful")

    # Calculate and display the runtime in hours
    days, hours, mins, secs = simbad.command_line.calculate_runtime(time_start, time.time())
    logger.info("All processing completed in %d days, %d hours, %d minutes, %d and seconds", days, hours, mins, secs)


if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        msg = "Error running main SIMBAD program: {0}".format(e.message)
        simbad.util.exit_util.exit_error(msg, sys.exc_info()[2])

