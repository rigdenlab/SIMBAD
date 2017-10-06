#!/usr/bin/env python

__author__ = "Felix Simkovic & Adam Simpkin"
__date__ = "06 Mar 2017"
__version__ = "0.1"

import argparse
import os
import sys

from pyjob.misc import StopWatch

import simbad.command_line
import simbad.exit
import simbad.util.pyrvapi_results

logger = None


def lattice_argparse():
    """Create the argparse options"""
    prep = argparse.ArgumentParser(add_help=False)
    prep.add_argument('-sg', dest="space_group", type=str, default=None,
                      help='The space group to use')
    prep.add_argument('-uc', dest="unit_cell", type=str, default=None,
                      help="The unit cell, format 'a,b,c,alpha,beta,gamma'")
    args, _ = prep.parse_known_args()

    p = argparse.ArgumentParser(
        parents=[prep], formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    simbad.command_line._argparse_core_options(p)
    simbad.command_line._argparse_job_submission_options(p)
    simbad.command_line._argparse_lattice_options(p)
    simbad.command_line._argparse_mtz_options(p)
    simbad.command_line._argparse_mr_options(p)
    if args.space_group and args.unit_cell:
        # Add to the namespace as we're looking for it later
        p.add_argument('-mtz', dest="mtz", default=None,
                       help=argparse.SUPPRESS)
    else:
        p.add_argument('mtz', help="The path to the input mtz file")
    return p


def main():
    """Main function to run SIMBAD's lattice search"""
    args = lattice_argparse().parse_args()

    args.work_dir = simbad.command_line.get_work_dir(
        args.run_dir, work_dir=args.work_dir, ccp4_jobid=args.ccp4_jobid
    )

    global logger
    log = os.path.join(args.work_dir, 'simbad.log')
    debug_log = os.path.join(args.work_dir, 'debug.log')
    logger = simbad.command_line.setup_logging(level=args.debug_lvl, logfile=log,
                                               debug_logfile=debug_log)

    gui = simbad.util.pyrvapi_results.SimbadOutput(
        args.rvapi_document, args.webserver_uri, args.display_gui, log, args.work_dir
    )

    simbad.command_line.print_header()
    logger.info("Running in directory: %s\n", args.work_dir)

    stopwatch = StopWatch()
    stopwatch.start()

    solution_found = simbad.command_line._simbad_lattice_search(args)
    if args.space_group and args.unit_cell:
        display_summary = False
    elif solution_found:
        logger.info(
            "Lucky you! SIMBAD worked its charm and found a lattice match for you.")
        display_summary = True
    else:
        logger.info("No results found - lattice search was unsuccessful")
        display_summary = True

    stopwatch.stop()
    logger.info("All processing completed in %d days, %d hours, %d minutes, and %d seconds",
                *stopwatch.time_pretty)

    gui.display_results(display_summary)
    if args.rvapi_document:
        gui.save_document()


if __name__ == "__main__":
    try:
        main()
    except Exception:
        simbad.exit.exit_error(*sys.exc_info())
