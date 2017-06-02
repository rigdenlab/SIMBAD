#!/usr/bin/env python

__author__ = "Adam Simpkin & Felix Simkovic"
__date__ = "16 May 2017"
__version__ = "0.1"

import argparse
import os
import sys

from mbkit.misc.stopwatch import StopWatch

import simbad.command_line
import simbad.exit
import simbad.util.pyrvapi_results

logger = None


def morda_argparse():
    """Create the argparse options"""
    p = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    simbad.command_line._argparse_core_options(p)
    simbad.command_line._argparse_job_submission_options(p)
    simbad.command_line._argparse_morda_options(p)
    simbad.command_line._argparse_rot_options(p)
    simbad.command_line._argparse_mtz_options(p)
    simbad.command_line._argparse_mr_options(p)
    p.add_argument('mtz', help="The path to the input mtz file")
    return p


def main():
    """Main function to run SIMBAD's MoRDa search"""
    args = morda_argparse().parse_args()

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
    global logger
    debug_log = os.path.join(args.work_dir, 'debug.log')
    logger = simbad.command_line.setup_logging(level=args.debug_lvl, logfile=debug_log)
    
    # GUI setup
    gui = simbad.util.pyrvapi_results.SimbadOutput()
    gui.display_results(args.webserver_uri, args.no_gui, debug_log, args.work_dir, summary=False)

    # Print some fancy info
    simbad.command_line.print_header()
    logger.info("Running in directory: %s\n", args.work_dir)

    # Start taking time
    stopwatch = StopWatch()
    stopwatch.start()
    
    # Perform the morda search
    solution_found = simbad.command_line._simbad_morda_search(args)
    if solution_found:
        logger.info("Solution found in SIMBAD's MoRDa search")
    else:
        logger.info("No results found - full search was unsuccessful")

    # Calculate and display the runtime in hours
    stopwatch.stop()
    runtime = StopWatch.convert(stopwatch.runtime)
    logger.info("All processing completed in %d days, %d hours, %d minutes, and %d seconds", *runtime)
    
    # Output summary in GUI
    gui.display_results(args.webserver_uri, args.no_gui, debug_log, args.work_dir, summary=True)


if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        msg = "Error running main SIMBAD program: {0}".format(e.message)
        simbad.exit.exit_error(*sys.exc_info())
