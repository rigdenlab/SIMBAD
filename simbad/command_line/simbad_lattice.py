#!/usr/bin/env python

__author__ = "Felix Simkovic & Adam Simpkin"
__date__ = "06 Mar 2017"
__version__ = "0.1"

import argparse
import os
import sys
import time

import simbad.command_line
import simbad.exit
import simbad.util.pyrvapi_results

logger = None


def lattice_argparse():
    """Create the argparse options"""
    p = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    simbad.command_line._argparse_core_options(p)
    simbad.command_line._argparse_job_submission_options(p)
    simbad.command_line._argparse_lattice_options(p)
    simbad.command_line._argparse_mtz_options(p)
    simbad.command_line._argparse_mr_options(p)
    p.add_argument('mtz', help="The path to the input mtz file")
    return p
    

def main():
    """Main function to run SIMBAD's lattice search"""
    args = lattice_argparse().parse_args()

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
    
    #GUI setup
    SR = simbad.util.pyrvapi_results.SimbadOutput()
    SR.display_results(args.webserver_uri, args.no_gui, debug_log, args.work_dir, summary=False)
    
    # Print some fancy info
    simbad.command_line.print_header()
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
    logger.info("All processing completed in %d days, %d hours, %d minutes, and %d seconds", days, hours, mins, secs)
    
    # Output summary in gui
    SR.display_results(args.webserver_uri, args.no_gui, debug_log, args.work_dir, summary=True)


if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        msg = "Error running main SIMBAD program: {0}".format(e.message)
        simbad.exit.exit_error(*sys.exc_info())

