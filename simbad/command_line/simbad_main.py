#!/usr/bin/env python

from __future__ import print_function

__author__ = "Adam Simpkin, and Felix Simkovic"
__contributing_authors__ = "Jens Thomas, and Ronan Keegan"
__credits__ = "Daniel Rigden, William Shepard, Charles Ballard, Villi Uski, and Andrey Lebedev"
__date__ = "17 May 2017"
__email__ = "hlasimpk@liv.ac.uk"
__version__ = "0.1"

import argparse
import os
import sys
import time

import simbad.command_line
import simbad.exit
import simbad.util.pyrvapi_results

logger = None


def simbad_argparse():
    """Create the argparse options"""
    p = argparse.ArgumentParser(
        description="SIMBAD: Sequence Independent Molecular replacement Based on Available Database",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    simbad.command_line._argparse_core_options(p)
    simbad.command_line._argparse_job_submission_options(p)
    simbad.command_line._argparse_contaminant_options(p)
    simbad.command_line._argparse_lattice_options(p)
    simbad.command_line._argparse_rot_options(p)
    simbad.command_line._argparse_mr_options(p)
    simbad.command_line._argparse_mtz_options(p)
    p.add_argument('mtz', help="The path to the input mtz file")
    return p


def main():
    """Main SIMBAD routine"""
    args = simbad_argparse().parse_args()

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
    
    # Account for the fact that argparse can't take bool
    if args.early_term.lower() == 'false':
        args.early_term = False
    
    # Logger setup
    global logger
    debug_log = os.path.join(args.work_dir, 'debug.log')
    logger = simbad.command_line.setup_logging(level=args.debug_lvl, logfile=debug_log)
    
    #GUI setup
    SR = simbad.util.pyrvapi_results.SimbadOutput()
    SR.display_results(args.webserver_uri, args.no_gui, debug_log, args.work_dir, summary=False)

    # Print some nice information
    simbad.command_line.print_header()
    logger.info("Running in directory: %s\n", args.work_dir)
    
    # Take the start time
    time_start = time.time()

    # Let's start searching old boyo
    end_of_cycle, solution_found = False, False
    while not (solution_found or end_of_cycle):
        # =====================================================================================
        # Perform the lattice search
        solution_found = simbad.command_line._simbad_lattice_search(args)
        if solution_found and args.early_term:
            logger.info("Lucky you! SIMBAD worked its charm and found a lattice match for you.")
            continue
        else:
            logger.info("No results found - lattice search was unsuccessful")
        
        SR.display_results(args.webserver_uri, args.no_gui, debug_log, args.work_dir, summary=False)
        
        # =====================================================================================
        # Perform the contaminant search
        solution_found = simbad.command_line._simbad_contaminant_search(args)   
        if solution_found:
            logger.info("Check you out, crystallizing contaminants! But don't worry, SIMBAD figured it out and found a solution.")
            continue
        else:
            logger.info("No results found - contaminant search was unsuccessful")
        
        SR.display_results(args.webserver_uri, args.no_gui, debug_log, args.work_dir, summary=False)
        # =====================================================================================
        # Make sure we only run the loop once for now
        end_of_cycle = True

    # Calculate and display the runtime in hours
    days, hours, mins, secs = simbad.command_line.calculate_runtime(time_start, time.time())
    logger.info("All processing completed in %d days, %d hours, %d minutes, and %d seconds", days, hours, mins, secs)
    
    # Output summary in GUI
    SR.display_results(args.webserver_uri, args.no_gui, debug_log, args.work_dir, summary=True)

if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        msg = "Error running main SIMBAD program: {0}".format(e.message)
        simbad.exit.exit_error(*sys.exc_info())
