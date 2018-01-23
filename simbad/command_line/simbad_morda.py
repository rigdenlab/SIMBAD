#!/usr/bin/env python

__author__ = "Adam Simpkin & Felix Simkovic"
__date__ = "16 May 2017"
__version__ = "0.1"

import argparse
import os
import sys

from pyjob.misc import StopWatch

import simbad.command_line
import simbad.exit
import simbad.util
import simbad.util.pyrvapi_results

logger = None


def morda_argparse():
    """Create the argparse options"""
    p = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
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

    args.work_dir = simbad.command_line.get_work_dir(
        args.run_dir, work_dir=args.work_dir
    )

    log_file = os.path.join(args.work_dir, 'simbad.log')
    debug_log_file = os.path.join(args.work_dir, 'debug.log')
    log_class = simbad.command_line.LogController()
    log_class.add_console(level=args.debug_lvl)
    log_class.add_logfile(log_file, level="info", format="%(message)s")
    log_class.add_logfile(debug_log_file, level="notset",
                          format="%(asctime)s\t%(name)s [%(lineno)d]\t%(levelname)s\t%(message)s")
    global logger
    logger = log_class.get_logger()

    if not os.path.isfile(args.amore_exe):
        raise OSError("amore executable not found")

    gui = simbad.util.pyrvapi_results.SimbadOutput(
        args.rvapi_document, args.webserver_uri, args.display_gui, log_file, args.work_dir, ccp4i2_xml=args.ccp4i2_xml
    )

    simbad.command_line.print_header()
    logger.info("Running in directory: %s\n", args.work_dir)

    stopwatch = StopWatch()
    stopwatch.start()

    solution_found = simbad.command_line._simbad_morda_search(args)
    if solution_found:
        logger.info("Solution found in SIMBAD's MoRDa search")
    else:
        logger.info("No results found - full search was unsuccessful")

    if args.output_pdb and args.output_mtz:
        run_dir = os.path.join(args.work_dir, 'morda')
        csv = os.path.join(run_dir, 'morda_mr.csv')
        result = simbad.util.result_by_score_from_csv(csv, 'final_r_free', ascending=True)
        simbad.util.output_files(run_dir, result, args.output_pdb, args.output_mtz)

    stopwatch.stop()
    logger.info("All processing completed in %d days, %d hours, %d minutes, and %d seconds",
                *stopwatch.time_pretty)

    gui.display_results(True, args.results_to_display)
    if args.rvapi_document:
        gui.save_document()
    log_class.close()


if __name__ == "__main__":
    try:
        main()
    except Exception:
        simbad.exit.exit_error(*sys.exc_info())
