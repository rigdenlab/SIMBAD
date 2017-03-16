#!/usr/bin/env python

__author__ = "Adam Simpkin"
__date__ = "08 Mar 2017"
__version__ = "0.1"

import argparse
import logging
import os

import simbad.constants
import simbad.rotsearch.amore_search
import simbad.util.mr_util

logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.DEBUG)

def main():
    p = argparse.ArgumentParser()
    p.add_argument('-amore_exe', default=os.path.join(os.environ["CCP4"], 'bin', 'amoreCCB2.exe'),
                   help="Path to amore exectutable")
    p.add_argument("-npic", default=50,
                   help="Number of peaks to output from the translation function map for each orientation")
    p.add_argument('-nproc', default=2, type=int,
                   help="Number of processors")
    p.add_argument('-max_to_keep', default=20,
                   help="The maximum number of results to return")
    p.add_argument('-min_solvent_content', default=30,
                   help="The minimum solvent content present in the unit cell with the input model")
    p.add_argument('-models_dir', default=simbad.constants.CONTAMINANT_MODELS,
                   help="The path to the directory containing the search models")
    p.add_argument('-mr_program', default='molrep',
                   help="The name of the molecular replacement program to use <molrep|phaser>")
    p.add_argument('mtz', help="The path to the input mtz file")
    p.add_argument('-output_dir', default=os.path.join(os.getcwd(), 'mr_output'),
                   help="The output directory")
    p.add_argument('-pklim', default=0.5,
                   help="Peak limit, output all peaks above <pklim>")
    p.add_argument('-refine_program', default='refmac5',
                   help="The name of the refinement program to use <refmac5>")
    p.add_argument('-rotastep', default=1.0,
                   help="Size of rotation step")
    p.add_argument('-shres', default=3.0,
                   help="Spherical harmonic resolution")
    p.add_argument('-work_dir', default=os.getcwd(),
                   help="The path to the directory where you want SIMBAD to run")
    args = p.parse_args()

    amore_exe = args.amore_exe
    mtz = args.mtz
    work_dir = args.work_dir
    max_to_keep = args.max_to_keep
    models_dir = args.models_dir
    logs_dir = os.path.join(work_dir, 'clogs')
    nproc = args.nproc
    shres = args.shres
    pklim = args.pklim
    npic = args.npic
    rotastep = args.rotastep
    min_solvent_content = args.min_solvent_content

    rotation_search = simbad.rotsearch.amore_search.AmoreRotationSearch(amore_exe, mtz, work_dir, max_to_keep)
    rotation_search.sortfun()
    rotation_search.amore_run(models_dir, logs_dir, nproc, shres, pklim, npic, rotastep, min_solvent_content)
    rotation_search.summarize()

    # MR with defaults for now
    model_dir = os.path.join(args.models_dir)
    molecular_replacement = simbad.util.mr_util.MrSubmit(args.mtz, args.mr_program, args.refine_program, model_dir,
                                                         args.output_dir)
    molecular_replacement.multiprocessing(rotation_search.search_results, nproc=args.nproc)
    molecular_replacement.summarize()

if __name__ == "__main__":
    main()