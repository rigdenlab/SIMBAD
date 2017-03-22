#!/usr/bin/env python

__author__ = "Felix Simkovic & Adam Simpkin"
__date__ = "06 Mar 2017"
__version__ = "0.1"

import argparse
import logging
import os

import simbad.constants
import simbad.lattice.search
import simbad.util.mtz_util
import simbad.util.mr_util

logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.DEBUG)


def main():
    p = argparse.ArgumentParser()
    p.add_argument('-d', '--database', default=simbad.constants.SIMBAD_LATTICE_DB,
                   help="The path to the LATTICE database")
    p.add_argument('-early_term', default=True,
                   help="Terminate the program early if a solution is found <True | False>")
    p.add_argument('-enam', default=False,
                   help="Trial enantimorphic space groups <True | False>")
    p.add_argument('-mr_program', default='molrep',
                   help="The name of the molecular replacement program to use <molrep|phaser>")
    p.add_argument('-refine_program', default='refmac5',
                   help="The name of the refinement program to use <refmac5>")
    p.add_argument('-working_dir', default=os.getcwd(),
                   help="The working directory")
    p.add_argument('-output_dir', default=os.path.join(os.getcwd(), 'output'),
                   help="The output directory")
    p.add_argument('-nproc', default=2, type=int,
                   help="Number of processors to run on")
    p.add_argument('mtz', help="The path to the input mtz file")
    args = p.parse_args()

    space_group, resolution, cell_parameters = simbad.util.mtz_util.crystal_data(args.mtz)
    lattice_database = args.database    

    lattice_search = simbad.lattice.search.LatticeSearch(cell_parameters, space_group, lattice_database)
    lattice_search.search()
    lattice_search.summarize()

    # Defaults put in just for now
    lattice_search.download_results(args.working_dir)
    model_dir = os.path.join(args.working_dir, 'lattice_input_models')
    molecular_replacement = simbad.util.mr_util.MrSubmit(args.mtz, args.mr_program, args.refine_program, model_dir,
                                                         args.output_dir, args.early_term, args.enam)
    molecular_replacement.multiprocessing(lattice_search.search_results, nproc=args.nproc)
    molecular_replacement.summarize()

if __name__ == "__main__":
    main()

