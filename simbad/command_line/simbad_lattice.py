#!/usr/bin/env python

__author__ = "Felix Simkovic"
__date__ = "06 Mar 2017"
__version__ = "0.1"

import argparse
import logging

import simbad.constants
import simbad.lattice.search
import simbad.util.mtz_util

logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.DEBUG)


def main():
    p = argparse.ArgumentParser()
    p.add_argument('-d', '--database', default=simbad.constants.SIMBAD_LATTICE_DB,
                   help="The path to the LATTICE database")
    p.add_argument('mtz', help="The path to the input mtz file")
    args = p.parse_args()

    space_group, resolution, cell_parameters = simbad.util.mtz_util.crystal_data(args.mtz)
    lattice_database = args.database    

    lattice_search = simbad.lattice.search.LatticeSearch(cell_parameters, space_group, lattice_database)
    lattice_search.search()
    lattice_search.summarize()


if __name__ == "__main__":
    main()

