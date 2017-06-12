#!/bin/bash

# Written by: Adam Simpkin
# Date:       18 May 2017
#
# This script is an example of how to run the a MoRDa search in 
# SIMBAD. Various command line flags could be set but below's 
# example is the most basic usage and illustrates its purpose
#.
# Note: the MoRDa database needs to be installed separately.

morda_db="<PATH TO MORDA DB>"
nproc=4

simbad-morda \
    -morda_db ${morda_db} \
    -nproc ${nproc} \
    input/1dtx.mtz
