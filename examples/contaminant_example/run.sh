#!/bin/bash

# Written by: Adam Simpkin
# Date:       18 May 2017
#
# This script is an example of how to run the contaminant
# search in SIMBAD. Various command line flags could be set
# but below's example is the most basic usage and illustrates
# its purporse.

simbad-contaminant \
    -nproc 10 \
    input/2fbb.mtz
