# Written by: Adam Simpkin
# Date:       21 July 2020
#
# This script is an example of how to run the contaminant
# search in SIMBAD. Various command line flags could be set
# but below's example is the most basic usage and illustrates
# its purporse.

%CCP4%\bin\simbad-contaminant.bat ^
-nproc 4 ^
-organism CHICK ^
input/2fbb.mtz