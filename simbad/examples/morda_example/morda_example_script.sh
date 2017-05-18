#!/usr/bin/env ccp4-python

"""Example script to run SIMBAD MoRDa search on 5UOU test case, note: need to provide location of local morda db install"""

morda_db="<input_path_to_morda_database>"
nproc=8

simbad-morda -morda_db ${morda_db} -nproc ${nproc} 5uou.mtz
