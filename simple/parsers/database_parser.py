'''
Functions for retrieving information from the PDB database
@author hlasimpk
'''

import os
import sys
import time

def printProgress (iteration, total, prefix = '', suffix = '', decimals = 2, barLength = 100):
    """
    Call in a loop to create terminal progress bar
    @params:
        iteration   - Required  : current iteration (Int)
        total       - Required  : total iterations (Int)
        prefix      - Optional  : prefix string (Str)
        suffix      - Optional  : suffix string (Str)
    """
    filledLength    = int(round(barLength * iteration / float(total)))
    percents        = round(100.00 * (iteration / float(total)), decimals)
    bar             = '#' * filledLength + '-' * (barLength - filledLength)
    sys.stdout.write('%s [%s] %s%s %s\r' % (prefix, bar, percents, '%', suffix)),
    sys.stdout.flush()
    if iteration == total:
        print("\n")

def get_pdb_dirs(directory):
    first = True
    pdb_dirs = []
    for e in os.walk(directory):
        if first:
            first = False
            continue
        else:
            pdb_dirs.append(e[0][-2:])
    return pdb_dirs

def get_pdb_names(directory):
    pdb_name = []
    for e in os.walk(directory):
        for f in e[2]:
            if ".pdb" in f:
                basename = os.path.splitext(f)
                pdb_name.append(basename[0])
    return pdb_name

def pdb_information(directory):
    pdb_name = get_pdb_names(directory)
    pdb_info = {}
    for pdb in pdb_name:
        pdb_files = []
        for e in os.walk(directory):
            for f in e[2]:
                if pdb in f:
                    pdb_files.append(f)
        pdb_info[pdb] = pdb_files
    return pdb_info

def directory_information(PDB_directory):
    pdb_dirs = get_pdb_dirs(PDB_directory)
    pdb_info = {}
    item = 0
    length = len(pdb_dirs)
    print "Creating dictionary of files for {0}".format(PDB_directory) + os.linesep
    time_start = time.time()
    printProgress(item, length, prefix = 'Progress:', suffix = 'Complete', barLength = 50)
    for i in pdb_dirs:
            directory = os.path.join(PDB_directory, i)
            dictionary = pdb_information(directory)
            pdb_info.update(dictionary)
            item += 1
            printProgress(item, length, prefix = 'Progress:', suffix = 'Complete', barLength = 50)
    
    time_stop = time.time()
    elapsed_time = time_stop - time_start
    run_in_min = elapsed_time / 60
    run_in_hours = run_in_min / 60
    msg = os.linesep + 'Dictionary created (in {0:6.2F} hours)'.format(run_in_hours) + os.linesep
    msg += '----------------------------------------' + os.linesep
    print msg
