'''
Created on 2 Dec 2014

@author: jmht
'''
import clipper
import logging
import os
import shutil
import sys

from iotbx import reflection_file_reader

import mbkit.dispatch.cexectools
import simbad.util.simbad_util

logger = logging.getLogger()
logger.setLevel(logging.DEBUG)


# TODO: Get rid of this function completely
def set_crystal_data(optd):
    """Set crystallographic parameters from mtz file"""

    space_group, resolution, cell_parameters = crystal_data(optd['mtz'])

    optd['space_group'] = space_group
    optd['resolution'] = resolution
    optd['cell_parameters'] = cell_parameters

    return optd


def crystal_data(mtz):
    """Set crystallographic parameters from mtz file

    Parameters
    ----------
    mtz : str
       The path to the mtz file

    Returns
    -------
    space_group : str
       The space group
    resolution : str
       The resolution
    cell_parameters : str
       The cell parameters

    """

    hkl_info = clipper.HKL_info()
    mtz_file = clipper.CCP4MTZfile()
    mtz_file.open_read(mtz)
    mtz_file.import_hkl_info(hkl_info)

    sg, cell = hkl_info.spacegroup(), hkl_info.cell()

    space_group = sg.symbol_hm()
    space_group = space_group.replace(" ", "")

    resolution = "%.2lf" % hkl_info.resolution().limit()
    cell_parameters = "%.2lf %.2lf %.2lf %.2lf %.2lf %.2lf" % (cell.a(),
                                                               cell.b(),
                                                               cell.c(),
                                                               cell.alpha_deg(),
                                                               cell.beta_deg(),
                                                               cell.gamma_deg())

    return space_group, resolution, cell_parameters


def del_column(file_name, column, overwrite=True):
    """Delete a column from an mtz file and return a path to the file"""
    mtzDel = simbad.util.simbad_util.filename_append(file_name, "d{0}".format(column))
    cmd = ["mtzutils", "hklin1", file_name, "hklout", mtzDel]
    stdin = "EXCLUDE 1 {0}".format(column)
    _ = mbkit.dispatch.cexectools.cexec(cmd, stdin=stdin)
    if overwrite:
        shutil.move(mtzDel,file_name)
        return file_name
    else:
        return mtzDel

def add_rfree(file_name,directory=None,overwrite=True):
    """Run uniqueify on mtz file to generate RFREE data column"""
    mtzUnique = simbad.util.simbad_util.filename_append(file_name, "uniqueify", directory=directory)
    cmd = ['uniqueify', file_name, mtzUnique]
    _ = mbkit.dispatch.cexectools.cexec(cmd)
    if overwrite:
        shutil.move(mtzUnique,file_name)
        return file_name
    else:
        return mtzUnique

def get_labels(file_name):
    """Return the F, FP and FREE column labels"""

    reflection_file = reflection_file_reader.any_reflection_file(file_name=file_name)
    if not reflection_file.file_type()=="ccp4_mtz":
        msg="File is not of type ccp4_mtz: {0}".format(file_name)
        logging.critical(msg)
        raise RuntimeError(msg)

    content=reflection_file.file_content()
    ctypes=content.column_types()
    clabels=content.column_labels()
    ftype='F'
    dtype='D'

    if ftype not in ctypes:
        msg = "Cannot find any structure amplitudes in: {0}".format(file_name)
        raise RuntimeError(msg)
    F=clabels[ctypes.index(ftype)]

    # FP derived from F
    FP='SIG'+F
    if FP not in clabels:
        msg = "Cannot find label {0} in file: {1}".format(FP,file_name)
        raise RuntimeError(msg)

    try:
        if dtype not in ctypes:
            msg = "Cannot find any structure amplitudes in: {0}".format(file_name)
            raise RuntimeError(msg)
            
        DANO = clabels[ctypes.index(dtype)]

        # SIGDANO derived from DANO
        SIGDANO = 'SIG' + DANO
        if SIGDANO not in clabels:
            msg = "Cannot find label {0} in file: {1}".format(SIGDANO, file_name)
            raise RuntimeError(msg)
    except RuntimeError:
        DANO, SIGDANO = None, None


    FREE=_get_rfree(content)

    return F,FP,DANO,SIGDANO,FREE

def get_rfree(file_name):
    """Return the Rfree label"""

    reflection_file = reflection_file_reader.any_reflection_file(file_name=file_name)
    if not reflection_file.file_type()=="ccp4_mtz":
        msg = "File is not of type ccp4_mtz: {0}".format(file_name)
        logging.critical(msg)
        raise RuntimeError(msg)
    
    # Read the file
    content=reflection_file.file_content()
    return _get_rfree(content)
    
def _get_rfree(content):
    rfree_label=None
    #print "GOT ",content.column_labels()
    for label in content.column_labels():
        if 'free' in label.lower():
            column = content.get_column(label=label)
            selection_valid = column.selection_valid()
            flags = column.extract_values()
            sel_0 = (flags == 0)
            # extract number of work/test reflections
            n0=( sel_0 & selection_valid).count(True)
            n1=(~sel_0 & selection_valid).count(True)
            #print "Number of 0 (work):",n0
            #print "Number of 1 (test):",n1
            #print float(n0)/float(n1)*100
            if n0>0 and n1>0:
                if rfree_label:
                    logger.warning("FOUND >1 RFREE label in file!")
                rfree_label=label
    return rfree_label

def to_hkl(mtz_file,hkl_file=None,directory=None,F=None,SIGF=None,FREE=None):
    
    if directory is None:
        directory=os.getcwd()
    
    if hkl_file is None:
        name=os.path.splitext(os.path.basename(mtz_file))[0]
        hkl_file=os.path.join(directory,name+".hkl")
        
    if F is None or SIGF is None or FREE is None:
        F,SIGF,DANO,SIGDANO,FREE=get_labels(mtz_file)
        
    cmd=['mtz2various','HKLIN',mtz_file,'HKLOUT', hkl_file]
    stdin  = """LABIN FP={0} SIGFP={1} FREE={2}
        OUTPUT SHELX
        FSQUARED
        END
        """
    stdin = stdin.format(F,SIGF,FREE)
    _ = mbkit.dispatch.cexectools.cexec(cmd, stdin=stdin)
    return hkl_file

def processReflectionFile(amoptd):
    """Make sure we have a valid mtz file. If necessary convert a given cif file.
       Set the mtz variable in the given amoptd to the reflection file to use
       Return True if it all worked or raise an exception if it failed
    """

    # Now have an mtz so check it's valid
    if not amoptd['mtz'] or not os.path.isfile( amoptd['mtz'] ):
        logger.critical("Cannot find MTZ file: %s", amoptd['mtz'])
        sys.exit(1)


    # Get column label info
    reflection_file = reflection_file_reader.any_reflection_file(file_name=amoptd['mtz'])
    if not reflection_file.file_type()=="ccp4_mtz":
        logger.critical("File is not of type ccp4_mtz: %s", amoptd['mtz'])
        sys.exit(1)

    # Read the file
    content=reflection_file.file_content()

    # Check any user-given flags

    for flag in ['F','SIGF','FREE']:
        if amoptd[flag] and amoptd[flag] not in content.column_labels():
            logger.critical("Cannot find flag %s label %s in mtz file %s", flag, amoptd[flag], amoptd['mtz'])
            sys.exit(1)

    # If any of the flags aren't given we set defaults based on what's in the file
    if not amoptd['F']:
        if 'F' not in content.column_types():
            logger.critical("Cannot find column type F for flag F in mtz file: %s", amoptd['mtz'])
            sys.exit(1)
        amoptd['F']  = content.column_labels()[content.column_types().index('F')]
    if not amoptd['SIGF']:
        l='SIG'+amoptd['F']
        if l not in content.column_labels():
            logger.critical("Cannot find column type %s for flag SIGF in mtz file: %s", l, amoptd['mtz'])
            sys.exit(1)
        amoptd['SIGF']  = l

    if amoptd['FREE']:
        # Check is valid
        rfree=_get_rfree(content)
        if not rfree or not rfree==amoptd['FREE']:
            logger.critical("Given RFREE label %s is not valid for mtz file: %s", amoptd['FREE'], amoptd['mtz'])
            sys.exit(1)
    else:
        # See if we can find a valid label in the file
        rfree=_get_rfree(content)
        if not rfree:
            # Need to generate RFREE
            logger.warning("Cannot find a valid FREE flag - running uniquefy to generate column with RFREE data.")
            amoptd['mtz'] = add_rfree( amoptd['mtz'], directory=amoptd['work_dir'],overwrite=False)

            # Check file and get new FREE flag
            rfree=get_rfree(amoptd['mtz'])
            if not rfree:
                logger.critical("Cannot find valid rfree flag in mtz file %s after running uniquiefy", amoptd['mtz'])
                sys.exit(1)
        amoptd['FREE']  = rfree

    # Return anomalous data labels if the columns exist in the file
    if not amoptd['DANO']:
        if 'D' not in content.column_types():
            logger.critical("Cannot find column type D for flag DANO in mtz file: %s", amoptd['mtz'])
        else:
            amoptd['DANO'] = content.column_labels()[content.column_types().index('D')]

    if not amoptd['SIGDANO'] and amoptd['DANO']:
        l = 'SIG' + amoptd['DANO']
        if not l in content.column_labels():
            logger.critical("Cannot find column type %s for flag SIGDANO in mtz file: %s", l, amoptd['mtz'])
        else:
            amoptd['SIGDANO'] = l

    return True

