import base64
import os
import zlib


def _from_dat(fhandle):
    """Little decompression/decoding function for SIMBAD .dat files"""
    return zlib.decompress(
        base64.b64decode(
            fhandle.read()
        )
    )


def _to_dat(fhandle):
    """Little compression/encoding function for SIMBAD .dat files"""
    return base64.b64encode(
        zlib.compress(
            fhandle.read()
        )
    )


def find_simbad_dat_files(directory):
    """Find all SIMBAD morda files

    Parameters
    ----------
    directory : str
       Path to the SIMBAD database

    Returns
    -------
    list
       A list of paths to the files
    """
    return [
        os.path.join(root, filename) for root, _, files in os.walk(directory)
        for filename in files if filename.endswith('.dat')
    ]


def convert_dat_to_pdb(infile, outfile):
    """Function to move dat file to pdb file

    Parameters
    ----------
    infile : str
        Path to the input dat file
    outfile : str
        Path to the output pdb file

    """
    with open(infile, 'rb') as f_in, open(outfile, 'w') as f_out:
        f_out.write(_from_dat(f_in))


def convert_pdb_to_dat(infile, outfile):
    """Function to move pdb file to dat file

    Parameters
    ----------
    infile : str
        Path to the input pdb file
    outfile : str
        Path to the output dat file

    """
    with open(output, "r") as f_in, open(final, "wb") as f_out:
        f_out.write(_to_dat(f_in))


def is_valid_dat(infile):
    """Validate a dat file for decompression/decoding

    Parameters
    ----------
    infile : str
       Path to the input pdb file

    Returns
    -------
    bool

    """
    with open(infile, "rb") as f_in:
        try:
            _from_dat(f_in)
            is_valid = True
        except:
            is_valid = False
    return is_valid


def read_dat(infile):
    """Read a SIMBAD .dat file

    Parameters
    ----------
    infile : str
       Path to the input pdb file

    Returns
    -------
    bool

    """
    with open(infile, "rb") as f_in:
        s = _from_dat(f_in)
    return s
