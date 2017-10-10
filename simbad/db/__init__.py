import base64
import zlib


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
        f_out.write(
            zlib.decompress(
                base64.b64decode(
                    f_in.read()
                )
            )
        )


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
        f_out.write(
            base64.b64encode(
                zlib.compress(
                    f_in.read()
                )
            )
        )


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
    with open(infile, "r") as f_in:
        try:
            zlib.decompress(base64.b64decode(f_in.read()))
            is_valid = True
        except:
            is_valid = False
    return is_valid
