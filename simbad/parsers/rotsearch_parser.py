"""Returns values for the top rotation peak from a rotsearch log file"""

__author__ = "Adam Simpkin"
__date__ = "08 May 2017"
__version__ = "0.1"

import simbad.parsers


class RotsearchParser(simbad.parsers._Parser):
    """Class to mine information from a rotsearch logfile"""
    
    def __init__(self, logfile):
        super(RotsearchParser, self).__init__(logfile)
        
        self.alpha = None
        self.beta = None
        self.gamma = None
        self.cc_f = None
        self.rf_f = None
        self.cc_i = None
        self.cc_p = None
        self.icp = None
        self.cc_f_z_score = None
        self.cc_p_z_score = None
        self.num_of_rot = None
        
        self.parse(logfile)
        
    def parse(self, logfile):
        """Parse information from the logfile"""
        for line in open(logfile):
            if line.startswith(" SOLUTIONRCD "):
                fields = line.split()
                if float(fields[-3]) > 0:
                    try:
                        self.alpha = float(fields[2])
                        self.beta = float(fields[3])
                        self.gamma = float(fields[4])
                        self.cc_f = float(fields[8])
                        self.rf_f = float(fields[9])
                        self.cc_i = float(fields[10])
                        self.cc_p = float(fields[11])
                        self.icp = float(fields[12])
                        self.cc_f_z_score = float(fields[-3])
                        self.cc_p_z_score = float(fields[-2])
                        self.num_of_rot = float(fields[-1])

                    except ValueError: 
                        self.alpha = float(fields[2])
                        self.beta = float(fields[3])
                        self.gamma = float(fields[4])
                        self.cc_f_z_score = float(fields[-3])
                        self.cc_p_z_score = float(fields[-2])
                        self.num_of_rot = float(fields[-1])
                    break