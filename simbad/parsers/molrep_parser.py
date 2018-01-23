"""Returns scores from a molrep log file"""

__author__ = "Adam Simpkin"
__date__ = "03 May 2017"
__version__ = "0.1"

import simbad.parsers


class MolrepParser(simbad.parsers._Parser):
    """Class to mine information from a molrep log file"""

    def __init__(self, logfile):
        super(MolrepParser, self).__init__(logfile)
        self.score = None
        self.tfscore = None
        self.time = None
        self.wrfac = None
        self.version = None
        self._parse()

    def _parse(self):
        with open(self.logfile) as f:
            line = f.readline()
            while line:
                if line.startswith(" ### CCP4") and "version" in line:
                    self.version = line.strip().split()[5]
                if "Nmon RF  TF   theta    phi     chi   tx     ty     tz     TF/sg  wRfac  Score" in line:
                    line = f.readline()
                    self.tfscore = float(line[58:65])
                    self.wrfac = float(line[65:72])
                    self.score = float(line[72:])

                if line.startswith("Times: User:"):
                    fields = line.strip().split()
                    time = fields[6]
                    m, s = time.split(":")
                    self.time = int(m) * 60 + int(s)
                line = f.readline()
