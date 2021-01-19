"""Returns scores from a molrep log file"""

__author__ = "Adam Simpkin"
__date__ = "03 May 2017"
__version__ = "0.1"

import simbad.parsers


class MolrepParser(simbad.parsers._Parser):
    """Class to mine information from a molrep log file"""

    def __init__(self, logfile):
        super(MolrepParser, self).__init__(logfile)
        self.contrast = None
        self.score = None
        self.tfscore = None
        self.time = None
        self.wrfac = None
        self.version = None
        self.parse()

    def parse(self):
        with open(self.fname) as f:
            line = f.readline()
            while line:
                if line.startswith(" ### CCP4") and "version" in line:
                    self.version = line.strip().split()[5]
                if "TF/sig" in line and "=" in line:
                    try:
                        self.tfscore = float(line.strip().split()[-1])
                    except ValueError:
                        continue
                if "Final CC" in line and "=" in line:
                    try:
                        self.score = float(line.strip().split()[-1])
                    except ValueError:
                        continue
                if "Contrast" in line and "=" in line:
                    try:
                        self.contrast = float(line.strip().split()[-1])
                    except ValueError:
                        continue

                if "Nmon RF  TF   theta    phi     chi   tx     ty     tz     TF/sg  wRfac  Score" in line:
                    line = f.readline()
                    try:
                        self.wrfac = float(line[65:72])
                    except ValueError:
                        continue

                if line.startswith("Times: User:"):
                    fields = line.strip().split()
                    time = fields[6]
                    m, s = time.split(":")
                    self.time = int(m) * 60 + int(s)
                line = f.readline()

    def summary(self):
        pass

    def check_input(self):
        pass

