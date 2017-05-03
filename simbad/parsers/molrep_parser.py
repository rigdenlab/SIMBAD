"""Returns scores from a molrep log file"""

__author__ = "Adam Simpkin"
__date__ = "03 May 2017"
__version__ = "0.1"


class MolrepParser(object):
    """Class to mine information from a molrep log file"""

    def __init__(self, logfile):
        self._logfile = None

        self.score = None
        self.tfscore = None
        self.time = None
        self.wrfac = None
        self.version = None

        self.logfile = logfile

        self.parse()

    @property
    def logfile(self):
        """Return logfile"""
        return self._logfile

    @logfile.setter
    def logfile(self, logfile):
        """Define logfile"""
        self._logfile = logfile

    def parse(self):
        """Parse information from the logfile"""
        with open(self.logfile) as f:
            line = f.readline()
            while line:
                if line.startswith(" ### CCP4") and "version" in line:
                    self.version = line.strip().split()[5]
                if "Nmon RF  TF   theta    phi     chi   tx     ty     tz     TF/sg  wRfac  Score" in line:
                    line = f.readline()
                    fields = line.strip().split()
                    self.tfscore = float(fields[9])
                    self.wrfac = float(fields[10])
                    self.score = float(fields[11])

                if line.startswith("Times: User:"):
                    fields = line.strip().split()
                    time = fields[6]
                    m, s = time.split(":")
                    self.time = int(m) * 60 + int(s)
                line = f.readline()
        return
