"""Returns values for the top anomalous peak from a anode lsa file"""

__author__ = "Adam Simpkin"
__date__ = "13 April 2018"
__version__ = "0.1"

import simbad.parsers


class AnodeParser(simbad.parsers._Parser):
    """Class to mine information from a anode lsa"""

    def __init__(self, logfile):
        super(AnodeParser, self).__init__(logfile)
        self.x = None
        self.y = None
        self.z = None
        self.peak_height = None
        self.nearest_atom = None
        self._parse()

    def _parse(self):
        with open(self.logfile, "r") as f:
            line = f.readline()
            while line:
                if "          X        Y        Z   Height(sig)  SOF     Nearest atom" in line:
                    f.readline()
                    line = f.readline()
                    fields = line.split()
                    self.x = fields[1]
                    self.y = fields[2]
                    self.z = fields[3]
                    self.peak_height = fields[4]
                    self.nearest_atom = fields[-1]
                line = f.readline()
