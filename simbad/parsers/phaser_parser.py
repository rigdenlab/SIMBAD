"""Returns scores from a phaser log file"""

__author__ = "Adam Simpkin"
__date__ = "03 May 2017"
__version__ = "0.1"

import simbad.parsers


class PhaserParser(simbad.parsers._Parser):
    def __init__(self, logfile):
        super(PhaserParser, self).__init__(logfile)
        self.llg = 0
        self.tfz = 0
        self.rfz = 0
        self.parse()

    def parse(self):
        with open(self.fname) as f:
            for line in f:
                if line.startswith("   SOLU SET") and "TFZ=" in line:
                    llist = line.split()
                    llist.reverse()
                    for i in llist:
                        if all(x not in i for x in ['*', '(', ')']) and "TFZ==" in i:
                            self.tfz = float(i.replace("TFZ==", ""))
                            break
                        if all(x not in i for x in ['*', '(', ')', 'TFZ==']) and "TFZ=" in i:
                            self.tfz = float(i.replace("TFZ=", ""))
                            break

                    for i in llist:
                        if "LLG==" in i:
                            self.llg = float(i.replace("LLG==", ""))
                            break
                        if "LLG=" in i and "LLG==" not in i:
                            self.llg = float(i.replace("LLG=", ""))
                            break

                    for i in llist:
                        if "RFZ==" in i:
                            self.rfz = float(i.replace("RFZ==", ""))
                            break
                        if "RFZ=" in i and "RFZ==" not in i:
                            self.rfz = float(i.replace("RFZ=", ""))
                            break
