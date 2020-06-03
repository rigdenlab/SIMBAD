"""Returns MTZ info"""

__author__ = "Adam Simpkin & Filomeno Sanchez Rodriguez"
__date__ = "27 May 2020"
__version__ = "0.1"


import re
import gemmi
from enum import Enum

import simbad.parsers


class MtzColumnLabels(Enum):
    """An enumerator that contains the regular expression used to detect the column labels of a given MTZ file"""

    free = re.compile(r"^.*?[Ff]?[Rr]?[Ee]?[Ee]?.*?")
    i = re.compile(r"^[Ii]")
    sigi = re.compile(r"^[Ss][Ii][Gg][Ii]")
    f = re.compile(r"^[Ff][Pp]?(?![Cc])(?![Ww][Tw])")
    sigf = re.compile(r"^[Ss][Ii][Gg][Ff][Pp]?")
    i_plus = re.compile(r"^[Ii].*(\(\+\)|[Pp][Ll][Uu][Ss])")
    sigi_plus = re.compile(r"^[Ss][Ii][Gg][Ii].*(\(\+\)|[Pp][Ll][Uu][Ss])")
    f_plus = re.compile(r"^[Ff][Pp]?.*(\(\+\)|[Pp][Ll][Uu][Ss])")
    sigf_plus = re.compile(r"^[Ss][Ii][Gg][Ff][Pp]?.*(\(\+\)|[Pp][Ll][Uu][Ss])?")
    i_minus = re.compile(r"^[Ii].*(\(-\)|[Mm][Ii][Nn][Uu][Ss])")
    sigi_minus = re.compile(r"^[Ss][Ii][Gg][Ii].*(\(-\)|[Mm][Ii][Nn][Uu][Ss])")
    f_minus = re.compile(r"^[Ff][Pp]?.*(\(-\)|[Mm][Ii][Nn][Uu][Ss])")
    sigf_minus = re.compile(r"^[Ss][Ii][Gg][Ff][Pp]?.*(\(-\)|[Mm][Ii][Nn][Uu][Ss])")
    dp = re.compile(r"^[Dd][An]?[Nn]?[Oo]?[Pp]?")
    sigdp = re.compile(r"^[Ss][Ii][Gg][Dd][Aa]?[Nn]?[Oo]?[Pp]?")


class MTZColumnTypes(Enum):
    """An enumerator that contains the different types expected for each column of a given MTZ file"""

    free = 'I'
    i = 'J'
    sigi = 'Q'
    f = 'F'
    sigf = 'Q'
    i_plus = 'K'
    sigi_plus = 'M'
    f_plus = 'G'
    sigf_plus = 'L'
    i_minus = 'K'
    sigi_minus = 'M'
    f_minus = 'G'
    sigf_minus = 'L'
    dp = 'D'
    sigdp = 'Q'


class MtzParser(simbad.parsers._Parser):
    def __init__(self, fname, logger=None):

        super(MtzParser, self).__init__(fname, logger=logger)

        self.reflection_file = None
        self.cell = None
        self.f = None
        self.sigf = None
        self.i = None
        self.sigi = None
        self.free = None
        self.dp = None
        self.sigdp = None
        self.f_plus = None
        self.sigf_plus = None
        self.i_plus = None
        self.sigi_plus = None
        self.f_minus = None
        self.sigf_minus = None
        self.i_minus = None
        self.sigi_minus = None
        self.resolution = None
        self.nreflections = None
        self.spacegroup_symbol = None
        self.read_reflections()

    @property
    def summary(self):
        """Tuple with all the parsed label names"""
        return (self.f, self.sigf, self.i, self.sigi, self.free, self.dp, self.sigdp, self.f_plus, self.sigf_plus,
                self.i_plus, self.sigi_plus, self.f_minus, self.sigf_minus, self.i_minus, self.sigi_minus)

    def read_reflections(self):
        """Read the reflection data data"""

        self.reflection_file = gemmi.read_mtz_file(self.fname)
        self.nreflections = self.reflection_file.nreflections
        self.spacegroup_symbol = self.reflection_file.spacegroup.hm
        self.resolution = self.reflection_file.resolution_high()
        self.cell = self.reflection_file.cell

    def parse(self):
        """Parse the input mtz file and retrieve the column labels"""

        if self.error:
            self.logger.warning("Previous errors prevent parsing mtz file!")
            return

        for label in MtzColumnLabels:
            label_subset = [col.label for col in self.reflection_file.columns if
                            col.type == MTZColumnTypes.__getattr__(label.name).value]
            matches = list(filter(label.value.match, label_subset))
            if any(matches):
                self.__setattr__(label.name, matches[0].encode('utf-8'))

        if not any([label for label in self.summary if label is not None]):
            self.logger.error('Cannot find any column names at %s' % self.fname)
            self.error = True