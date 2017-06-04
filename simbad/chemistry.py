"""Storage for chemical information"""

__author__ = "Felix Simkovic & Adam Simpkin"
__date__ = "26 Apr 2017"
__version__ = "0.1"

__all__ = ['atomic_composition', 'periodic_table']


class AtomicComposition(object):

    def __init__(self):
        acids = [
            ('ALA', 'A', {'H': 5, 'C': 3, 'N': 1, 'O': 1}),
            ('CYS', 'C', {'H': 5, 'C': 3, 'N': 1, 'O': 1, 'S': 1}),
            ('ASP', 'D', {'H': 4, 'C': 4, 'N': 1, 'O': 3}),
            ('GLU', 'E', {'H': 6, 'C': 5, 'N': 1, 'O': 3}),
            ('PHE', 'F', {'H': 9, 'C': 9, 'N': 1, 'O': 1}),
            ('GLY', 'G', {'H': 3, 'C': 2, 'N': 1, 'O': 1}),
            ('HIS', 'H', {'H': 8, 'C': 6, 'N': 1, 'O': 1}),
            ('ILE', 'I', {'H': 11, 'C': 6, 'N': 1, 'O': 1}),
            ('LYS', 'K', {'H': 13, 'C': 6, 'N': 2, 'O': 1}),
            ('LEU', 'L', {'H': 11, 'C': 6, 'N': 1, 'O': 1}),
            ('MET', 'M', {'H': 9, 'C': 5, 'N': 1, 'O': 1, 'S': 1}),
            ('MSE', 'M', {'H': 9, 'C': 5, 'N': 1, 'O': 1, 'S': 1}),
            ('ASN', 'N', {'H': 6, 'C': 4, 'N': 2, 'O': 2}),
            ('PRO', 'P', {'H': 7, 'C': 5, 'N': 1, 'O': 1}),
            ('GLN', 'Q', {'H': 8, 'C': 5, 'N': 2, 'O': 2}),
            ('ARG', 'R', {'H': 13, 'C': 6, 'N': 4, 'O': 1}),
            ('SER', 'S', {'H': 5, 'C': 3, 'N': 1, 'O': 2}),
            ('THR', 'T', {'H': 7, 'C': 4, 'N': 1, 'O': 2}),
            ('VAL', 'V', {'H': 9, 'C': 5, 'N': 1, 'O': 1}),
            ('TRP', 'W', {'H': 10, 'C': 11, 'N': 2, 'O': 1}),
            ('TYR', 'Y', {'H': 9, 'C': 9, 'N': 1, 'O': 2}),
        ]
        self._aadict = {}
        for three, one, prop in acids:
            self._aadict[three] = self._aadict[one] = _AminoAcidComposition(**prop)

    def __getitem__(self, k):
        if k.upper() in self._aadict:
            return self._aadict[k.upper()]
        return None


class _AminoAcidComposition(object):

    __slots__ = ['H', 'C', 'N', 'O', 'S']

    def __init__(self, **kwargs):
        for k, v in kwargs.items():
            setattr(self, k, v)

    def __repr__(self):
        return "{0}({1})".format(
            self.__class__.__name__, ", ".join(["{0}={1}".format(k, v) for k, v in self.__dict__.items()])
        )


class PeriodicTable(object):

    def __init__(self):
        atoms = [
            ('HYDROGEN', 'H', {'atomic_number': 1, 'atomic_mass': 1.008, 'group': 'Non-metal'}),
            ('HELIUM', 'HE', {'atomic_number': 2, 'atomic_mass': 4.003, 'group': 'Nobel gas'}),
            ('LITHIUM', 'LI', {'atomic_number': 3, 'atomic_mass': 6.941, 'group': 'Alkali metal'}),
            ('BERYLLIUM', 'BE', {'atomic_number': 4, 'atomic_mass': 9.012, 'group': 'Alkaline earth'}),
            ('BORON', 'B', {'atomic_number': 5, 'atomic_mass': 10.811, 'group': 'Semi-metal'}),
            ('CARBON', 'C', {'atomic_number': 6, 'atomic_mass': 12.011, 'group': 'Non-metal'}),
            ('NITROGEN', 'N', {'atomic_number': 7, 'atomic_mass': 14.007, 'group': 'Non-metal'}),
            ('OXYGEN', 'O', {'atomic_number': 8, 'atomic_mass': 15.999, 'group': 'Non-metal'}),
            ('FLUORINE', 'F', {'atomic_number': 9, 'atomic_mass': 18.998, 'group': 'Halogen'}),
            ('NEON', 'NE', {'atomic_number': 10, 'atomic_mass': 20.180, 'group':  'Nobelgas'}),
            ('SODIUM', 'NA', {'atomic_number': 11, 'atomic_mass': 22.990, 'group': 'Alkali metal'}),
            ('MAGNESIUM', 'MG', {'atomic_number': 12, 'atomic_mass': 24.305, 'group': 'Alkaline earth'}),
            ('ALUMINIUM', 'AL', {'atomic_number': 13, 'atomic_mass': 26.982, 'group': 'Basic metal'}),
            ('SILICON', 'SI', {'atomic_number': 14, 'atomic_mass': 28.086, 'group': 'Semi-metal'}),
            ('PHOSPHORUS', 'P', {'atomic_number': 15, 'atomic_mass': 30.974, 'group': 'Non-metal'}),
            ('SULPHUR', 'S', {'atomic_number': 16, 'atomic_mass': 32.066, 'group': 'Non-metal'}),
            ('CHLORINE', 'CL', {'atomic_number': 17, 'atomic_mass': 35.453, 'group': 'Halogen'}),
            ('ARGON', 'AR', {'atomic_number': 18, 'atomic_mass': 39.948, 'group': 'Nobel gas'}),
            ('POTASSIUM', 'K', {'atomic_number': 19, 'atomic_mass': 39.098, 'group': 'Alkali metal'}),
            ('CALCIUM', 'CA', {'atomic_number': 20, 'atomic_mass': 40.078, 'group': 'Alkaline earth'}),
            ('SCANDIUM', 'SC', {'atomic_number': 21, 'atomic_mass': 44.956, 'group': 'Transition metal'}),
            ('TITANIUM', 'TI', {'atomic_number': 22, 'atomic_mass': 47.880, 'group': 'Transition metal'}),
            ('VANADIUM', 'V', {'atomic_number': 23, 'atomic_mass': 50.942, 'group': 'Transition metal'}),
            ('CHROMIUM', 'CR', {'atomic_number': 24, 'atomic_mass': 51.996, 'group': 'Transition metal'}),
            ('MANGANESE', 'MN', {'atomic_number': 25, 'atomic_mass': 54.938, 'group': 'Transition metal'}),
            ('IRON', 'FE', {'atomic_number': 26, 'atomic_mass': 55.933, 'group': 'Transition metal'}),
            ('COBALT', 'CO', {'atomic_number': 27, 'atomic_mass': 58.933, 'group': 'Transition metal'}),
            ('NICKLE', 'NI', {'atomic_number': 28, 'atomic_mass': 58.693, 'group': 'Transition metal'}),
            ('COPPER', 'CU', {'atomic_number': 29, 'atomic_mass': 63.546, 'group': 'Transition metal'}),
            ('ZINC', 'ZN', {'atomic_number': 30, 'atomic_mass': 65.390, 'group': 'Transition metal'}),
            ('GALLIUM', 'GA', {'atomic_number': 31, 'atomic_mass': 69.732, 'group': 'Basic metal'}),
            ('GERMANIUM', 'GE', {'atomic_number': 32, 'atomic_mass': 72.610, 'group': 'Semi-metal'}),
            ('ARSENIC', 'AS', {'atomic_number': 33, 'atomic_mass': 74.922, 'group': 'Semi-metal'}),
            ('SELENIUM', 'SE', {'atomic_number': 34, 'atomic_mass': 78.972, 'group': 'Non-metal'}),
            ('BROMINE', 'BR', {'atomic_number': 35, 'atomic_mass': 79.904, 'group': 'Halogen'}),
            ('KRYPTON', 'KR', {'atomic_number': 36, 'atomic_mass': 84.800, 'group': 'Nobel gas'}),
            ('RUBIDIUM', 'RB', {'atomic_number': 37, 'atomic_mass': 84.468, 'group': 'Alkali metal'}),
            ('STRONTIUM', 'SR', {'atomic_number': 38, 'atomic_mass': 87.620, 'group': 'Alkaline earth'}),
            ('YTTRIUM', 'Y', {'atomic_number': 39, 'atomic_mass': 88.906, 'group': 'Transition metal'}),
            ('ZIRCONIUM', 'ZR', {'atomic_number': 40, 'atomic_mass': 91.224, 'group': 'Transition metal'}),
            ('NIOBIUM', 'NB', {'atomic_number': 41, 'atomic_mass': 92.906, 'group': 'Transition metal'}),
            ('MOLYBDENUM', 'MO', {'atomic_number': 42, 'atomic_mass': 95.950, 'group': 'Transition metal'}),
            ('TECHNETIUM', 'TC', {'atomic_number': 43, 'atomic_mass': 98.907, 'group': 'Transition metal'}),
            ('RUTHENIUM', 'RU', {'atomic_number': 44, 'atomic_mass': 101.070, 'group': 'Transition metal'}),
            ('RHODIUM', 'RH', {'atomic_number': 45, 'atomic_mass': 102.906, 'group': 'Transition metal'}),
            ('PALLADIUM', 'PD', {'atomic_number': 46, 'atomic_mass': 106.420, 'group': 'Transition metal'}),
            ('SILVER', 'AG', {'atomic_number': 47, 'atomic_mass': 107.868, 'group': 'Transition metal'}),
            ('CADMIUM', 'CD', {'atomic_number': 48, 'atomic_mass': 112.411, 'group': 'Transition metal'}),
            ('INDIUM', 'IN', {'atomic_number': 49, 'atomic_mass': 114.818, 'group': 'Basic metal'}),
            ('TIN', 'SN', {'atomic_number': 50, 'atomic_mass': 118.710, 'group': 'Basic metal'}),
            ('ANTIMONY', 'SB', {'atomic_number': 51, 'atomic_mass': 121.760, 'group': 'Semi-metal'}),
            ('TELLURIUM', 'TE', {'atomic_number': 52, 'atomic_mass': 127.600, 'group': 'Semi-metal'}),
            ('IODINE', 'I', {'atomic_number': 53, 'atomic_mass': 126.904, 'group': 'Halogen'}),
            ('XENON', 'XE', {'atomic_number': 54, 'atomic_mass': 131.290, 'group': 'Nobel gas'}),
            ('CESIUM', 'CS', {'atomic_number': 55, 'atomic_mass': 132.905, 'group': 'Alkali metal'}),
            ('BARIUM', 'BA', {'atomic_number': 56, 'atomic_mass': 137.327, 'group': 'Alkaline earth'}),
            ('LANTHANUM', 'LA', {'atomic_number': 57, 'atomic_mass': 138.906, 'group': 'Lanthanide'}),
            ('CERIUM', 'CE', {'atomic_number': 58, 'atomic_mass': 140.115, 'group': 'Lanthanide'}),
            ('PRASEODYMIUM', 'PR', {'atomic_number': 59, 'atomic_mass': 140.908, 'group': 'Lanthanide'}),
            ('NEODYMIUM', 'ND', {'atomic_number': 60, 'atomic_mass': 144.240, 'group': 'Lanthanide'}),
            ('PROMETHIUM', 'PM', {'atomic_number': 61, 'atomic_mass': 144.913, 'group': 'Lanthanide'}),
            ('SAMARIUM', 'SM', {'atomic_number': 62, 'atomic_mass': 150.360, 'group': 'Lanthanide'}),
            ('EUROPIUM', 'EU', {'atomic_number': 63, 'atomic_mass': 151.966, 'group': 'Lanthanide'}),
            ('GADOLINIUM', 'GD', {'atomic_number': 64, 'atomic_mass': 157.250, 'group': 'Lanthanide'}),
            ('TERBIUM', 'TB', {'atomic_number': 65, 'atomic_mass': 158.925, 'group': 'Lanthanide'}),
            ('DYSPROSIUM', 'DY', {'atomic_number': 66, 'atomic_mass': 162.50, 'group': 'Lanthanide'}),
            ('HOLMIUM', 'HO', {'atomic_number': 67, 'atomic_mass': 164.930, 'group': 'Lanthanide'}),
            ('ERBIUM', 'ER', {'atomic_number': 68, 'atomic_mass': 167.260, 'group': 'Lanthanide'}),
            ('THULIUM', 'TM', {'atomic_number': 69, 'atomic_mass': 168.934, 'group': 'Lanthanide'}),
            ('YTTERBIUM', 'YB', {'atomic_number': 70, 'atomic_mass': 173.04, 'group': 'Lanthanide'}),
            ('LUTETIUM', 'LU', {'atomic_number': 71, 'atomic_mass': 174.967, 'group': 'Lanthanide'}),
            ('HAFNIUM', 'HF', {'atomic_number': 72, 'atomic_mass': 178.490, 'group': 'Transition metal'}),
            ('TANTALUM', 'TA', {'atomic_number': 73, 'atomic_mass': 180.948, 'group': 'Transition metal'}),
            ('TUNGSTEN', 'W', {'atomic_number': 74, 'atomic_mass': 183.850, 'group': 'Transition metal'}),
            ('RHENIUM', 'RE', {'atomic_number': 75, 'atomic_mass': 186.207, 'group': 'Transition metal'}),
            ('OSMIUM', 'OS', {'atomic_number': 76, 'atomic_mass': 190.230, 'group': 'Transition metal'}),
            ('IRIDIUM', 'IR', {'atomic_number': 77, 'atomic_mass': 192.220, 'group': 'Transition metal'}),
            ('PLATINUM', 'PT', {'atomic_number': 78, 'atomic_mass': 195.08, 'group': 'Transition metal'}),
            ('GOLD', 'AU', {'atomic_number': 79, 'atomic_mass': 196.967, 'group': 'Transition metal'}),
            ('MERCURY', 'HG', {'atomic_number': 80, 'atomic_mass': 200.590, 'group': 'Transition metal'}),
            ('THALLIUM', 'TL', {'atomic_number': 81, 'atomic_mass': 204.383, 'group': 'Basic metal'}),
            ('LEAD', 'PB', {'atomic_number': 82, 'atomic_mass': 207.200, 'group': 'Basic metal'}),
            ('BISMUTH', 'BI', {'atomic_number': 83, 'atomic_mass': 208.980, 'group': 'Basic metal'}),
            ('POLONIUM', 'PO', {'atomic_number': 84, 'atomic_mass': 208.982, 'group': 'Semi-metal'}),
            ('ASTATINE', 'AT', {'atomic_number': 85, 'atomic_mass': 209.987, 'group': 'Halogen'}),
            ('RADON', 'RN', {'atomic_number': 86, 'atomic_mass': 222.018, 'group': 'Nobel gas'}),
            ('FRANCIUM', 'FR', {'atomic_number': 87, 'atomic_mass': 223.020, 'group': 'Alkali metal'}),
            ('RADIUM', 'RA', {'atomic_number': 88, 'atomic_mass': 226.025, 'group': 'Alkaline earth'}),
            ('ACTINIUM', 'AC', {'atomic_number': 89, 'atomic_mass': 227.028, 'group': 'Actinide'}),
            ('THORIUM', 'TH', {'atomic_number': 90, 'atomic_mass': 232.038, 'group': 'Actinide'}),
            ('PROTACTINIUM', 'PA', {'atomic_number': 91, 'atomic_mass': 231.036, 'group': 'Actinide'}),
            ('URANIUM', 'U', {'atomic_number': 92, 'atomic_mass': 238.029, 'group': 'Actinide'}),
            ('NEPTUNIUM', 'NP', {'atomic_number': 93, 'atomic_mass': 237.048, 'group': 'Actinide'}),
            ('PLUTONIUM', 'PU', {'atomic_number': 94, 'atomic_mass': 244.064, 'group': 'Actinide'}),
            ('AMERICIUM', 'AM', {'atomic_number': 95, 'atomic_mass': 243.061, 'group': 'Actinide'}),
            ('CURIUM', 'CM', {'atomic_number': 96, 'atomic_mass': 247.070, 'group': 'Actinide'}),
            ('BERKELIUM', 'BK', {'atomic_number': 97, 'atomic_mass': 247.070, 'group': 'Actinide'}),
            ('CALIFORNIUM', 'CF', {'atomic_number': 98, 'atomic_mass': 251.080, 'group': 'Actinide'}),
            ('EINSTEINIUM', 'ES', {'atomic_number': 99, 'atomic_mass': 254.000, 'group': 'Actinide'}),
            ('FERMIUM', 'FM', {'atomic_number': 100, 'atomic_mass': 257.095, 'group': 'Actinide'}),
            ('MENDELEVIUM', 'MD', {'atomic_number': 101, 'atomic_mass': 258.100, 'group': 'Actinide'}),
            ('NOBELIUM', 'NO', {'atomic_number': 102, 'atomic_mass': 259.101, 'group': 'Actinide'}),
            ('LAWRENCIUM', 'LR', {'atomic_number': 103, 'atomic_mass': 262.00, 'group': 'Actinide'}),
            ('RUTHERFORDIUM', 'RF', {'atomic_number': 104, 'atomic_mass': 261.000, 'group': 'Transition metal'}),
            ('DUBNIUM', 'DB', {'atomic_number': 105, 'atomic_mass': 262.000, 'group': 'Transition metal'}),
            ('SEABORGIUM', 'SG', {'atomic_number': 106, 'atomic_mass': 266.000, 'group': 'Transition metal'}),
            ('BOHRIUM', 'BH', {'atomic_number': 107, 'atomic_mass': 264.000, 'group': 'Transition metal'}),
            ('HASSIUM', 'HS', {'atomic_number': 108, 'atomic_mass': 269.000, 'group': 'Transition metal'}),
            ('MEITNERIUM', 'MT', {'atomic_number': 109, 'atomic_mass': 268.000, 'group': 'Transition metal'}),
            ('DARMSTADTIUM', 'DS', {'atomic_number': 110, 'atomic_mass': 269.000, 'group': 'Transition metal'}),
            ('ROENTGENIUM', 'RG', {'atomic_number': 111, 'atomic_mass': 272.000, 'group': 'Transition metal'}),
            ('COPERNICIUM', 'CN', {'atomic_number': 112, 'atomic_mass': 277.000, 'group': 'Transition metal'}),
            ('UNUNTRIUM', 'UUT', {'atomic_number': 113, 'atomic_mass': None, 'group': 'Basic metal'}),
            ('FLEROVIUM', 'FL', {'atomic_number': 114, 'atomic_mass': 289.000, 'group': 'Basic metal'}),
            ('UNUNPENTIUM', 'UUP', {'atomic_number': 115, 'atomic_mass': None, 'group': 'Basic metal'}),
            ('LIVERMORIUM', 'LV', {'atomic_number': 116, 'atomic_mass': 298.000, 'group': 'Basic metal'}),
            ('UNUNSEPTIUM', 'UUS', {'atomic_number': 117, 'atomic_mass': None, 'group': 'Halogen'}),
            ('UNUNOCTIUM', 'UUO', {'atomic_number': 118, 'atomic_mass': None, 'group': 'Nobel gas'})
        ]

        self._atomdict = {}
        for name, code, prop in atoms:
            self._atomdict[name] = self._atomdict[code] = _AtomComposition(**prop)

    def __getitem__(self, k):
        if k.upper() in self._atomdict:
            return self._atomdict[k.upper()]
        return None


class _AtomComposition(object):

    __slots__ = ['atomic_number', 'atomic_mass', 'group']

    def __init__(self, **kwargs):
        for k, v in kwargs.items():
            setattr(self, k, v)

    def __repr__(self):
        return "{0}({1})".format(
            self.__class__.__name__, ", ".join(["{0}={1}".format(k, v) for k, v in self.__dict__.items()])
        )

# Instantiate some stuff here so we can call it immediately
atomic_composition = AtomicComposition()
periodic_table = PeriodicTable()