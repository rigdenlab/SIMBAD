'''
Script to set up configeration options for SIMBAD

Largely borrowed from AMPLE

@author: hlasimpk\
'''

import logging
import os

from simbad.constants import SIMBAD_CONFIG_FILE
from simbad.util import version

try:
    import configparser as ConfigParser
except ImportError:
    import ConfigParser

logger = logging.getLogger(__name__)

_SECTIONS_REFERENCE = {"SIMBAD_info": ["simbad_version",
                                       "ccp4_version",
                                       "cmdline_flags"],

                       "Queueing": ["submit_array",
                                    "submit_cluster",
                                    "submit_qtype",
                                    "submit_queue"],

                       "Databases": ["PDB",
                                     "Contaminant_database",
                                     "MoRDa_database"],

                       "AMORE_options": ["SHRES",
                                         "PKLIM",
                                         "NPIC",
                                         "ROTASTEP"],

                       "Molecular_replacement": ["ncyc_contam",
                                                 "ncyc_full"],

                       "Files": ["mtz"],

                       "Unspecified": [],
                       }


class SIMBADConfigOptions(object):
    def __init__(self):

        self.d = {}

        self.cmdlin_opts = {}
        self.debug = False

        self.webserver_uri = {
            'submit_cluster': True,
            'submit_max_array': 10,
            'submit_qtype': "SGE",
            'submit_queue': "all.q",
        }

    def populate(self, cmdline_opts):

        # Convert Namespace to Dictionary
        self.cmdline_opts = cmdline_opts

        # Identify which config file to use
        config_file = self._get_config_file(cmdline_opts['config_file'])

        # Read the configuration file
        self._read_config_file(config_file)

        # Read the command line arguments
        self._read_cmdline_opts(cmdline_opts)

        # Set further options
        self._process_options()
        return

    def _get_config_file(self, cmd_file=None):
        config_file = os.path.abspath(cmd_file) if cmd_file else SIMBAD_CONFIG_FILE
        if not os.path.isfile(config_file):
            msg = "Cannot find configuration file: {0} - terminating...".format(config_file)
            logger.critical(msg)
            raise RuntimeError(msg)

        logger.debug("Using configuration file: %s", config_file)
        return config_file

    def _process_options(self):

        """
        Handle any top-level options that affect the overall running of SIMBAD.

        Notes
        ----
        Any specific processing of options should be handled in simbad/util/options_processor.py/process_options

        See also
        ----
        options_processor

        :return:
        """

        self.d['simbad_version'] = version.__version__

        if "run_dir" in self.d and not self.d["run_dir"]:
            self.d["run_dir"] = os.getcwd()

        # Set full file paths
        for k, v in self.d.iteritems():
            if k in _SECTIONS_REFERENCE["Files"] and v:
                self.d[k] = os.path.abspath(v)

        # Check is using any preset options
        if self.d['webserver_uri']: self._preset_options('webserver_uri')

        return

    def _preset_options(self, mode):
        assert hasattr(self, mode), "Unknown mode: {0}".format(mode)
        logger.info("Using preset mode: %s", mode)
        for k, v in getattr(self, mode).iteritems():
            if 'cmdline_flags' in self.d and k in self.d['cmdline_flags']:
                if self.d[k] == v:
                    msg = 'WARNING! {0} flag {1} > {2} was duplicated on the command line!'.format(mode, v, k)
                else:
                    msg = 'WARNING! Overriding {0} setting: {1} > {2} with {3}'.format(mode, k, v, self.d[k])
                logger.critical(msg)
            elif k in self.d:
                logger.debug("%s overriding default setting: %s -> %s with %s", mode, k, v, self.d[k])
                self.d[k] = v
            else:
                logger.debug("%s setting: %s -> %s", mode, k, v)
                self.d[k] = v
        return

    def _read_config_file(self, config_file):
        config = ConfigParser.SafeConfigParser()
        # We need to make sure that the keys aren't converted to lower case on reading
        config.optionxform = str
        config.read(config_file)

        for section in config.sections():

            if section not in _SECTIONS_REFERENCE:
                _SECTIONS_REFERENCE[section] = []

            # Basic switch statement to determine the type of variable
            for k, v in config.items(section):
                if v.lower() == "none":
                    self.d[k] = None

                elif v.lower == "true":
                    self.d[k] = True

                elif v.lower == "false":
                    self.d[k] = False

                elif section.lower() == "databases":
                    self.d[k] = os.path.abspath(v)

                elif section.lower() == "executables":
                    self.d[k] = os.path.abspath(v)

                elif section.lower() == "files":
                    self.d[k] = os.path.abspath(v)

                elif v.isdigit():
                    self.d[k] = int(v)

                elif self._isfloat(v):
                    self.d[k] = float(v)

                else:
                    self.d[k] = v

                _SECTIONS_REFERENCE[section].append(k)

        return

    def _read_cmdline_opts(self, cmdline_opts):
        cmdline_flags = []

        for k, v in cmdline_opts.iteritems():
            if v is not None: cmdline_flags.append(k)
            if isinstance(v, str):
                if v.lower() == "true":
                    v = True
                elif v.lower() == "false":
                    v = False
                elif v.lower() == "none":
                    v = None

            if k not in self.d:
                self.d[k] = v
            elif v is not None:
                logger.debug("Cmdline setting %s: %s -> %s", k, self.d[k], v)
                self.d[k] = v

        self.d['cmdline_flags'] = cmdline_flags
        return

    def _isfloat(self, value):
        try:
            float(value)
            return True
        except ValueError:
            return False

    def prettify_parameters(self):
        """Return the parameters nicely formated as a list of strings suitable
        for writing out to a file"""
        pstr = 'Parameters Used in this Run\n\n'
        for k, v in sorted(self.d.items()):
            pstr += "{0} : {1}\n".format(k, v)
        return pstr

    def write_config_file(self, config_file=None):
        config = ConfigParser.SafeConfigParser()
        # We need to make sure that the keys aren't converted to lower case on writing
        config.optionxform = str
        self._update_config(config)
        if config_file is None:
            # Can be None for testing
            config_file = os.path.join(self.d['work_dir'], self.d['name'] + ".ini")
        # Write config to job specific directory
        self.d["out_config_file"] = config_file
        logger.info("SIMBAD configuration written to: %s", config_file)
        with open(config_file, "w") as out: config.write(out)
        return

    def _update_config(self, config_parser):
        # Add all sections to the configparser
        for section in sorted(_SECTIONS_REFERENCE.keys()):
            if section.lower() == "no_config": continue
            config_parser.add_section(section)

        # Place all entries in our dictionary in the corresponding section in
        # the configparser
        for option in sorted(self.d.keys()):
            # Extract the section in which the entry needs to go
            sections = [k for (k, v) in _SECTIONS_REFERENCE.items() \
                        if any(entry.lower() == option.lower() for entry in v)]

            # Make sure we only have each option assigned to a single section
            section = "Unspecified" if len(sections) != 1 else sections[0]

            # We do not want to re-use files or at least not by default.
            # Comment those specifically out to avoid any errors
            if section.lower() == "no_config":
                continue
            elif section.lower() == "simbad_info" or \
                            section.lower() == "files" or \
                            section.lower() == "unspecified":
                config_parser.set(section, "#" + option, str(self.d[option]))
            else:
                config_parser.set(section, option, str(self.d[option]))

        return
