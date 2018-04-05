#!/usr/bin/env ccp4-python

"""Module containing a framework for unittesting of SIMBAD modules"""

from __future__ import print_function

__author__ = "Felix Simkovic"
__date__ = "22 Mar 2016"
__version__ = "1.0"

import glob
import logging
import os
import sys

from unittest import TestLoader, TextTestRunner, TestSuite

SIMBAD_DIR = os.path.join(os.path.dirname(__file__),  "..", "simbad")
PACKAGES = ["db", "lattice", "mr", "parsers", "rotsearch", "util"]


def get_cli_args():
    import argparse
    parser = argparse.ArgumentParser(prog="run_tests.py")
    parser.add_argument(
        '-b', dest='buffer', action="store_false", default=True, help="debugging by printing print messages")
    parser.add_argument('-v', dest="verbosity", default=2, type=int, help="level of verbosity [default: 2]")
    parser.add_argument('test_cases', nargs='*', help="[ {0} ]".format(" | ".join(PACKAGES)))
    return parser.parse_args()


class SIMBADUnittestFramework(object):
    """Framework to run SIMBAD unittesting"""

    def run(self, buffer=False, cases=None, pattern="test*.py", verbosity=2):
        """Main routine for running the test cases"""
        suite = SuiteLoader().load_suite(SIMBAD_DIR, cases=cases, pattern=pattern)
        if int(suite.countTestCases()) <= 0:
            msg = 'Could not find any tests to run in directory: {0}'.format(SIMBAD_DIR) + os.linesep
            sys.stderr.write(msg)
            sys.exit(1)
        logging.disable(logging.CRITICAL)
        TextTestRunner(verbosity=verbosity, buffer=buffer).run(suite)
        logging.disable(logging.NOTSET)


class SuiteLoader(object):
    """Loader designed to obtain all test cases in a package"""

    def load_suite(self, directory, pattern="test*.py", cases=[]):
        """Load a unittest test suite"""
        if len(cases) < 1:
            search_pattern = os.path.join(directory, "*")
            cases = [
                os.path.basename(folder) for folder in glob.iglob(search_pattern)
                if os.path.isdir(folder) and os.path.basename(folder) in PACKAGES
            ]
        return self._load_suite(cases, pattern, directory)

    def _load_suite(self, cases, pattern, directory):
        suite = TestSuite()
        for case in cases:
            path = os.path.join(directory, case, "tests")
            try:
                _suite = TestLoader().discover(path, pattern=pattern, top_level_dir=directory)
                suite.addTests(_suite)
                del _suite
            except ImportError:
                print("*** not a package: {0} ***".format(path))
        return suite


if __name__ == "__main__":
    args = get_cli_args()
    m = SIMBADUnittestFramework()
    m.run(buffer=args.buffer, cases=args.test_cases, verbosity=args.verbosity)
