#!/usr/bin/env python

"""Script to compute some SIMBAD statistics --- work in progress"""
from __future__ import print_function

__author__ = "Felix Simkovic"
__date__ = "21 Jul 2019"
__version__ = "2.0"

import logging

logging.basicConfig(level=logging.INFO)
import os
import re
import requests
import urllib2

LOG = logging.getLogger(__name__)
BASE_URL = "http://www.rcsb.org/pdb/rest/"
PAYLOAD = """
<orgPdbQuery>
<queryType>org.pdb.query.simple.AdvancedKeywordQuery</queryType>
<description>Text Search for: simbad</description>
<keywords>simbad</keywords>
</orgPdbQuery>
"""


def update_readme(n_sols):
    readme = os.path.join(os.path.dirname(__file__), "..", "README.rst")
    with open(readme, "r") as fh:
        text = fh.read()
    text = re.sub(r"solution\%20count-\d+-blue", "solution%20count-{}-blue".format(n_sols), text)
    with open(readme, "w") as fh:
        fh.write(text)


if __name__ == "__main__":
    header = {"Content-Type": "application/x-www-form-urlencoded"}
    url = requests.compat.urljoin(BASE_URL, "search")

    response = requests.post(url, data=PAYLOAD, headers=header)
    response.raise_for_status()

    entries = response.content.splitlines()

    n_solutions = len(entries)
    LOG.info("Found %s solutions", n_solutions)
    update_readme(n_solutions)
