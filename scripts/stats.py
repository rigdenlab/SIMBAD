#!/usr/bin/env python

"""Script to compute some SIMBAD statistics --- work in progress"""
from __future__ import print_function

__author__ = "Felix Simkovic & Adam Simpkin"
__date__ = "21 Jul 2019"
__version__ = "2.0"

import json
import logging
import os
import re
import requests

LOG = logging.getLogger(__name__)
URL = 'https://search.rcsb.org/rcsbsearch/v2/query?json={"query":{"type":"terminal","service":"full_text","parameters":{"value":"simbad"}},"return_type":"entry"}'


def update_readme(n_sols):
    readme = os.path.join(os.path.dirname(__file__), "..", "README.rst")
    with open(readme, "r") as fh:
        text = fh.read()
    text = re.sub(
        r"solution\%20count-\d+-blue", "solution%20count-{}-blue".format(n_sols), text
    )
    with open(readme, "w") as fh:
        fh.write(text)


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)

    response = requests.get(URL)
    out = json.loads(response.text)

    n_solutions = out['total_count']
    LOG.info("Found %s solutions", n_solutions)
    update_readme(n_solutions)
