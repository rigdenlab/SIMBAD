#!/usr/bin/env ccp4-python

"""Script to compute some SIMBAD statistics --- work in progress"""

from __future__ import print_function

__author__ = "Felix Simkovic"
__date__ = "30 Mar 2018"
__version__ = "1.0"

import os
import re
import requests

url = 'http://www.rcsb.org/pdb/rest/search'
query_text = """
<orgPdbQuery>
<queryType>org.pdb.query.simple.AdvancedKeywordQuery</queryType>
<description>Text Search for: simbad</description>
<keywords>simbad</keywords>
</orgPdbQuery>
"""
header = {'Content-Type': 'application/x-www-form-urlencoded'}
response = requests.post(url, data=query_text, headers=header)
if response.status_code == 200:
    count = len(response.text.split())
    print("SIMBAD solved {} structures".format(count))
else:
    raise Exception("Could not retrieve data")
readme = os.path.join(os.path.dirname(__file__), "..", "README.rst")
with open(readme, "r") as f_in:
    data = f_in.read()
    data = re.sub(r"solution\%20count-\d+-blue", "solution%20count-{}-blue".format(count), data)
with open(readme, "w") as f_out:
    f_out.write(data)
