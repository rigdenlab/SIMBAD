**********************************************************************
Sequence Independent Molecular Replacement Based on Available Database
**********************************************************************

Installation
++++++++++++

.. code:: bash
   
   git clone https://github.com/rigdenlab/SIMBAD.git
   cd SIMBAD
   ccp4-python setup.py build --script-python-path ccp4-python install --install-scripts $CCP4/bin --install-lib $CCP4/lib/py2/site-packages

This will install SIMBAD into your CCP4 installation. On top of installing the source code, an executable script ``simbad`` should be automatically installed.

SIMBAD uses a modified version of AMORE to perform rotation searches. This version of amore will be added to the CCP4
distribution soon however in the meantime you will need to manually link this to your CCP4 installation:

.. code:: bash

   ln -s SIMBAD/static/amoreCCB2.exe $CCP4/bin/amoreCCB2.exe

Documentation & Usage
+++++++++++++++++++++
Please refer to `SIMBAD's documentation <http://simbad.readthedocs.io/en/latest/>`_

Contributing
++++++++++++
There are two ways by which you can contribute to SIMBAD:

1. Submit any suggestions to the `GitHub Issue Tracker <https://github.com/rigdenlab/simbad/issues>`_, or
2. Fork this repository, commit your changes and submit a pull request.


Contributors
++++++++++++

- Adam Simpkin
- Ronan Keegan
- Felix Simkovic
- Daniel Rigden
