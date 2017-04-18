.. _installation:

Installation
============

.. code-block:: bash
   
   git clone https://github.com/rigdenlab/SIMBAD.git
   cd SIMBAD
   ccp4-python setup.py build --script-python-path ccp4-python install --install-scripts $CCP4/bin --install-lib $CCP4/lib/py2/site-packages

This will install SIMBAD into your CCP4 installation. On top of installing the source code, an executable script ``simbad`` should be automatically installed.

SIMBAD uses a modified version of AMORE to perform rotation searches. This version of amore will be added to the CCP4
distribution soon however in the meantime you will need to manually link this to your CCP4 installation:

.. code-block:: bash

   ln -s SIMBAD/static/amoreCCB2.exe $CCP4/bin/amoreCCB2.exe

