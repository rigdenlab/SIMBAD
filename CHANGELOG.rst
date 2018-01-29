
Changelog
=========

[unreleased]
------------
Added
~~~~~
- ``run_tests.py`` script to execute all unittests
- PDB-redo download for structures
- Test cases for pyrvapi metadata object added
- Test case for ``latticesearch.pdb_in_results`` added
- ``-tab_prefix`` option added for JScoFe
- ``CCP4`` container for related information

Changed
~~~~~~~
- Removed reference to deprecated module ``iotbx.pdb.mining``
- Bug fix in ``simbad.lattice.latticescore`` string representation
- Bug fixes to all unittests 
- Bug fix plus added test cases for ``simbad.parsers.molrep_parser``
- Standardised parsers internal structure
- Bug fix in ``simbad.command_line.simbad_morda`` and ``simbad.command_line.simbad_full`` to fix missing ccp4i2 argument 
- Bug fix in ``simbad.lattice.latticesearch`` for duplicate entries from alternative unit cells
- Bug fix for logging and error message handling prior to logger initialisation
- Bug fix in ``simbad.util.pdb_util`` variable name 
- ``ccp4_root`` function changed to ``CCP4RootDirectory`` class
- ``ccp4_version`` function changed to ``CCP4Version`` class and version extracted from official CCP4 release file

0.1.0
-----
- Initial release
