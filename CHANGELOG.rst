
Changelog
=========

[unreleased]
------------
Added
~~~~~
- ``CCP4`` container for related information
- ``standardize`` function in ``simbad.util.pdb_util`` to remove hydrogen and hetatm atoms from downloaded PDB. This resolves a bug in refmac5 where unknown ligands are bound to a pdb. 
- Test cases for ``simbad.mr.molrep_mr`` added
- ``.bat`` files so that simbad can be run in windows
- ``mtz_util.change_space_group`` function to change the space group of an mtz.

Changed
~~~~~~~
- ``ccp4_root`` function changed to ``CCP4RootDirectory`` class
- ``ccp4_version`` function changed to ``CCP4Version`` class and version extracted from official CCP4 release file
- ``-enant`` flag replaced by ``sga``
- ``simbad.mr.molrep_mr`` and ``simbad.mr.phaser_mr`` altered to check all space groups
- ``simbad.mr.molrep_mr`` modified so that if an alternative space group is found the input mtz space group will be changed accordingly. 
- ``simbad.mr.anomalous`` DANO map calculation modified and the scores reported have been changed
- ``simbad.util.mtz_util`` Altered how converted miller arrays are handled and how the R-free column label is identified
- ``simbad.lattice.latticesearch`` modified to use standardise function
- ``i2`` code updated to use ``sga``

0.1.10
------
Added
~~~~~
- ``run_tests.py`` script to execute all unittests
- PDB-redo download for structures
- Test cases for pyrvapi metadata object added
- Test case for ``latticesearch.pdb_in_results`` added
- ``-tab_prefix`` option added for JScoFe

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

0.1.0
-----
- Initial release
