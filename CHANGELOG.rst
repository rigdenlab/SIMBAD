
Changelog
=========

[unreleased]
------------

Added
~~~~~
- Check for ``SSL_CERT_FILE`` environmental variable in all command line scripts
- eRMSD calculation added into the phaser rotation search, default set to 70% ID for now but this may change

Changed
~~~~~~~

0.1.11
------

Added
~~~~~
- ``CCP4`` container for related information
- ``standardize`` function in ``simbad.util.pdb_util`` to remove hydrogen and hetatm atoms from downloaded PDB. This resolves a bug in refmac5 where unknown ligands are bound to a pdb. 
- Test cases for ``simbad.mr.molrep_mr`` added
- ``.bat`` files so that simbad can be run in windows
- Fix for bug calling multiple programs from a single script in Windows.
- ``mtz_util.change_space_group`` function to change the space group of an mtz.
- Test case for ``simbad.util.pdb_util`` added
- Added in a function to check if there is a solution in the rotation search given a sufficiently high peak height
- Perform a cell content analysis prior to the AMORE search in order to rank search models by molecular weight
- Phaser rotation search module 
- ``parsers.anode_parser`` added

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
- Updated lattice and mtz_util test cases
- Altered default MR program used in i2 to molrep
- ``simbad.mr.phaser_mr`` altered to use phaser python interface
- Reduced the number of refinement cycles for the lattice search
- Reduced the max penalty score in the lattice search from 12 to 7 to speed up the search
- Reduce the max lattice results from 50 to 20 to speed up the search
- Fixed bug when standardising files in the lattice search
- Updated ccp4i2 files to reflect recent changes made to ccp4i2
- ``simbad.rotsearch.amore_search`` moved to ``simbad.rotsearch.__init__.py`` in addition to phaser module
- Fixed test cases and parsers affected by change to rotation search code
- Altered anomalous fourier calulcation to use ANODE
- Refactored the rotsearch module and the scoring classes

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
