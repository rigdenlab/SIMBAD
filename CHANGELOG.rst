
Changelog
=========

[unreleased]
------------

Added
~~~~~
- When using phaser, a LLG/TFZ > 120/8 was added as an additional criteria for early termination
- Added a process all flag to amore/phaser rotation functions so that they skip the early termination funciton

Changed
~~~~~~~
- Default sequence identity changed from 30 percent to 70
- Fixed a test case that was broken
- Fixed several bugs when running jobs as multiple chunks
- Changed process to skip successive chunks if a solution is found

0.1.16
------

Added
~~~~~
- Added link to CCP4 SW 2019 presentation

Changed
~~~~~~~
- Changed Phaser settings so it doesn't check ensemble deviation. We want to run all the ensembles even if they're poor. 
- Fixed bugs with how ctruncate was being called in mrbump
- Refactored the code to assign mtz labels and fixed a bug caused by changes to the latest CCTBX version  

0.1.15
------

Added
~~~~~
- Added code to Pyrvapi to return metadata when performing a lattice search with a space group and unit cell

Changed
~~~~~~~
- Changed phaser to just look for 1 molecule in the ASU in full MR search due to time constraints. 
- Fixed a bug that was providing the wrong number of processors to PyJob
- Fixed bug in the simbad-database code

0.1.14
------

Added
~~~~~
- Added code to Pyrvapi to allow SIMBAD to work on CCP4online
- Added reference manager code
- Added bibtex file containing SIMBAD references
- Added citation tab to pyvapi

Changed
~~~~~~~
- Changed `-morda_db` argument in `simbad.command_line.simbad-database` and `simbad.command_line` so that it now has a default location and so that after generating the database the `-morda_db` flag doesn't need to be called every time. 
- Fixed a bug in how pointless was called that was changing the a,b,c,alpha,beta,gamma order. 
- tidied up clean up function

0.1.13
------

Added
~~~~~
- Added in code that uses pointless to reindex mtz files
- EXE_EXT calls to all executable commands to allow for windows compatibility
- Completely reformatted all pyjob calls to use version 0.2.0 and updated dependencies list.
- Function to get the sequence from a pdb file
- Function to generate ensemble database


Changed
~~~~~~~
- Replaced CCTBX code that checked the columns in the input MTZ with MrBump code as CCTBX was giving errors for certain input MTZs. 
- Replaced `simbad.util.mtz_util.get_labels` with a class called `simbad.util.mtz_util.GetLabels`. This returns more types of input labels and simplifies how the labels are passed into other functions.
- Updated test cases affected by changes
- `simbad.util.mtz_util.GetLabels` was changed to use miller arrays and is therefore more robust when handling columns with non-standard names
- `simbad.rotsearch.phaser_search` changed to rank models by RFZ rather than LLG.

0.1.12
------

Added
~~~~~
- Check for ``SSL_CERT_FILE`` environmental variable in all command line scripts
- eRMSD calculation added into the phaser rotation search, default set to 70% ID for now but this may change
- Added SIMBAD paper to docs
- Added output_files directory to store all output files

Changed
~~~~~~~
- Fixed bug in phaser that fixed a problem in check all/enant spacegroups
- Changed molrep/phaser to output an hkl, this is needed for changes of basis as a result of all/enant searches
- Altered the cleanup algorithm to completely remove the mr_search directory and instead rely on the newly created output_files directory
- Altered logging to use enum, closes #81
- Removed eRMSD calculation and changed to use the seqenence identity directly through phaser, therefore using their equation directly. 
- Altered phaser so that it outputs the input MTZ with a basis change instead of the phaser output mtz, due to missing r-free columns 
- Fixed bug in rotation search solution check due to missing argument

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
