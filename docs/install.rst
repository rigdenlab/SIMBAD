.. _installation:

Installation of SIMBAD
======================

.. code-block:: bash
   
   git clone https://github.com/rigdenlab/SIMBAD.git
   cd SIMBAD
   ccp4-python setup.py build --script-python-path ccp4-python install --install-scripts $CCP4/bin --install-lib $CCP4/lib/py2/site-packages

This will install SIMBAD into your CCP4 installation. On top of installing the source code, an executable script ``simbad`` should be automatically installed.

----

Installation of SIMBAD databases
================================

This section of the documentation relates to the installation or updating of SIMBAD-specific databases. We advise you update your databases regularly because your chances of success increase with every newly deposited structure.

SIMBAD currently requires four databases, although only three can be updated/installed manually. The contaminant database is "as is" and for updated versions, please contact us directly at ``ccp4[at]stfc.ac.uk``. However, do not that the database is automatically updated with SIMBAD updates in CCP4, i.e. manual updates will probably not be necessary.

The remaining databases required by SIMBAD are
    - :ref:`Lattice parameter database`
    - :ref:`MoRDa-like database with domain coordinates`
    - :ref:`Custom database`

In the following sections, we will explain how to install/update each of these databases. Note, all databases require an active internet connection!

Lattice parameter database
~~~~~~~~~~~~~~~~~~~~~~~~~~

The lattice parameter database ships by default with SIMBAD. However, you might want to update this database regularly.

To do so, use the ``simbad-create-db`` command with the ``lattice`` subcommand. If you are in a Unix terminal, use the following code to update:

.. code-block:: bash

   $ simbad-create-db lattice

Hit the ``<Enter>`` key and your default database will be updated automatically.

If you do not have write permissions to the CCP4 installation or would prefer to keep a separate copy of the update lattice parameter database, you can use the ``-latt_db`` flag with a path to your preferred location. For example, your command could instead look like this:

.. code-block:: bash

   $ simbad-create-db lattice -latt_db $HOME/Documents/simbad_lattice_db.npz

.. note::
   If you create a custom copy of the lattice parameter database, make sure to point SIMBAD to that. Provide the ``-latt_db`` flag when invoking relevant scripts.

MoRDa-like database with domain coordinates
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The MoRDa-like database in SIMBAD is a requirement for running the full rotation function search. It does not ship by default and its creation relies on the `MoRDa <http://www.biomexsolutions.co.uk/morda/>`_ database. It uses the curated single domain files and relevant scripts to generate its own copy. **As a result, this database and associated features in SIMBAD are currently limited to Unix systems.**

The most basic command to install this database is

.. code-block:: bash

   $ simbad-create-db morda $HOME/Documents/simbad_db

The previous command will download, and install the SIMBAD database to the directory ``simbad_db``. It will create a PDB-like substructure of folders to group relevant entries. **Each file is encoded to save disk space, for instructions on how to decode it, please contact us at ccp4[at]stfc.ac.uk``.**

The creation of the MoRDa-like SIMBAD database scales with the number of processors available. If you have more available, provide the ``-nproc`` command line argument and we will make use of as many processors as you provide. **If you are installing SIMBAD on a computing cluster, make use of the ``-submit_cluster`` option.**

.. code-block:: bash

   $ simbad-create-db morda -nproc 10 $HOME/Documents/simbad_db

This database will currently require ~3Gb of disk space. If you do not have much more available, you might want to consider providing the ``-chunk_size`` argument to the script call. By default, this value is ``5000`` meaning that 5000 domains are processed at the same time. However, this will require ~100Gb of temporary disk space to be available. If you do not have this space available, reduce this number accordingly [``-chunk_size 100`` does not exceed ~10Gb].

.. code-block:: bash

   $ simbad-create-db morda -chunk_size 100 $HOME/Documents/simbad_db

After the first installation of this database, we do not need to process every domain again in consecutive runs. If you want to update your database in the future, you can simply run the same commands as before, and point the script to the root of the simbad database. I.e., if we created the database with the command ``simbad-create-db morda -chunk_size 100 $HOME/Documents/simbad_db``, we can update it now with the following:

.. code-block:: bash

   $ simbad-create-db morda -chunk_size 100 $HOME/Documents/simbad_db

The installation procedure will determine any new files in the MoRDa database, and only process them.

Custom database
~~~~~~~~~~~~~~~

Alternatively you may wish to run SIMBAD using a custom database. In order to do this, first the database must be converted into a SIMBAD compatible format.

SIMBAD provides a script to create a SIMBAD compatible database from a database of PDB files. The command to run this is:

.. code-block:: bash

   $ simbad-create-db custom -input_db $HOME/Documents/input_db -custom_db $HOME/Documents/custom_db

.. note::
   If you create a custom database, make sure to point SIMBAD to that. Provide the ``-contaminant_db`` or ``-morda_db`` flags when invoking relevant scripts.