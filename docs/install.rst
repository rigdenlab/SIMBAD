.. _installation:

Installation of SIMBAD
======================

SIMBAD is officially shipped with CCP4 v7.0.46.

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

To do so, use the ``simbad-database`` command with the ``lattice`` subcommand. If you are in a Unix terminal, use the following code to update:

.. code-block:: bash

   $ simbad-database lattice

Hit the ``<Enter>`` key and your default database will be updated automatically.

If you do not have write permissions to the CCP4 installation or would prefer to keep a separate copy of the update lattice parameter database, you can use the ``-latt_db`` flag with a path to your preferred location. For example, your command could instead look like this:

.. code-block:: bash

   $ simbad-database lattice -latt_db $HOME/Documents/simbad_lattice_db.npz

.. note::
   If you create a custom copy of the lattice parameter database, make sure to point SIMBAD to that. Provide the ``-latt_db`` flag when invoking relevant scripts.

Contaminant database
~~~~~~~~~~~~~~~~~~~~

The contaminant database is shipped with SIMBAD by default. A script to update the contaminant database will be made available in a future update.

MoRDa-like database with domain coordinates
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The MoRDa-like database in SIMBAD is a requirement for running the MoRDa database search. As the database is quite large (~3Gb) it does not ship with SIMBAD by default.
The `MoRDa <http://www.biomexsolutions.co.uk/morda/>`_ database is derived from the PDB and contains a compact description of non-redundant protein chains, domains, homo- and hetero-oligomers. Instructions to install the MoRDa through CCP4 are available `here <http://www.ccp4.ac.uk/html/morda_installation.html>`_. In order to use the MoRDa database in SIMBAD, the database must be reformatted, thus creating our MoRDa-like database.
**MoRDa is not currently available on Windows, therefore this database and associated features in SIMBAD are currently limited to Unix systems.**

The most basic command to generate the MoRDa-like database is:

.. code-block:: bash

   $ simbad-database morda $HOME/Documents/simbad_db

The previous command will install the MoRDa-like database to the directory ``simbad_db``. It will create a PDB-like substructure of folders to group relevant entries. **Each file is encoded to save disk space, for instructions on how to decode it, please contact us at ccp4[at]stfc.ac.uk``.**

If MoRDa is installed, SIMBAD will use the associated MoRDa database to generate the MoRDa-like database. Otherwise, the MoRDa package will be temporarily downloaded.

The creation of the MoRDa-like SIMBAD database scales with the number of processors available. If you have more available, provide the ``-nproc`` command line argument and we will make use of as many processors as you provide. **If you are installing SIMBAD on a computing cluster, make use of the ``-submit_cluster`` option.**

.. code-block:: bash

   $ simbad-database morda -nproc 10 $HOME/Documents/simbad_db

After the first installation of this database, we do not need to process every domain again in consecutive runs. If you want to update your database in the future, you can simply run the same commands as before, and point the script to the root of the simbad database. I.e., if we created the database with the command ``simbad-database morda $HOME/Documents/simbad_db``, we can update it now with the following:

.. code-block:: bash

   $ simbad-database morda $HOME/Documents/simbad_db

The installation procedure will determine any new files in the MoRDa database, and only process them.

Custom database
~~~~~~~~~~~~~~~

Alternatively you may wish to run SIMBAD using a custom database. In order to do this, first the database must be converted into a SIMBAD compatible format.

SIMBAD provides a script to create a SIMBAD compatible database from a database of PDB files. The command to run this is:

.. code-block:: bash

   $ simbad-database custom $HOME/Documents/custom_db $HOME/Documents/input_db

.. note::
   If you create a custom database, make sure to point SIMBAD to that. Provide the ``-cont_db`` or ``-morda_db`` flags when invoking relevant scripts.
