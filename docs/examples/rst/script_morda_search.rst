.. _script_morda_search:

Searching the MoRDa database with SIMBAD
----------------------------------------

.. note::

   Data used throughout this example can be found in ``<ROOT>/examples/morda_example``. If SIMBAD is part of your CCP4 installation,
   then the example files can be downloaded as part of the `GitHub repository <https://github.com/rigdenlab/SIMBAD>`_.

.. warning::

   You need to have a full copy of the `MoRDa database <http://www.biomexsolutions.co.uk/morda/>`_ installed locally. We also recommend that this search is not run on local machines, but clusters instead.


0. Command line options
^^^^^^^^^^^^^^^^^^^^^^^
Check out this page explaining the :ref:`simbad-morda <simbad_morda_options>` script command line options.

1. Running the script
^^^^^^^^^^^^^^^^^^^^^
In this example, the ``simbad-morda`` script simply takes the crystallographic data file in MTZ format and the path to the local copy of the MoRDa database. It then runs the rotation function search with all models in the database against your crystallographic data.

.. literalinclude:: /../examples/morda_example/run.sh
   :language: bash
   :lines: 12-18
