.. _script_lattice_search:

Performing a lattice search with SIMBAD
---------------------------------------

.. note::

   Data used throughout this example can be found in ``<ROOT>/examples/contaminant_example``. If SIMBAD is part of your CCP4 installation,
   then the example files can be downloaded as part of the `GitHub reposity <https://github.com/rigdenlab/SIMBAD>`_.


0. Command line options
^^^^^^^^^^^^^^^^^^^^^^^
Check out this page explaining the :ref:`simbad-lattice <simbad_lattice_options>` script command line options.

1. Running the script
^^^^^^^^^^^^^^^^^^^^^
In this example, the ``simbad-lattice`` script simply takes the crystallographic data file in MTZ format, and runs the lattice search on your local machine.

.. literalinclude:: /../examples/lattice_example/run.sh
   :language: bash
   :lines: 10-11 

