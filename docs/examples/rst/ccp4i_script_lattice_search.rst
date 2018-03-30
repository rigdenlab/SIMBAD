.. _ccp4i_script_lattice_search:

Performing a lattice search with SIMBAD
---------------------------------------

.. note::

   Data used throughout this example can be found in ``<ROOT>/examples/lattice_example``. If SIMBAD is part of your CCP4 installation,
   then the example files can be downloaded as part of the `GitHub repository <https://github.com/rigdenlab/SIMBAD>`_.


0. Command line options
^^^^^^^^^^^^^^^^^^^^^^^
Check out this page explaining the :ref:`simbad-lattice <simbad_lattice_options>` script command line options.

1. Running the SIMBAD lattice search
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
The lattice parameters for a crystallised protein are often unique and therefore provide a quick and easy route to identify previously solved structures.
The SIMBAD lattice search compares the lattice parameters of an input MTZ file to all the structures in the PDB.

In this example, the ``simbad-lattice`` script simply takes the crystallographic data file in MTZ format, and runs the lattice search followed by Molecular Replacement on your local machine.

SIMBAD can be found under the Molecular Replacement menu in the CCP4i GUI:

.. figure:: ../images/ccp4i_simbad.png
   :width: 50%
   :align: center

Opening SIMBAD will bring you to the following menu:

.. figure:: ../images/ccp4i_simbad_run.png
   :width: 50%
   :align: center

To run the lattice search, all that is needed is an MTZ file:

.. figure:: ../images/ccp4i_simbad_lattice.png
   :width: 50%
   :align: center

SIMBAD Output
-------------
On starting SIMBAD a separate window will appear summarising the progress of the SIMBAD lattice search and any results found.
The window will contain up to three tabs, the contents of which are explained below:

.. contents:: Output Tabs
   :depth: 1
   :local:

Log File
^^^^^^^^
This displays the text output by SIMBAD as it is running. Any problems or errors will be displayed here.

.. figure:: ../images/ccp4i_lattice_log.png
   :width: 50%
   :align: center

------------------------------------------------------------------


Lattice Parameter Search Results
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
The Lattice Parameter Search Results tab contains 4 different sections. Below you can find information about each:

.. contents:: Sections
   :depth: 1
   :local:

.. figure:: ../images/ccp4i_lattice_results.png
   :width: 50%
   :align: center

------------------------------------------------------------------

Lattice Parameter Search Results
================================
This shows the results from the Lattice Parameter Search. The columns of the table are:

* **PDB_code:** The 4 letter code representing the protein in the protein data bank
* **alt:** Alternative Niggli cell, denoted by a *
* **a:** Lattice parameter a
* **b:** Lattice parameter b
* **c:** Lattice parameter c
* **alpha:** Lattice parameter alpha
* **beta:** Lattice parameter beta
* **gamma:** Lattice parameter gamma
* **length_penalty:** The sum of the differences between lattice parameters a, b and c for the model and the target
* **angle_penalty:** The sum of the differences between lattice parameters alpha, beta and gamma for the model and the target
* **total_penalty:** The sum of the length penalty and the angle penalty
* **Probability_score:** The probability that a structure giving a total penalty score will provide a solution

The structures are scored by total_penalty score where a lower score is better.

Molecular Replacement Search Results
====================================
Molecular replacement is performed on the top 20 structures identified by the Lattice Parameter Search. This section displays the results of that molecular replacement.

By default SIMBAD runs Molecular replacement using MOLREP. If run the following columns are added to the table:

* **molrep_score:** MOLREP score for the Molecular Replacement solution
* **molrep_tfscore:** MOLREP translation function score for the Molecular Replacement solution

Alternatively SIMBAD can run Molecular replacement using PHASER. If run the following columns are added to the table:

* **phaser_llg:** PHASER Log-likelihood gain for the Molecular Replacement solution
* **phaser_tfz:** PHASER Translation Function Z-score for the Molecular Replacement solution
* **phaser_rfz:** PHASER Rotational Function Z-score for the Molecular Replacement solution

Following Molecular replacement, refinement is run using REFMAC. This add the following columns are added to the table:

* **final_r_fact:** R-fact score for REFMAC refinement of the Molecular Replacement solution
* **final_r_free:** R-free score for REFMAC refinement of the Molecular Replacement solution

.. note::

   Typically a result with a final_r_fact and a final_r_free below 0.45 is indicative of a solution.

Additionally if there is anomalous signal in your data set SIMBAD will try to validate the quality of the molecular replacement solution using by plotting the peaks from an anomalous fourier map. If run the following columns are added to the table:

* **dano_peak_height:** The highest anomalous peaks found
* **dano_z_score:** DANO peak Z-score

Top 10 Lattice Parameter Search Downloads
=========================================
This section contains the refined placed model and mtz for the top 10 solutions in the lattice parameter search (as ranked by final_r_free)

.. note::

   This section may contain less than 10 solutions if for example a solution was found and the molecular replacement search ended early

Top 10 Lattice Parameter Search Log Files
=========================================
This section contains the molecular replacement and refinement logs for the top 10 solutions in the lattice parameter search (as ranked by final_r_free)

.. note::

   This section may contain less than 10 solutions if for example a solution was found and the molecular replacement search ended early

Summary
^^^^^^^
The summary tab contains three different sections. Below you can find information about each:

.. contents:: Sections
   :depth: 1
   :local:


.. figure:: ../images/ccp4i_lattice_summary.png
   :width: 50%
   :align: center

------------------------------------------------------------------

SIMBAD Summary
==============
This details the best model found by SIMBAD and reports the final_r_fact and final_r_free scores found

Best SIMBAD result Download
===========================
This section contains the refined placed model and mtz for best solution found by the lattice parameter search (as ranked by final_r_free)

Best SIMBAD result Log Files
============================
This section contains the molecular replacement and refinement logs for best solution found by the lattice parameter search (as ranked by final_r_free)
