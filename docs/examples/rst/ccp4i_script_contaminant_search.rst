.. _ccp4i_script_contaminant_search:

Searching for contaminants with SIMBAD
--------------------------------------

.. note::

   Data used throughout this example can be found in ``<ROOT>/examples/contaminant_example``. If SIMBAD is part of your CCP4 installation,
   then the example files can be downloaded as part of the `GitHub repository <https://github.com/rigdenlab/SIMBAD>`_.


0. Command line options
^^^^^^^^^^^^^^^^^^^^^^^
Check out this page explaining the :ref:`simbad-contaminant <simbad_contaminant_options>` script command line options.

1. Running the SIMBAD contaminant script
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
In this example, the ``simbad-contaminant`` script simply takes the crystallographic data file in MTZ format, and runs the contaminant search on your local machine.

SIMBAD can be found under the Molecular Replacement menu in the CCP4i GUI:

.. figure:: ../images/ccp4i_simbad.png
   :width: 50%
   :align: center

Opening SIMBAD will bring you to the following menu:

.. figure:: ../images/ccp4i_simbad_run.png
   :width: 50%
   :align: center

To run the lattice search, all that is needed is an MTZ file:

.. figure:: ../images/ccp4i_simbad_contaminant.png
   :width: 50%
   :align: center

SIMBAD Output
-------------
On starting SIMBAD a separate window will appear summarising the progress of the SIMBAD contaminant search and any results found.
The window will contain up to three tabs, the contents of which are explained below:

.. contents:: Output Tabs
   :depth: 1
   :local:

Log File
^^^^^^^^
This displays the text output by SIMBAD as it is running. Any problems or errors will be displayed here.

.. figure:: ../images/ccp4i_contaminant_log.png
   :width: 50%
   :align: center

------------------------------------------------------------------


Contaminant Search Results
^^^^^^^^^^^^^^^^^^^^^^^^^^
The Contaminant Search Results tab contains 5 different sections. Below you can find information about each:

.. contents:: Sections
   :depth: 1
   :local:

------------------------------------------------------------------

Contaminant database AMORE Rotation Search Results
==================================================

.. figure:: ../images/ccp4i_contaminant_results.png
   :width: 50%
   :align: center

This shows the results from the AMORE Rotation Search carried out on the contaminant database. The columns of the table are:

* **PDB_code:** The 4 letter code representing the protein in the protein data bank
* **ALPHA:** Lattice parameter alpha
* **BETA:** Lattice parameter beta
* **GAMMA:** Lattice parameter gamma
* **CC_F:** he correlation coefficient between the observed amplitudes for the crystal and the calculated amplitudes for the model
* **RF_F:** The classic R factor between the observed amplitudes for the crystal and the calculated amplitudes for the model
* **CC_I:** The correlation coefficient between the observed intensities for the crystal and the sum of calculated intensities for all symmetry equivalents of the model
* **CC_P:** The Patterson correlation coefficient between the crystal and the model Pattersons evaluated within the defined sphere centered on the Patterson origin
* **Icp:**
* **CC_F_Z_score:** Z-score of CC_F peaks
* **CC_P_Z_score:** Z-score of CC_P peaks
* **Number_of_rotation_searches_producing_peak:** Number of rotations searches which produce each peak [out of 5]

The structures are scored by CC_F_Z_score score where a higher score is better.

Molecular Replacement Search Results
====================================
Molecular replacement is performed on the top 20 structures identified by the Contaminant Search. This section displays the results of that molecular replacement.

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

Top 10 Contaminant Search Downloads
===================================
This section contains the refined placed model and mtz for the top 10 solutions in the contaminant database search (as ranked by final_r_free)

.. note::

   This section may contain less than 10 solutions if for example a solution was found and the molecular replacement search ended early

Top 10 Contaminant Search Log Files
===================================
This section contains the molecular replacement and refinement logs for the top 10 solutions in the contaminant database search (as ranked by final_r_free)

.. note::

   This section may contain less than 10 solutions if for example a solution was found and the molecular replacement search ended early

Summary
^^^^^^^
The summary tab contains three different sections. Below you can find information about each:

.. contents:: Sections
   :depth: 1
   :local:


.. figure:: ../images/ccp4i_contaminant_summary.png
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
