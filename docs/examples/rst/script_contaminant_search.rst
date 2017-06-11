.. _script_contaminant_search:

Searching for contaminants with SIMBAD
--------------------------------------

.. note::
   
   Data used throughout this example can be found in ``<ROOT>/examples/contaminant_example``. If SIMBAD is part of your CCP4 installation,
   then the example files can be downloaded as part of the `GitHub repository <https://github.com/rigdenlab/SIMBAD>`_.


0. Command line options
^^^^^^^^^^^^^^^^^^^^^^^
Check out this page explaining the :ref:`simbad-contaminant <simbad_contaminant_options>` script command line options.

1. Running the script
^^^^^^^^^^^^^^^^^^^^^
In this example, the ``simbad-contaminant`` script simply takes the crystallographic data file in MTZ format, and runs the contaminant search on your local machine.

.. literalinclude:: /../examples/contaminant_example/run.sh
   :language: bash
   :lines: 11-12

SIMBAD Output
-------------
On starting SIMBAD a separate window will appear summarising the progress of the SIMBAD contaminant search and any results found.
The window will contain up to three tabs, the contents of which are explained below:

Log File
^^^^^^^^
This displays the text output by SIMBAD as it is running. Any problems or errors will be displayed here.

.. figure:: ../images/contaminant_log.png
   :align: center

------------------------------------------------------------------


Contaminant Search Results
^^^^^^^^^^^^^^^^^^^^^^^^^^
The Contaminant Search Results tab contains 5 different sections. Below you can find information about each:

Contaminant database AMORE Rotation Search Results
==================================================

.. figure:: ../images/contaminant_amore_table.png
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
Molecular replacement is performed on the top 20 structures identified by the contaminant database AMORE Rotation search. This section displays the results of that molecular replacement.

.. figure:: ../images/contaminant_mr_table.png
   :align: center

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

Additionally if there is anomalous signal in your dataset SIMBAD will try to validate the quality of the molecular replacement solution using by plotting the peaks from a phased amonalous fourier map. If run the following columns are added to the table:

* **peaks_over_6_rms:** Anomalous peaks over 6 RMS
* **peaks_over_6_rms_within_4a_of_model:** Anomalous peaks over 6 RMS within 4 Angstroms of the Molecular Replacement solution
* **peaks_over_9_rms:** Anomalous peaks over 9 RMS
* **peaks_over_9_rms_within_4a_of_model:** Anomalous peaks over 9 RMS within 4 Angstroms of the Molecular Replacement solution

Molecular Replacement Search Graphs
===================================
Graphs showing the relationship between the final R-Free and various MR parameters are also presented alongside the Molecular Replacement Search results. These are:

* **R-Fact/R-Free Vs. Rank (by R-free):**

.. figure:: ../images/contaminant_mr_graph.png
   :align: center

.. note::

   Here only 1 results is shown because the contaminant search was run with early terminate set to true. These graphs provide more information for borderline cases. 

If using MOLREP:
* **MOLREP score Vs. Rank (by R-free):**
* **MOLREP TF/sig Vs. Rank (by R-free):**

if using Phaser:
* **PHASER TFZ Vs. Rank (by R-free):**
* **PHASER LLG Vs. Rank (by R-free):**
* **PHASER RFZ Vs. Rank (by R-free):**

Top 10 Contaminant Search Downloads
===================================
This section contains the refined placed model and mtz for the top 10 solutions in the contaminant database search (as ranked by final_r_free)

.. figure:: ../images/contaminant_mr_downloads.png
   :align: center

.. note::

   This section may contain less than 10 solutions if for example a solution was found and the molecular replacement search ended early

Top 10 Contaminant Search Log Files
===================================
This section contains the molecular replacement and refinement logs for the top 10 solutions in the contaminant database search (as ranked by final_r_free)

.. figure:: ../images/contaminant_log_downloads.png
   :align: center

.. note::

   This section may contain less than 10 solutions if for example a solution was found and the molecular replacement search ended early

Summary
^^^^^^^
The summary tab contains three different sections. Below you can find information about each:

.. figure:: ../images/contaminant_summary.png
   :align: center

------------------------------------------------------------------

SIMBAD Summary
==============
This details the best model found by SIMBAD and reports the final_r_fact and final_r_free scores found

Best SIMBAD result Download
===========================
This section contains the refined placed model and mtz for best solution found by the contaminant search (as ranked by final_r_free)

Best SIMBAD result Log Files
============================
This section contains the molecular replacement and refinement logs for best solution found by the contaminant search (as ranked by final_r_free)