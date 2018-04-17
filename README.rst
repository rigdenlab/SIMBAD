**********************************************************************
Sequence Independent Molecular Replacement Based on Available Database
**********************************************************************

.. image:: https://readthedocs.org/projects/simbad/badge/?version=latest
   :target: http://simbad.readthedocs.io/en/latest/?badge=latest
   :alt: Documentation Status

.. image:: https://landscape.io/github/rigdenlab/SIMBAD/master/landscape.svg?style=flat
   :target: https://landscape.io/github/rigdenlab/SIMBAD/master
   :alt: Code Health

.. image:: https://img.shields.io/badge/solution%20count-4-blue.svg?style=flat
   :alt: Solution count

About
+++++

SIMBAD is a sequence independant molecular replacement pipeline developed by the by the group of `Daniel Rigden <https://www.liverpool.ac.uk/integrative-biology/staff/daniel-rigden/>`_ at the at the University of Liverpool.
SIMBAD provides an alternate strategy to identify Molecular Replacement search models in a sequence-independent manner.
This makes it suited to:

* Solve cases of contaminant crystallisation, and other mishaps such as mistaken identity (swapped crystallisation trays),
* Solving unsequenced targets
* Provide a brute-force approach where sequence-dependent search model identification could be non-trivial e.g. because of conformational diversity among identifiable homologues.

Flowchart
+++++++++

SIMBAD implements a three-step pipeline to efficiently identify a suitable search model in a database of known structures.

1. SIMBAD performs a lattice search against the entire `Protein Data Bank <https://www.rcsb.org/>`_, rapidly determining whether or not a homologue exists in the same crystal form.
2. SIMBAD screens the target data for the presence of a crystallised contaminant, a not uncommon occurrence in macromolecular crystallography. To catch for this eventuality, SIMBAD rapidly screens the data against a database of known contaminant structures. This database is compiled to include entries from `ContaBase <https://strube.cbrc.kaust.edu.sa/contaminer/contabase>`_ and `Dimple <https://github.com/ccp4/dimple>`_.
3. The final step in SIMBAD can be invoked to perform a brute-force search of a non-redundant PDB database provided by the `MoRDa MR software <http://www.biomexsolutions.co.uk/morda/>`_.

The following flowchart provides a quick overview of the SIMBAD approach:

.. figure:: _static/flowchart.png
   :width: 50%
   :align: left

Found a Bug?
++++++++++++
Please report all bugs to `CCP4 Help Desk <ccp4@stfc.ac.uk>`_.


.. _GitHub Issue Tracker: https://github.com/rigdenlab/simbad/issues
