**********************************************************************
Sequence Independent Molecular Replacement Based on Available Database
**********************************************************************

.. image:: https://readthedocs.org/projects/simbad/badge/?version=latest
   :target: http://simbad.readthedocs.io/en/latest/?badge=latest
   :alt: Documentation Status

.. image:: https://travis-ci.com/rigdenlab/SIMBAD.svg?branch=master
   :target: https://travis-ci.com/rigdenlab/SIMBAD
   :alt: CI Status

.. image:: https://img.shields.io/badge/solution%20count-8-blue.svg?style=flat
   :alt: Solution count

.. image:: https://img.shields.io/badge/DOI-10.1107/S2059798318005752-blue.svg
   :target: https://doi.org/10.1107/S2059798318005752
   :alt: SIMBAD paper

.. image:: https://codecov.io/gh/rigdenlab/SIMBAD/branch/master/graph/badge.svg
  :target: https://codecov.io/gh/rigdenlab/SIMBAD


About
+++++

SIMBAD is a sequence independant molecular replacement pipeline developed by the group of `Daniel Rigden <https://www.liverpool.ac.uk/integrative-biology/staff/daniel-rigden/>`_ at the University of Liverpool.
SIMBAD provides an alternate strategy to identify Molecular Replacement search models in a sequence-independent manner.

This makes SIMBAD suited to:

* Solve cases of contaminant crystallisation, and other mishaps such as mistaken identity (swapped crystallisation trays),
* Solving unsequenced targets
* Provide a brute-force approach where sequence-dependent search model identification could be non-trivial e.g. because of conformational diversity among identifiable homologues.

For an overview of SIMBAD, watch `this video <https://www.youtube.com/watch?v=HYGe7541qeQ>`_.

Flowchart
+++++++++

SIMBAD implements a three-step pipeline to efficiently identify a suitable search model in a database of known structures.

1. SIMBAD performs a lattice search against the entire `Protein Data Bank <https://www.rcsb.org/>`_, rapidly determining whether or not a homologue exists in the same crystal form.

2. SIMBAD screens the target data for the presence of a crystallised contaminant, a not uncommon occurrence in macromolecular crystallography. To catch for this eventuality, SIMBAD rapidly screens the data against a database of known contaminant structures. This database is compiled to include entries from `ContaBase <https://strube.cbrc.kaust.edu.sa/contaminer/contabase>`_ and `Dimple <https://github.com/ccp4/dimple>`_.

3. The final step in SIMBAD can be invoked to perform a brute-force search of a non-redundant PDB database provided by the `MoRDa MR software <http://www.biomexsolutions.co.uk/morda/>`_.

The following flowchart provides a quick overview of the SIMBAD approach:

.. image:: https://github.com/rigdenlab/SIMBAD/raw/master/docs/_static/flowchart.png
   :width: 50%
   :align: center

