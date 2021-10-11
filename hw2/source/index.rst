.. BINF690 documentation master file, created by
   sphinx-quickstart on Sun Sep 20 06:58:03 2020.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.


BINF730 Biological Sequence Analysis
====================================

| Alexander Smith
| School of Systems Biology
| George Mason University
| Fall 2021


Phylogenetic Trees
==================

See the "Source Code" for implementation details.

All source code and documentation can also be found on Github:
https://github.com/ajsmith/sequence-analysis


Distance Functions
==================

Jukes-Cantor or HKY are two methods which can be used to calculate
pairwise distances between two sequences.

This project implements Jukes-Cantor.


Clustering Functions
====================



Source Code
===========

..  literalinclude:: ../../src/coolseq/phylo.py
    :caption: src/coolseq/phylo.py


References
==========

 1. Course materials & textbooks

 2. https://en.wikipedia.org/wiki/UPGMA

 3. https://en.wikipedia.org/wiki/WPGMA

 4. https://en.wikipedia.org/wiki/Neighbor_joining

 5. https://en.wikipedia.org/wiki/Models_of_DNA_evolution#JC69_model_(Jukes_and_Cantor_1969)

 6. https://www.megasoftware.net/web_help_7/hc_jukes_cantor_distance.htm
