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

Algorithms such as WPGMA, UPGMA, and Neighbor Distance are useful for
estimating the evolutionary relationships between different taxa. The
information these algorithms produce can be used to produce
phylogenetic trees, which provide a graphical representation of those
relationships.

See the "Source Code" for implementation details.

This document is includes Python doctests which can be used to verify
the correctness of the code examples shown below. Additional tests
used are included in the "Unit Tests" section below.

All source code and documentation can also be found on Github:
https://github.com/ajsmith/sequence-analysis


Distance Functions
==================

Jukes-Cantor or HKY are two methods which can be used to calculate
pairwise distances between two sequences.

This project implements Jukes-Cantor.


Sample Data
===========

Example 1 is the distance matrix example for WPGMA and UPGMA.

    >>> from coolseq.phylo import get_examples
    >>> from coolseq.utils import print_matrix
    >>> examples = get_examples()
    >>> matrix, names = examples['example1']
    >>> print_matrix(matrix)
    [0, 6, 10, 10, 10]
    [6, 0, 10, 10, 10]
    [10, 10, 0, 2, 6]
    [10, 10, 2, 0, 6]
    [10, 10, 6, 6, 0]
    >>> print(names)
    ['A', 'B', 'C', 'D', 'E']

Example 2 is the distance matrix used in the Neighbor Joining example.

    >>> matrix, names = examples['example2']
    >>> print_matrix(matrix)
    [0, 5, 4, 7, 6, 8]
    [5, 0, 7, 10, 9, 11]
    [4, 7, 0, 7, 6, 8]
    [7, 10, 7, 0, 5, 9]
    [6, 9, 6, 5, 0, 8]
    [8, 11, 8, 9, 8, 0]
    >>> print(names)
    ['A', 'B', 'C', 'D', 'E', 'F']

Example 3 is the distance matrix for primate data.

    >>> matrix, names = examples['example3']
    >>> print_matrix(matrix)
    [0, 79, 92, 144, 162]
    [79, 0, 95, 154, 169]
    [92, 95, 0, 150, 169]
    [144, 154, 150, 0, 169]
    [162, 169, 169, 169, 0]
    >>> print(names)
    ['Human', 'Chimp', 'Gorilla', 'Orangutan', 'Gibbon']

Example 4 is a collection of rRNA 5S subunit sequences from various
organisms. To process the FASTA file, we use Biopython.

    >>> from Bio import SeqIO
    >>> from pkg_resources import resource_filename
    >>> file_path = resource_filename('coolseq', 'samples/rRNA-5S.fasta')
    >>> rrna_samples = list(SeqIO.parse(file_path, 'fasta'))
    >>> for (i, sample) in enumerate(rrna_samples):
    ...     print(f'{i}: {sample.description}')
    0: NR_131385.1 Caenorhabditis elegans 5S ribosomal RNA (rrn-4.16), rRNA
    1: NR_033366.1 Cavia porcellus 5S ribosomal RNA (LOC100379616), ribosomal RNA
    2: NR_001870.2 Drosophila melanogaster 5S ribosomal RNA (5SrRNA:CR33375), rRNA
    3: NR_023378.1 Homo sapiens RNA, 5S ribosomal 16 (RNA5S16), ribosomal RNA
    4: NR_046144.1 Mus musculus nuclear encoded rRNA 5S 136 (n-R5s136), ribosomal RNA
    5: NR_033176.2 Rattus norvegicus 5S RNA (Rn5s), ribosomal RNA
    >>> matrix, names = examples['example4']
    >>> print_matrix(matrix)
    [0, 0.51369, 0.83499, 0.51369, 0.51369, 0.70835]
    [0.51369, 0, 0.37697, 0, 0.00831, 0.06921]
    [0.83499, 0.37697, 0, 0.37697, 0.38932, 0.44084]
    [0.51369, 0, 0.37697, 0, 0.00831, 0.06921]
    [0.51369, 0.00831, 0.38932, 0.00831, 0, 0.07833]
    [0.70835, 0.06921, 0.44084, 0.06921, 0.07833, 0]
    >>> print(names)
    ['NR_131385.1', 'NR_033366.1', 'NR_001870.2', 'NR_023378.1', 'NR_046144.1', 'NR_033176.2']

To build a distance matrix from the samples, first we generate
pairwise alignments for each pair combination using
Waterman-Smith-Beyer. Then we apply Jukes-Cantor to each alignment to
get a distance, which we enter into the distance matrix. All of this
is packed into the `generate_distance_matrix()` function seen in the
program source code.

Additional details about the aligned sample data are provided in
Appendix I.


Clustering Functions
====================

This project implements WPGMA, UPGMA, and Neighbor Joining. All
implementations return a tuple representation of the tree. WPGMA and
UPGMA also provide linkage data in Scipy format which is used to draw
dendrograms of the phylogenetic trees.

For WPGMA and UPGMA, the tuple representation is such that a branch is
represented as a triple of (HEIGHT, L_TREE, R_TREE), and a leaf is
represented as a pair of (TAXA_INDEX, TAXA_LABEL).

For Neighbor Joining, the tuple representation of a branch is a nested
tuple of ((LEFT_DISTANCE, RIGHT_DISTANCE), LEFT_TREE, RIGHT_TREE). A
leaf is still represented as a pair of (TAXA_INDEX, TAXA_LABEL).


WPGMA
-----

WPGMA Example 1
+++++++++++++++

This tree seems to be consistent with the example from lecture.

    >>> from coolseq.phylo import print_wpgma_ex1
    >>> print_wpgma_ex1()
    (5.0,
      (3.0,
        (0, 'A'),
        (1, 'B')),
      (3.0,
        (4, 'E'),
        (1.0,
          (2, 'C'),
          (3, 'D'))))

..  plot::
    :include-source:

    from coolseq.phylo import plot_wpgma_ex1
    plot_wpgma_ex1()


WPGMA Example 3
+++++++++++++++

The WPGMA results for the primate data seem similar to that of the
minimum evolution tree from the original study, both in terms of tree
structure and distances. There are some differences however.

    >>> from coolseq.phylo import print_wpgma_ex3
    >>> print_wpgma_ex3()
    (84.0625,
      (4, 'Gibbon'),
      (74.75,
        (3, 'Orangutan'),
        (46.75,
          (2, 'Gorilla'),
          (39.5,
            (0, 'Human'),
            (1, 'Chimp')))))

..  plot::
    :include-source:

    from coolseq.phylo import plot_wpgma_ex3
    plot_wpgma_ex3()


WPGMA Example 4
+++++++++++++++

    >>> from coolseq.phylo import print_wpgma_ex4
    >>> print_wpgma_ex4()
    (0.3615025,
      (0, 'NR_131385.1'),
      (0.20599625,
        (2, 'NR_001870.2'),
        (0.036885,
          (5, 'NR_033176.2'),
          (0.004155,
            (4, 'NR_046144.1'),
            (0.0,
              (1, 'NR_033366.1'),
              (3, 'NR_023378.1'))))))

..  plot::
    :include-source:

    from coolseq.phylo import plot_wpgma_ex4
    plot_wpgma_ex4()


UPGMA
-----

UPGMA Example 1
+++++++++++++++

This tree seems to be consistent with the example from lecture.

    >>> from coolseq.phylo import print_upgma_ex1
    >>> print_upgma_ex1()
    (5.0,
      (3.0,
        (0, 'A'),
        (1, 'B')),
      (3.0,
        (4, 'E'),
        (1.0,
          (2, 'C'),
          (3, 'D'))))

..  plot::
    :include-source:

    from coolseq.phylo import plot_upgma_ex1
    plot_upgma_ex1()


UPGMA Example 3
+++++++++++++++

The UPGMA results for the primate data seem similar to that of the
minimum evolution tree from the original study, both in terms of tree
structure and distances. There are some differences however.

    >>> from coolseq.phylo import print_upgma_ex3
    >>> print_upgma_ex3()
    (83.625,
      (4, 'Gibbon'),
      (74.66666666666667,
        (3, 'Orangutan'),
        (46.75,
          (2, 'Gorilla'),
          (39.5,
            (0, 'Human'),
            (1, 'Chimp')))))

..  plot::
    :include-source:

    from coolseq.phylo import plot_upgma_ex3
    plot_upgma_ex3()


UPGMA Example 4
+++++++++++++++

Caenorhabditis is branched from the root without touching other
clusters, suggesting it is evolutionarily distant from the other
samples. The 2nd most distant is drosophilia. The 3rd most distant is
brown rat. The 4th most distant is mouse. Guinea pig and human share
the innermost cluster.

    >>> from coolseq.phylo import print_upgma_ex4
    >>> print_upgma_ex4()
    (0.30844099999999997,
      (0, 'NR_131385.1'),
      (0.19801250000000004,
        (2, 'NR_001870.2'),
        (0.036125,
          (5, 'NR_033176.2'),
          (0.004155,
            (4, 'NR_046144.1'),
            (0.0,
              (1, 'NR_033366.1'),
              (3, 'NR_023378.1'))))))

..  plot::
    :include-source:

    from coolseq.phylo import plot_upgma_ex4
    plot_upgma_ex4()


Neighbor Joining
----------------

This tree seems to be consistent with the example from lecture.

Neighbor Joining Example 2
++++++++++++++++++++++++++

    >>> from coolseq.phylo import print_nj_ex2
    >>> print_nj_ex2()
    ((0.5, 0.5),
      ((3.0, 2.0),
        (3, 'D'),
        (4, 'E')),
        ((5.0, 1.0),
          (5, 'F'),
          ((2.0, 1.0),
            (2, 'C'),
            ((1.0, 4.0),
              (0, 'A'),
              (1, 'B')))))


Neighbor Joining Example 3
++++++++++++++++++++++++++

To the right of the root human and chimp are clustered together, and
to the left is gorilla clustered with a sub-cluster which contains
orangutan and gibbon. This structure and the distances are similar to
that of the minimum evolution tree results.

    >>> from coolseq.phylo import print_nj_ex3
    >>> print_nj_ex3()
    ((3.0625, 3.0625),
      ((36.625, 42.375),
        (0, 'Human'),
        (1, 'Chimp')),
      ((47.875, 27.125),
        (2, 'Gorilla'),
        ((75.83333333333333, 93.16666666666667),
          (3, 'Orangutan'),
          (4, 'Gibbon'))))


Neighbor Joining Example 4
++++++++++++++++++++++++++

Brown rat and guinea pig are clustered close together, and their
cluster is connected to another node with a branch to human and the
other branch leading to a cluster containing mouse. Drosophilia and
caenorhabditis are clustered the furthest away.

    >>> from coolseq.phylo import print_nj_ex4
    >>> print_nj_ex4()
    ((0.0038040624999999984, 0.0038040624999999984),
      ((-0.007608125, 0.07681812499999999),
        (1, 'NR_033366.1'),
        (5, 'NR_033176.2')),
      ((-0.007608124999999997, 0.008675625000000003),
        (3, 'NR_023378.1'),
        ((-0.002631666666666671, 0.036641666666666656),
          (4, 'NR_046144.1'),
          ((0.50066, 0.33433),
            (0, 'NR_131385.1'),
            (2, 'NR_001870.2')))))


References
==========

 1. Course materials & textbooks

 2. https://en.wikipedia.org/wiki/UPGMA

 3. https://en.wikipedia.org/wiki/WPGMA

 4. https://en.wikipedia.org/wiki/Neighbor_joining

 5. https://en.wikipedia.org/wiki/Models_of_DNA_evolution#JC69_model_(Jukes_and_Cantor_1969)

 6. https://www.megasoftware.net/web_help_7/hc_jukes_cantor_distance.htm

 7. https://www.scipy.org/

 8. https://matplotlib.org/


Source Code
===========

..  literalinclude:: ../../src/coolseq/phylo.py
    :caption: src/coolseq/phylo.py


Unit Tests
==========

..  literalinclude:: ../../tests/phylo.rst
    :caption: tests/phylo.rst


Appendix I
==========

Verbose details for sample sequence alignment and distance.

..  literalinclude:: alignment.output
    :caption: alignment.output
