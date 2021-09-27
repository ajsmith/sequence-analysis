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


Pairwise Sequence Alignment
===========================

All source code and documentation can be found here:
https://github.com/ajsmith/sequence-analysis


Examples By Hand
================

Below are small examples by hand using Needleman-Wunsch,
Waterman-Smith-Beyer, and Smith-Waterman.

.. image:: images/binf730-hw1-nw.jpg

.. image:: images/binf730-hw1-wsb.jpg

.. image:: images/binf730-hw1-sw.jpg


Programmatic Examples
=====================

Implementations of Needleman-Wunsch and Waterman-Smith-Beyer
algorithms are provided. Below are code examples and program output
which demonstrates their behavior.


Needleman-Wunsch
----------------

The `nw_align()` function aligns a pair of sequences using
Needleman-Wunsch by similarity.

    >>> from coolseq.align.pairwise import nw_align, print_alignment
    >>> nw_options = {
    ...     'match': 1,
    ...     'mismatch': -1,
    ...     'gap_extend': -1,
    ... }
    >>> result = nw_align('at', 'aagt', nw_options)
    >>> print_alignment(result)
    a--t
    |  |
    aagt
    >>> result = nw_align('gattaca', 'gcatgcu', nw_options)
    >>> print_alignment(result)
    g-attaca
    | ||  |
    gcatg-cu
    >>> result = nw_align('atgc', 'attgagc', nw_options)
    >>> print_alignment(result)
    at-g--c
    || |  |
    attgagc

Now we adjust the scoring options and re-run the alignments.

    >>> nw_options = {
    ...     'match': 1,
    ...     'mismatch': -10,
    ...     'gap_extend': -5,
    ... }
    >>> result = nw_align('at', 'aagt', nw_options)
    >>> print_alignment(result)
    a--t
    |  |
    aagt
    >>> result = nw_align('gattaca', 'gcatgcu', nw_options)
    >>> print_alignment(result)
    g-atta-ca-
    | ||   |
    gcat--gc-u
    >>> result = nw_align('atgc', 'attgagc', nw_options)
    >>> print_alignment(result)
    at-g--c
    || |  |
    attgagc

Change the options again.

    >>> nw_options = {
    ...     'match': 10,
    ...     'mismatch': -2,
    ...     'gap_extend': -10,
    ... }
    >>> result = nw_align('at', 'aagt', nw_options)
    >>> print_alignment(result)
    a--t
    |  |
    aagt
    >>> result = nw_align('gattaca', 'gcatgcu', nw_options)
    >>> print_alignment(result)
    gattaca
    |  | |
    gcatgcu
    >>> result = nw_align('atgc', 'attgagc', nw_options)
    >>> print_alignment(result)
    at-g--c
    || |  |
    attgagc


Waterman-Smith-Beyer
--------------------

The `wsb_align()` function aligns a pair of sequences using
Waterman-Smith-Beyer by distance.

    >>> from coolseq.align.pairwise import wsb_align
    >>> wsb_options = {
    ...     'match': 0,
    ...     'mismatch': 1,
    ...     'gap_start': 1,
    ...     'gap_extend': 1,
    ... }
    >>> result = wsb_align('at', 'aagt', wsb_options)
    >>> print_alignment(result)
    a--t
    |  |
    aagt
    >>> result = wsb_align('gattaca', 'gcatgcu', wsb_options)
    >>> print_alignment(result)
    gattaca
    |  | |
    gcatgcu
    >>> result = wsb_align('atgc', 'attgagc', wsb_options)
    >>> print_alignment(result)
    at---gc
    ||   ||
    attgagc

With different scoring options.

    >>> wsb_options = {
    ...     'match': 0,
    ...     'mismatch': 5,
    ...     'gap_start': 10,
    ...     'gap_extend': 3,
    ... }
    >>> result = wsb_align('at', 'aagt', wsb_options)
    >>> print_alignment(result)
    a--t
    |  |
    aagt
    >>> result = wsb_align('gattaca', 'gcatgcu', wsb_options)
    >>> print_alignment(result)
    gattaca
    |  | |
    gcatgcu
    >>> result = wsb_align('atgc', 'attgagc', wsb_options)
    >>> print_alignment(result)
    at---gc
    ||   ||
    attgagc

Change the options again.

    >>> wsb_options = {
    ...     'match': 0,
    ...     'mismatch': 10,
    ...     'gap_start': 2,
    ...     'gap_extend': 5,
    ... }
    >>> result = wsb_align('at', 'aagt', wsb_options)
    >>> print_alignment(result)
    a--t
    |  |
    aagt
    >>> result = wsb_align('gattaca', 'gcatgcu', wsb_options)
    >>> print_alignment(result)
    g-attaca
    | ||  |
    gcatg-cu
    >>> result = wsb_align('atgc', 'attgagc', wsb_options)
    >>> print_alignment(result)
    at---gc
    ||   ||
    attgagc


Demonstration With NCBI Sequences
=================================

Ribosomal RNA from saccharomyces cerevisiae (yeast) and drosophila
melanogaster (fruit fly) were chosen for samples.

..  literalinclude:: ../../src/coolseq/align/samples.fasta
    :caption: Sample Sequences

EMBOSS Needle on Fedora Linux was used to produce an alignment to
compare against.

Needle command::

    $ needle -asequence scere.fasta -bsequence dmela.fasta -outfile result.needle
    Needleman-Wunsch global alignment of two sequences
    Gap opening penalty [10.0]: 10.0
    Gap extension penalty [0.5]: 0.5

Which produced the following alignment:

..  literalinclude:: result.needle
    :caption: EMBOSS Needle Result


The following code examples demonstrates the alignment algorithms
provided by this project.

1. Load the samples

    >>> from Bio import SeqIO
    >>> from pkg_resources import resource_filename
    >>> file_path = resource_filename('coolseq', 'align/samples.fasta')
    >>> sample1, sample2 = list(SeqIO.parse(file_path, 'fasta'))
    >>> print(sample1.description)
    D25212.1 Saccharomyces cerevisiae gene for 26S rRNA, partial sequence, strain: IFO 2376
    >>> print(sample2.description)
    AY319386.1 Drosophila melanogaster 28S ribosomal RNA, partial sequence

2. Align samples using Needleman-Wunsch

    >>> result = nw_align(sample1.seq, sample2.seq)
    >>> print_alignment(result)
    ACTTG-GATATGGATTC--TTCA-CGGTAACGTAACT--GA--ATGTGG-AGACGTCG-GC-GCGAGCCCTG
     |  | || | | |  |  |||| |  ||| || |||  |   || ||| | | || | |  || ||   ||
    -CGAGCGAAAAGAAAACAGTTCAGCACTAA-GTCACTTTGTCTATATGGCAAATGT-GAGATGC-AG---TG
    <BLANKLINE>
    GGA-GGAGTTATCTTTTCTTCTTAACAGCTTATCACCCCGGAATTG--GTTTATCCGGAGATGGGGTCTTAT
      | ||||   ||  |  |||| |   | | || |    | ||||   | || |    | | |   || | |
    T-ATGGAGCG-TCAATA-TTCT-A---G-T-ATGA----GAAATTAACGATT-T----A-A-G---TCCT-T
    <BLANKLINE>
    GGCTGGAA-GAGGCCAGC-ACCTTT-GCTGGCT-CCGGTGCGCTTGTGACGGCCCGTGAAA-ATCCAC-AGG
      ||  || |||||||   |||  | |  || | || | || | | | |||     | ||  ||  || ||
    --CTTAAATGAGGCCATTTACCCATAGAGGG-TGCCAG-GCCCGTATAACGT----T-AATGATT-ACTAG-
    <BLANKLINE>
    AAGGAATAGTTTTCAT-GC-TAG-GT--C--G-TAC-T-
    | | | | |||| ||  |  | | ||  |  | ||| |
    ATG-A-T-GTTTCCAAAGAGTCGTGTTGCTTGATACGTG

3. Align samples using Waterman-Smith-Beyer

    >>> result = wsb_align(sample1.seq, sample2.seq)
    >>> print_alignment(result)
    ACTTGGATATGGATTC--TTCA-CGGTAACGTAACTGAATGTGGA-GACGTCGGCGCGA-GCCCTGGGAGGA
         || | | |  |  |||| |  ||| || |||   | |  | | |    | | || ||  ||   |||
    CGAGCGAAAAGAAAACAGTTCAGCACTAA-GTCACTTTGTCTATATGGCAAATGTGAGATGCAGTGTATGGA
    <BLANKLINE>
    GT-TATCTTTTCTTCTTAACAGCTTATCACCCCGGAATTGGTTTATCCGGAGATGGGGTCTTATGGC--TGG
    |  |   | ||||  | |  ||   || | |   || ||   |  | |  | ||| || | | |  |  | |
    GCGTCAATATTCTAGT-ATGAGAA-ATTAAC---GATTTAAGTCCTTCTTAAATGAGGCCATTTACCCATAG
    <BLANKLINE>
    AAGAGGCCAGCACCTTTGCTGGCTCCGGTGCGCTTGTGACGGCCCGTGAAAATCCACAGGAAGGAATAGTTT
    | |  |||||  || |     | |   |    || |  | |     ||    |||| ||  ||   | || |
    AGGGTGCCAGGCCCGTATAACGTTAATGATTACTAG--ATGA----TGTT--TCCAAAG--AG---TCGTGT
    <BLANKLINE>
    TCATGCTAGGT-CGTACT
    |   ||| | | |||
    T---GCTTGATACGTG--


Source Code
===========

..  literalinclude:: ../../src/coolseq/align/pairwise.py
    :caption: src/coolseq/align/pairwise.py


References
==========

 1. Course materials & textbooks

 2. https://en.wikipedia.org/wiki/Needleman%E2%80%93Wunsch_algorithm

 3. https://en.wikipedia.org/wiki/Smith%E2%80%93Waterman_algorithm
