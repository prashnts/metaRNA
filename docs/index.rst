.. metaRNA documentation master file, created by
   sphinx-quickstart on Fri Apr 22 00:21:17 2016.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

metaRNA, find target sites for the miRNAs
=========================================

metaRNA finds potential target sites for the microRNAs in genomic
sequences.

* Written in Python
* Built on miRanda_.

It is built on miRanda, an algorithm for detection and
ranking of the targets of microRNA.

.. _miRanda: http://www.microrna.org/microrna/getDownloads.do

Quickstart
----------

.. code:: python

    from metarna.target_scan import scan, free_energy

    gene_sequence = (
        "ACAAGATGCCATTGTCCCCCGGCCTCCTGCTGCTGCTGCTCTCCGGGGCCACGGCCACCGCTGCCCTGCC"
        "CCTGGAGGGTGGCCCCACCGGCCGAGACAGCGAGCATATGCAGGAAGCGGCAGGAATAAGGAAAAGCAGC"
        "CTCCTGACTTTCCTCGCTTGGTGGTTTGAGTGGACCTCCCAGGCCAGTGCCGGGCCCCTCATAGGAGAGG"
    )

    mirna_sequence = "UGGCGAUUUUGGAACUCAAUGGCA"

    # Get free Energy value:
    delta_g = free_energy(gene_sequence, mirna_sequence)

    # Get full targets information:
    targets = scan(gene_sequence, mirna_sequence)

    # Specifying Calculation Parameters
    targets = scan(gene_sequence, mirna_sequence, scale=5.0)


Contents

.. toctree::
   :maxdepth: 2

   installation
   parameters
   getting-started


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

