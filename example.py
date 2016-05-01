#!/usr/bin/env python
# -*- coding: utf-8 -*-
# metaRNA

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
