#!/usr/bin/env python
# -*- coding: utf-8 -*-
# metaRNA
'''Target Scan Report Generation

Each hits contain the following three attributes:
 - aln_mirna
 - aln_map
 - aln_utr
Where aln_utr is the sequence from the cDNA, aln_map are the base pairs, and
aln_mirna is the miRNA sequence.

Expected Render:
  acgguaacucaagguuuUAGCGGu
                   ||:|:|
  -------acaagatgccATTGTCc
'''

