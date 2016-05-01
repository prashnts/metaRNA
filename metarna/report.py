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

HIT_TEMPLATE = '''

  ΔG    {energy}

{query_end:>6}
     {aln_mirna}
     {aln_map}
     {aln_utr}
{ref_end:>6}   {ref_start}

'''

def get_report(targets):
  hits = targets.get('hits', [])
  hits.sort(key=lambda x: x['score'], reverse=True)
  out = ""
  for hit in hits:
    out += HIT_TEMPLATE.format(**hit)
  return out
