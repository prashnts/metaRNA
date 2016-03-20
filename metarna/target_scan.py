#!/usr/bin/env python
# -*- coding: utf-8 -*-
#.--. .-. ... .... -. - ... .-.-.- .. -.

import json
from metarna.pymiranda import find_targets

def scan(target_sequence, mirna_sequence, **kwa):
  """
  Find targets of miRNA on the target based on the mirna and target sequences.

  Args:
    target_sequence (str): The target genomic sequence.
    mirna_sequence (str): The microRNA sequence.

  Optional Args: These parameters are passed on to the pymiranda algorithm.
    scale (float): The 5' miRNA scaling parameter. (default: 4.0)
    strict (int): Perform a Strict Seed search when set to 1. (default: 0)
    gap_open (float): Gap-open Penalty (default: -9.0)
    gap_extend (float): Gap-extend Penalty (default: -4.0)
    score_threshold (float): Score Threshold for reporting hits (default: 50.0)
    energy_threshold (float): Energy Threshold for reporting hits (default: 1.0)
    length_5p_for_weighting (int): The 5' sequence length to be weighed except
                                   for the last residue (default: 8).
    temperature (int): Used while calculating Free Energy (default: 30)
    alignment_len_threshold (int): Minimum alignment. (default: 8)

  Return (dict): A dictionary containing the stats of each possible target
                 site and the pairing obtained.
  """
  mirna_sequence_rev = mirna_sequence[::-1]
  enc_result = find_targets(target_sequence, mirna_sequence_rev, **kwa)
  return json.loads(enc_result)

def free_energy(target_sequence, mirna_sequence, **kwa):
  """
  Shortcut to find the Maximum Free energy in the given pair.

  See metarna.pymiranda.target_scan.scan for details on the arguments.

  Return (float): The free energy value.
  """
  try:
    doc = scan(target_sequence, mirna_sequence, **kwa)
    dat = max(doc.get('hits', []), key= lambda x: x['score'])
    return dat['energy']
  except KeyError:
    raise ValueError("No targets found.")
