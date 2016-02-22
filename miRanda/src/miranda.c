/* miranda.c */
/* -------------------------------------------------------------------
 * miRanda- An miRNA target scanner, aims to predict mRNA targets for microRNAs,
 * using dynamic-programming alignment and thermodynamics
 *
 * Copyright (C) (2003) Memorial Sloan-Kettering Cancer Center, New York
 *
 * Distributed under the GNU Public License (GPL)
 * See the files 'COPYING' and 'LICENSE' for details
 *
 * Authors: Anton Enright, Bino John, Chris Sander and Debora Marks
 * Email: mirnatargets (at) cbio.mskcc.org - reaches all authors
 *
 * Written By: Anton Enright
 *
 * Please send bug reports to: miranda (at) cbio.mskcc.org
 *
 * If you use miRanda in your research please cite:
 * Enright AJ, John B, Gaul U, Tuschl T, Sander C and Marks DS;
 * (2003) Genome Biology; 5(1):R1.
 *
 * This software will be further developed under the open source model,
 * coordinated by Anton Enright and Chris Sander:
 * miranda (at) cbio.mskcc.org (reaches both).
 *
 * Copyright (C) (2003) Memorial Sloan-Kettering Cancer Center
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either
 * version 2 of the License, or (at your option) any later
 * version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 * -------------------------------------------------------------------
 */

#include <config.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "miranda.h"

void initialize_globals() {
  initialize_bases(); /* Prepare the generic base lookup array*/
  initialize_scores();
  initialize_seqio_buffers();
}

void destroy_globals() {
  destroy_seqio_buffers();
}

int main (int argc, char* argv[]) {
  FILE* fpout = stdout;
  int total_pairs = 0;
  pair_struct* pairs = 0;
  /* Set Default Parameter Values*/
  length_5p_for_weighting = 8;  /* The 5' sequence length to be weighed  except for the last residue*/
  scale = 4.0;      /* The 5' miRNA scaling parameter*/
  strict = 0;     /* Strict seed model on/off*/
  debug = 0;      /* Debugging mode on/off*/
  key_value_pairs = 0;
  gap_open = -9.0;    /* Gap-open Penalty*/
  gap_extend = -4.0;    /* Gap-extend Penalty*/
  score_threshold = 100.0;  /* SW Score Threshold for reporting hits*/
  energy_threshold = 1.0;   /* Energy Threshold (DG) for reporting hits*/
  verbosity = 0;      /* Verbose mode on/off*/
  outfile = 0;      /* Dump to file on/off*/
  truncated = 0;      /* Truncate sequences on/off*/
  no_energy = 0;      /* Turn off Vienna Energy Calcs - FASTER*/
  restricted = 0;     /* Perform restricted search space*/
  initialize_globals();
  if (restricted && verbosity) {
    printf("Performing Restricted Scan on:%d pairs\n", total_pairs);
  }
  char* gene_seq = "GCTACAGTTTTTATTTAGCATGGGGATTGCAGAGTGACCAGCACACTGGACTCCGAGGTGGTTCAGACAAGACAGAGGGGAGCAGTGGCCATCATCCTCCCGCCAGGAGCTTCTTCGTTCCTGCGCATATAGACTGTACATTATGAAGAATACCCAGGAAGACTTTGTGACTGTCACTTGCTGCTTTTTCTGCGCTTCAGTAACAAGTGTTGGCAAACGAGACTTTCTCCTGGCCCCTGCCTGCTGGAGATCAGCATGCCTGTCCTTTCAGTCTGATCCATCCATCTCTCTCTTGCCTGAGGGGAAAGAGAGATGGGCCAGGCAGAGAACAGAACTGGAGGCAGTCCATCTA";
  char* mirna_seq = "GCGGUUAUAAAUGCACGACGAU";
  int gene_len = 352, mirna_len = 22;
  find_targets(gene_seq, mirna_seq, gene_len, mirna_len, fpout);
  destroy_globals();
  if (outfile) fclose(fpout);
  return 0;
}
