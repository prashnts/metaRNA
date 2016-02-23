/* output.c */
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
#include <math.h>
#include <string.h>
#include "miranda.h"

extern double scale;
extern int strict;
extern int debug;
extern double gap_open;
extern double gap_extend;
extern double score_threshold;
extern double energy_threshold;
extern int length_5p_for_weighting;
extern int length_3p_for_weighting;
extern int key_value_pairs;
extern int no_energy;
extern int verbosity;
extern int truncated;
extern int restricted;

void printhit(int query_length, hit_struct* hit, double energy,
		int keyval_mode, ExpString* outjson) {
	double similarity = 0;
	double identity = 0;
	int alignment_length = 0;
	int i = 0;
	alignment_length = strlen(hit->alignment[0]);

	for (i = 0; i < alignment_length; i++) {
		if (hit->alignment[1][i] == '|') {
			similarity++;
			identity++;
		}
		if (hit->alignment[1][i] == ':') {
			similarity++;
		}
	}

	similarity = (similarity / (double)alignment_length) * 100;
	identity = (identity / (double)alignment_length) * 100;

	revstring(hit->alignment[0]);
	revstring(hit->alignment[1]);
	revstring(hit->alignment[2]);

	char temp[512];
	sprintf(temp,
		"{\"score\": %f, \"energy\": %f, \"query_start\": %d, \"query_end\": %d, "
		"\"ref_start\": %d, \"ref_end\": %d, \"aln_length\": %d, "
		"\"identity\": %f, \"similarity\": %f, \"aln_mirna\": \"%s%s%s\", "
		"\"aln_map\": \"%s%s%s\", \"aln_utr\": \"%s%s%s\"}",
			hit->score, energy, (query_length - hit->query_end + 1),
			(query_length - hit->query_start + 1), hit->ref_start + 1,
			hit->ref_end + 1, alignment_length, identity, similarity,
			hit->rest[0], hit->alignment[0], hit->rest[3],
			hit->rest[2], hit->alignment[1], hit->rest[5],
			hit->rest[1], hit->alignment[2], hit->rest[4]
		);
	append_string_ExpString(outjson, temp);
}
