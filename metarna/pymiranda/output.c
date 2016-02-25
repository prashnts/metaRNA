/**
 * Adapted from miRanda.
 *
 * Refactored by: Prashant Sinha (prashant@ducic.ac.in) on 24 Feb 2016
 *
 * Original Authors: Anton Enright, Bino John, Chris Sander and Debora Marks
 * Copyright (C) (2003) Memorial Sloan-Kettering Cancer Center, New York
 * Distributed under the GNU Public License (GPL)
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "miranda.h"

void printhit(int query_length, hit_struct* hit, double energy,
		ExpString* outjson) {
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
