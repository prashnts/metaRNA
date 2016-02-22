/* pairs.c */
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
#include <math.h>
#include "miranda.h"

/* Load in a flat-file containing miRNA / UTR pairs to scan*/
int readinpairs(FILE* fp , pair_struct** pairs) {
	char BUFFER[400];
	int n_pairs = 0;
	int N_pairs = 16;
	if (!(*pairs = malloc(N_pairs * sizeof(pair_struct)))) {
		fprintf(stderr, "Cannot obtain %ld pair structs\n", (long)N_pairs);
		exit(1);
	}
	while (fgets(BUFFER, 400, fp)) {
		if (n_pairs >= N_pairs) {
			N_pairs *= 1.44;
			if (!(*pairs = realloc(*pairs, N_pairs * sizeof(pair_struct)))) {
				fprintf(stderr, "Cannot realloc to %ld pair structs\n", (long)N_pairs);
				exit(1);
			}
		}
		sscanf(BUFFER, "%s\t%s", (*pairs)[n_pairs].identifier1, (*pairs)[n_pairs].identifier2);
		n_pairs++;
	}
	return n_pairs;
}


int sort_pairs(const void* s1, const void* s2) {
	pair_struct* a = (pair_struct*)s1;
	pair_struct* b = (pair_struct*)s2;
	if (strcmp(a->identifier1, b->identifier1) == 0) {
		return strcmp(a->identifier2, b->identifier2);
	} else {
		return strcmp(a->identifier1, b->identifier1);
	}
}

int load_pairs(FILE* fp_pairs, pair_struct** pairs) {
	int total_pairs = readinpairs(fp_pairs, pairs);
	qsort(*pairs, total_pairs, sizeof(pair_struct), sort_pairs);
	return total_pairs;
}

int find_pair(char* ident1, char* ident2, int total_pairs, pair_struct* pairs) {
	pair_struct search_key;
	strcpy(search_key.identifier1, ident1);
	strcpy(search_key.identifier2, ident2);
	return bsearch(&search_key, pairs, total_pairs, sizeof(pair_struct), sort_pairs) != 0;
}
