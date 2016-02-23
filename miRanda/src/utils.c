/* utils.c */
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
#include <string.h>
#include <ctype.h>
#include "miranda.h"

int cmpscores(const void* p1, const void* p2) {
	score_struct* s1 = (score_struct*)p1;
	score_struct* s2 = (score_struct*)p2;
	/* primary key: score */
	if (s1->score < s2->score) return 1;
	if (s1->score > s2->score) return -1;
	/* secondary key: earlier end in reference */
	if (s1->reference_trace_end > s2->reference_trace_end) return 1;
	if (s1->reference_trace_end < s2->reference_trace_end) return -1;
	return 0;
}

void clear_hit(hit_struct* hit, int query_length, int reference_length) {
	hit->score = 0;
	hit->query_start = 0;
	hit->query_end = 0;
	hit->ref_start = 0;
	hit->ref_end = 0;
	memset(hit->alignment[0], '\0', query_length + reference_length);
	memset(hit->alignment[1], '\0', query_length + reference_length);
	memset(hit->alignment[2], '\0', query_length + reference_length);
	memset(hit->rest[0], '\0', query_length + 10);
	memset(hit->rest[1], '\0', query_length + 10);
	memset(hit->rest[2], '\0', query_length + 10);
	memset(hit->rest[3], '\0', query_length + 10);
	memset(hit->rest[4], '\0', query_length + 10);
	memset(hit->rest[5], '\0', query_length + 10);
}

void revstring(char s[]) {
	int c, i, j;
	for (i = 0, j = strlen(s) - 1; i < j; i++, j--) {
		c = s[i];
		s[i] = s[j];
		s[j] = c;
	}
}

void string_toupper(char *s) {
	while (*s) {
		*s = toupper(*s);
		s++;
	}
}

