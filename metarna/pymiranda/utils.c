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

