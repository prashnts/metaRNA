/* ExpString.c */
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

/* expanding_string - a string module which auto-allocates memory */
const int ExpString_START_MEM = 64;
const int ExpString_GROW_NUMERATOR = 3;
const int ExpString_GROW_DENOMINATOR = 2;

struct ExpStringT {
	int char_capacity;
	char *top;	/* points to string (and memory) beginning */
	char *bottom;	/* points to string end (null terminator) */
};

void create_ExpString(ExpString **es) {
	if (*es != 0) {
		fprintf(stderr,"Error: ExpString created on non-null pointer (overwrite/memory leak)\n");
		exit(1);
	}
	*es = malloc(sizeof(ExpString));
	if (*es == 0) {
		fprintf(stderr,"Error: memory allocation failed\n");
		exit(1);
	}
	(*es)->char_capacity = ExpString_START_MEM;
	(*es)->top = malloc((*es)->char_capacity * sizeof(char));
	if ((*es)->top == 0) {
		fprintf(stderr,"Error: memory allocation failed\n");
		exit(1);
	}
	(*es)->bottom = (*es)->top;
	*((*es)->top) = '\0';
}

void destroy_ExpString(ExpString **es) {
	if (*es == 0) return;
	free((*es)->top);
	free(*es);
	*es = 0;
}

void clear_ExpString(ExpString *es) {
	es->bottom = es->top;
	*(es->top) = '\0';
}

int length_ExpString(ExpString *es) {
	return (es->bottom - es->top) / sizeof(char);
}

void append_char_ExpString(ExpString *es, char c) {
	int start_capacity = es->char_capacity;
	if ((es->bottom - es->top + 1) / sizeof(char) + sizeof(char) > es->char_capacity) {
		es->char_capacity = (es->char_capacity * ExpString_GROW_NUMERATOR) / ExpString_GROW_DENOMINATOR;
		es->top = realloc(es->top,es->char_capacity * sizeof(char));
		if (es->top == 0) {
			fprintf(stderr,"Error: memory allocation failed\n");
			exit(1);
		}
		es->bottom = es->top + (start_capacity - 1) * sizeof(char);
	}
	*(es->bottom) = c;
	es->bottom += sizeof(char);
	*(es->bottom) = '\0';
}

void append_string_ExpString(ExpString *es, char *c) {
	while (*c != '\0') append_char_ExpString(es,*c++);
}

void append_int_ExpString(ExpString *es, int i) {
	static char int_to_char[10] = {'0','1','2','3','4','5','6','7','8','9'};
	int q = i / 10;
	int r = i % 10;
	if (i < 0) {
		append_char_ExpString(es,'-');
		r = i - q * 10;
		r = -r;
		q = -q;
	}
	if (q > 0) {
		append_int_ExpString(es, q);
	}
	append_char_ExpString(es,int_to_char[r]);
}

char *access_ExpString(ExpString *es) {
	return es->top;
}
