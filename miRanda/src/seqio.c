/* seqio.c */
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

/* invalid nucleotide matrix */
/* '!' indicates an illegal character */

char nt_invalid[256];

int initialize_nt_invalid() {
	int i;
	for (i = 0; i < 256; i++) nt_invalid[i] = '\1';
	nt_invalid['a'] = '\0';
	nt_invalid['A'] = '\0';
	nt_invalid['c'] = '\0';
	nt_invalid['C'] = '\0';
	nt_invalid['g'] = '\0';
	nt_invalid['G'] = '\0';
	nt_invalid['u'] = '\0';
	nt_invalid['U'] = '\0';
	nt_invalid['t'] = '\0';
	nt_invalid['T'] = '\0';
	nt_invalid['n'] = '\0';
	nt_invalid['N'] = '\0';
	nt_invalid['r'] = '\0';
	nt_invalid['R'] = '\0';
	nt_invalid['y'] = '\0';
	nt_invalid['Y'] = '\0';
	nt_invalid['k'] = '\0';
	nt_invalid['K'] = '\0';
	nt_invalid['m'] = '\0';
	nt_invalid['M'] = '\0';
	nt_invalid['s'] = '\0';
	nt_invalid['S'] = '\0';
	nt_invalid['w'] = '\0';
	nt_invalid['W'] = '\0';
	nt_invalid['b'] = '\0';
	nt_invalid['B'] = '\0';
	nt_invalid['d'] = '\0';
	nt_invalid['D'] = '\0';
	nt_invalid['h'] = '\0';
	nt_invalid['H'] = '\0';
	nt_invalid['v'] = '\0';
	nt_invalid['V'] = '\0';
	nt_invalid['x'] = '\0';
	nt_invalid['X'] = '\0';
	return 0;
}

ExpString *sequence_buffer = 0;
ExpString *description_buffer = 0;
ExpString *name_buffer = 0;

void initialize_seqio_buffers() {
	initialize_nt_invalid();
	create_ExpString(&sequence_buffer);
	create_ExpString(&description_buffer);
	create_ExpString(&name_buffer);
}

void destroy_seqio_buffers() {
	destroy_ExpString(&sequence_buffer);
	destroy_ExpString(&description_buffer);
	destroy_ExpString(&name_buffer);
}

/* returns 0 if EOF is reached, otherwise 1 */
int skip_leading_whitespace(FILE* fp) {
	for (;;) {
		int cbuf = fgetc(fp);
		switch (cbuf) {
			case '\n': case '\r': case '\t': case ' ': break;
			case EOF: return 0;
			default: ungetc((char)cbuf, fp); return 1;
		}
	}
	return 0; /* never reached */
}

int skip_leading_space(FILE* fp) {
	for (;;) {
		int cbuf = fgetc(fp);
		switch (cbuf) {
			case '\r': case '\t': case ' ': break;
			case EOF: return 0;
			default: ungetc((char)cbuf, fp); return 1;
		}
	}
	return 0; /* never reached */
}

/* requires file position to be at '>' */
void read_fasta_name_into_buffer(FILE* fp, ExpString *name_buf) {
	int cbuf = fgetc(fp);
	if (cbuf != '>') {
		fprintf(stderr, "Error: illegal format in sequence file (record does not start with '>')\n");
		exit(1);
	}
	clear_ExpString(name_buf);
	for (;;) {
		cbuf = fgetc(fp);
		if (cbuf == EOF) {
			fprintf(stderr, "Error: illegal format in sequence file (EOF while reading name)\n");
			exit(1);
		}
		if (cbuf == '\n' || cbuf == '\r' || cbuf == '\t' || cbuf == ' ') {
			ungetc((char)cbuf, fp);
			break;
		}
		append_char_ExpString(name_buf, (char)cbuf);
	}
}

void read_fasta_description_into_buffer(FILE* fp, ExpString *description_buf) {
	if (skip_leading_space(fp) == 0) {
		fprintf(stderr, "Error: illegal format in sequence file (EOF while reading description)\n");
		exit(1);
	}
	clear_ExpString(description_buf);
	for (;;) {
		int cbuf = fgetc(fp);
		if (cbuf == EOF) {
			fprintf(stderr, "Error: illegal format in sequence file (EOF while reading description)\n");
			exit(1);
		}
		if (cbuf == '\n') return;
		append_char_ExpString(description_buf, (char)cbuf);
	}
}

void read_fasta_sequence_into_buffer(FILE* fp, ExpString *sequence_buf, ExpString *name_buf) {
	clear_ExpString(sequence_buf);
	for (;;) {
		int cbuf = fgetc(fp);
		if (cbuf == EOF) return;
		if (cbuf == '\n' || cbuf == '\r' || cbuf == '\t' || cbuf == ' ') continue;
		if (cbuf == '>') {
			ungetc((char)cbuf,fp);
			return;
		}
		if (nt_invalid[cbuf]) {
			fprintf(stderr, "Error: illegal character '%c' in sequence file in sequence: %s\n", (char)cbuf, access_ExpString(name_buf));
			exit(1);
		}
		append_char_ExpString(sequence_buf, (char)cbuf);
	}
}

int readinseq(FILE* fp, char** sequence, char** description, char** name, int* seqlen) {
	if (skip_leading_whitespace(fp) == 0) return 0; /* end of file reached */
	read_fasta_name_into_buffer(fp, name_buffer);
	read_fasta_description_into_buffer(fp, description_buffer);
	read_fasta_sequence_into_buffer(fp, sequence_buffer, name_buffer);
	if (length_ExpString(name_buffer) == 0) {
		fprintf(stderr, "Error: sequence file contains a record with a length 0 name\n");
		exit(1);
	}
	if (length_ExpString(sequence_buffer) == 0) {
		fprintf(stderr, "Error: sequence file contains a record with a length 0 sequence\n");
		exit(1);
	}
	*name = realloc(*name, (length_ExpString(name_buffer) + 1) * sizeof(char));
	*description = realloc(*description, (length_ExpString(description_buffer) + 1) * sizeof(char));
	*sequence = realloc(*sequence, (length_ExpString(sequence_buffer) + 1) * sizeof(char));
	if (*name == 0 || *description == 0 || *sequence == 0) {
		fprintf(stderr,"Error: memory allocation failed\n");
		exit(1);
	}
	strcpy(*name, access_ExpString(name_buffer));
	strcpy(*description, access_ExpString(description_buffer));
	strcpy(*sequence, access_ExpString(sequence_buffer));
	*seqlen = length_ExpString(sequence_buffer);
	return 1;
}
