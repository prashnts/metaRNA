/* miranda.h
 * -------------------------------------------------------------------
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

/* Adjustable Algorithm Parameters*/
double scale;
int strict;
int debug;
double gap_open;
double gap_extend;
double score_threshold;
double energy_threshold;
int length_5p_for_weighting;
int length_3p_for_weighting;
int key_value_pairs;
int no_energy;
int verbosity;
int outfile;
int truncated;
int restricted;

/* Expandable string type (allocates memory as needed) */
struct ExpStringT;
typedef struct ExpStringT ExpString;
void create_ExpString(ExpString **es);
void destroy_ExpString(ExpString **es);
void clear_ExpString(ExpString *es);
int length_ExpString(ExpString *es);
void append_char_ExpString(ExpString *es, char c);
void append_string_ExpString(ExpString *es, char *c);
void append_int_ExpString(ExpString *es, int i);
char *access_ExpString(ExpString *es);

/* Structure to store individual hit information*/
typedef struct hit_struct
{
	double score;
	int query_start;
	int query_end;
	int ref_start;
	int ref_end;
	char* alignment[3];
	char* rest[6];
} hit_struct;

/* Score structure for non-optimal path detection*/
typedef struct score_struct
{
	int score;
	int path;
	int query_trace_end;
	int reference_trace_end;
} score_struct;

/* Structure for pair-wise restriction*/
typedef struct pair_struct
{
	char identifier1[200];
	char identifier2[200];
} pair_struct;

/* Functions Declarations */
/*   in scan.c */
int find_targets(FILE*, FILE*, pair_struct*, int, char*);
/*   in swat.c */
void traceback(int**, int***, char*, char*, int, int, hit_struct*, double);
int testfor_overlap(int*, int*, int*, int, int);
void build_matrix(int**, int***, int**, int**, int**, int**, char*, char*, int, int, score_struct*, int*);
void get_nt_nt_seq_scores(int**, char*, char*, int, int);
void initialize_bases();
void initialize_scores();
/*   in seqio.c */
int readinseq(FILE*, char**, char**, char**, int *);
void initialize_seqio_buffers();
void destroy_seqio_buffers();
/*   in utils.c */
void revstring(char s[]);
void clear_hit(hit_struct*, int, int);
int cmpscores(const void*, const void*);
void string_toupper(char*);
/*   in output.c */
int parse_command_line(int, char**, char*, char*, char*, char*);
void printhit(char*, int, char*, hit_struct*, double, int, FILE*);
void print_parameters(char*, char*, FILE*);
/*   in pairs.c */
int load_pairs(FILE*, pair_struct**);
int find_pair(char*, char*, int, pair_struct*);
/*   in thermo.c */
double get_energy(hit_struct*);
