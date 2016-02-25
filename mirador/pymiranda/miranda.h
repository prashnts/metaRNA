/**
 * Adapted from miRanda.
 *
 * Refactored by: Prashant Sinha (prashant@ducic.ac.in) on 24 Feb 2016
 *
 * Original Authors: Anton Enright, Bino John, Chris Sander and Debora Marks
 * Copyright (C) (2003) Memorial Sloan-Kettering Cancer Center, New York
 * Distributed under the GNU Public License (GPL)
 */

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
void find_targets(char*, char*, ExpString*);

/*   in swat.c */
void traceback(int**, int***, char*, char*, int, int, hit_struct*, double);
int testfor_overlap(int*, int*, int*, int, int);
void build_matrix(int**, int***, int**, int**, int**, int**, char*,
		char*, int, int, score_struct*, int*);
void get_nt_nt_seq_scores(int**, char*, char*, int, int);
void initialize_bases();
void initialize_scores();

/*   in utils.c */
void revstring(char s[]);
void clear_hit(hit_struct*, int, int);
int cmpscores(const void*, const void*);
void string_toupper(char*);

/*   in output.c */
void printhit(int, hit_struct*, double, ExpString*);

/*   in thermo.c */
double get_energy(hit_struct*);
