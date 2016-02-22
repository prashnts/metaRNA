/* swat.c */
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
#include <ctype.h>
#include <string.h>
#include <limits.h>
#include "miranda.h"

/* Master list of Nucleotide Bases*/
int bases[256];

void initialize_bases() {
	int i;
	for (i = 0; i < 256; i++) bases[i] = INT_MAX;
	bases['c'] = 0;
	bases['C'] = 0;
	bases['g'] = 1;
	bases['G'] = 1;
	bases['a'] = 2;
	bases['A'] = 2;
	bases['u'] = 3;
	bases['U'] = 3;
	bases['t'] = 3;
	bases['T'] = 3;
	bases['n'] = 4;
	bases['N'] = 4;
	bases['r'] = 4;
	bases['R'] = 4;
	bases['y'] = 4;
	bases['Y'] = 4;
	bases['k'] = 4;
	bases['K'] = 4;
	bases['m'] = 4;
	bases['M'] = 4;
	bases['s'] = 4;
	bases['S'] = 4;
	bases['w'] = 4;
	bases['W'] = 4;
	bases['b'] = 4;
	bases['B'] = 4;
	bases['d'] = 4;
	bases['D'] = 4;
	bases['h'] = 4;
	bases['H'] = 4;
	bases['v'] = 4;
	bases['V'] = 4;
	bases['x'] = 4;
	bases['X'] = 4;
}

const int ALIGNMENT_MATCH_SCORE = 5;
const int ALIGNMENT_MISMATCH_SCORE = -3;
const int ALIGNMENT_WOBBLE_SCORE = 1;
const int ALIGNMENT_MASKED_SCORE = -1;
const char ALIGNMENT_MATCH_REPRESENTATION = '|';
const char ALIGNMENT_MISMATCH_REPRESENTATION = ' ';
const char ALIGNMENT_WOBBLE_REPRESENTATION = ':';
const char ALIGNMENT_MASKED_REPRESENTATION = ' ';

/* Scoring Matrix - Only Takes Normal Nucleotides*/
int ali_score[5][5];
int ali_score5p[5][5];
char ali_rep[5][5];

int bases_are_complementary(char nt1, char nt2) {
	if (nt1 == 'C') return nt2 == 'G';
	if (nt1 == 'G') return nt2 == 'C';
	if (nt1 == 'A') return nt2 == 'U';
	if (nt1 == 'U') return nt2 == 'A';
	return 0;
}

int bases_are_wobble_pair(char nt1, char nt2) {
	if (nt1 == 'G') return nt2 == 'U';
	if (nt1 == 'U') return nt2 == 'G';
	return 0;
}

void initialize_scores() {
	char scored_bases[5] = {'C','G','A','U','X'};
	int i;
	int j;
	for (i = 0; i < sizeof(scored_bases) / sizeof(char); i++) {
		for (j = 0; j < sizeof(scored_bases) / sizeof(char); j++) {
			if (scored_bases[i] == 'X' || scored_bases[j] == 'X') {
				ali_rep[i][j] = ALIGNMENT_MASKED_REPRESENTATION;
				ali_score[i][j] = ALIGNMENT_MASKED_SCORE;
				ali_score5p[i][j] = scale * ALIGNMENT_MASKED_SCORE;
				continue;
			}
			if (bases_are_complementary(scored_bases[i],scored_bases[j])) {
				ali_rep[i][j] = ALIGNMENT_MATCH_REPRESENTATION;
				ali_score[i][j] = ALIGNMENT_MATCH_SCORE;
				ali_score5p[i][j] = scale * ALIGNMENT_MATCH_SCORE;
				continue;
			}
			if (bases_are_wobble_pair(scored_bases[i],scored_bases[j])) {
				ali_rep[i][j] = ALIGNMENT_WOBBLE_REPRESENTATION;
				ali_score[i][j] = ALIGNMENT_WOBBLE_SCORE;
				ali_score5p[i][j] = scale * ALIGNMENT_WOBBLE_SCORE;
				continue;
			}
			ali_rep[i][j] = ALIGNMENT_MISMATCH_REPRESENTATION;
			ali_score[i][j] = ALIGNMENT_MISMATCH_SCORE;
			ali_score5p[i][j] = scale * ALIGNMENT_MISMATCH_SCORE;
		}
	}
}

int score(char nt1, char nt2) {
	return ali_score[bases[(int)nt1]][bases[(int)nt2]];
}

int score5p(char nt1, char nt2) {
	return ali_score5p[bases[(int)nt1]][bases[(int)nt2]];
}

void get_nt_nt_seq_scores(int** nt_nt_scores, char* sequence1, char* sequence2, int query_length, int reference_length) {
	int i, j;
	if (query_length <= 3) {
		/* all positions have 0 nt-nt score */
		for (i = 1; i<= query_length; i++) {
			for (j = 1; j <= reference_length; j++) {
				nt_nt_scores[i][j] = 0;
			}
		}
		return;
	}
	for (j = 1; j <= reference_length; j++) {
		nt_nt_scores[1][j] = 0;	/* scores of positions 1, n-1 and n in the miRNA are scored as 0*/
		nt_nt_scores[2][j] = 0; /* note that the miRNA is reversed in sequence1, so this is positions 1 2 and n*/
		nt_nt_scores[query_length][j] = 0;
	}
	for (i = 3; i <= query_length - 1; i++) {
		for (j = 1; j <= reference_length; j++) {
			if (i <= length_3p_for_weighting || i == query_length) {
				nt_nt_scores[i][j] = score(sequence1[i - 1], sequence2[j - 1]);
			} else {
				nt_nt_scores[i][j] = score5p(sequence1[i - 1], sequence2[j - 1]);
			}
		}
	}
}

void nullify_overlaps(int total_elements, score_struct* scores, int overlap_limit, int *scores_length) {
	int index;
	int outer_ref_end;
	int inner_index;
	for (index = 0; index <= total_elements; index++) {
		if (scores[index].score == 0) {continue;}
		(*scores_length)++;
		/* use query end to figure where ref offset relative to fictional query offset 0 would be */
		outer_ref_end = scores[index].reference_trace_end - scores[index].query_trace_end;
		for (inner_index = total_elements; inner_index > index; inner_index--) {
			int inner_ref_end = scores[inner_index].reference_trace_end -
					scores[inner_index].query_trace_end;
			int diff;
			diff = abs(inner_ref_end - outer_ref_end);
			if (diff <= overlap_limit) {
				scores[inner_index].score = 0;
			}
		}
	}
}

int max_finder_and_track_fourstates(int a, int b, int c, int *track) {
	int max;
	if (a >= b && a >= c) {
		max = a;
		*track = 1;
	} else if (b >= c) {
		max = b;
		*track = 2;
	} else {
		max = c;
		*track = 3;
	}
	return max;
}

int max_finder_and_track_threestates(int a, int b, int c, int *track) {
	int max;
	if (a >= b && a >= c) {
		max = a;
		*track = 1;
	} else if (b >= c) {
		max = b;
		*track = 2;
	} else {
		max = c;
		*track = 3;
	}
	if (max <= 0) {
		max = 0;
		*track = 0;
	}
	return max;
}

void build_matrix(int** best, int*** track, int** a_nt_nt, int** b_gap_nt, int** c_nt_gap, int** nt_nt_score,
		char* sequence1, char* sequence2, int query_length, int reference_length, score_struct* scores, int *scores_length) {
	int i, j, max, in_seed, open, ext, good_count, tmp_a, tmp_b, tmp_c;
	good_count = 0;

	for (i = 1; i <= query_length; i++) {
		if (i <= length_3p_for_weighting || i == query_length) {
			open = gap_open;
			ext = gap_extend;
			in_seed = 0;
		} else {
			open = scale * gap_open;
			ext = scale * gap_extend;
			in_seed = 1;
		}
		for (j = 1; j <= reference_length; j++) {
			max = max_finder_and_track_fourstates(a_nt_nt[i-1][j-1], b_gap_nt[i-1][j-1], c_nt_gap[i-1][j-1], &(track[1][i][j]));
			a_nt_nt[i][j] = max + nt_nt_score[i][j];
			if (a_nt_nt[i][j] <= 0){
				a_nt_nt[i][j] = 0;
				track[1][i][j] = 0;
			}
			tmp_a = (a_nt_nt[i][j-1] + open);
			tmp_b = (b_gap_nt[i][j-1] + ext);
			if (tmp_b > tmp_a) {
				b_gap_nt[i][j] = tmp_b;
				track[2][i][j] = 2;
			} else {
				b_gap_nt[i][j] = tmp_a;
				track[2][i][j] = 1;
			}
			if (!in_seed) {
				tmp_a = (a_nt_nt[i-1][j] + open);
				tmp_c = (c_nt_gap[i-1][j] + ext);
				if (tmp_c > tmp_a) {
					c_nt_gap[i][j] = tmp_c;
					track[3][i][j] = 3;
				} else {
					c_nt_gap[i][j] = tmp_a;
					track[3][i][j] = 1;
				}
			} else {
				c_nt_gap[i][j] = -1;
				track[3][i][j] = 1;
			}
			best[i][j] = max_finder_and_track_threestates(a_nt_nt[i][j], b_gap_nt[i][j], c_nt_gap[i][j], &(track[0][i][j]));
			if (best[i][j] >= score_threshold) {
				/* reject alignments ending with an unaligned reference base */
				if (track[0][i][j] != 2) {
					scores[good_count].score = best[i][j];
					scores[good_count].path = good_count;
					scores[good_count].query_trace_end = i;
					scores[good_count].reference_trace_end = j;
					good_count++;
				}
			}
		}
	}
	qsort(scores, good_count, sizeof(score_struct), cmpscores);
	nullify_overlaps(good_count - 1, scores, 6, scores_length);
	qsort(scores, good_count, sizeof(score_struct), cmpscores);

}

void traceback(int** best, int*** track, char* sequence1, char* sequence2, int i, int j,
		hit_struct* hit_ptr, double hit_score) {
	int length, track_array;
	length = 0;
	hit_ptr->query_end = i;
	hit_ptr->ref_end = j;
	hit_ptr->score = hit_score;
	track_array = track[0][i][j];
	while (best[i][j] > 0) {
		if (track_array == 1) {
			length++;
			hit_ptr->alignment[0][length - 1] = sequence1[i - 1];
			hit_ptr->alignment[2][length - 1] = sequence2[j - 1];
			hit_ptr->alignment[1][length - 1] = ali_rep[(int)bases[(int)sequence1[i - 1]]][bases[(int)sequence2[j - 1]]];
			track_array = track[track_array][i][j];
			i--;
			j--;
		} else if (track_array == 2) {
			length++;
			hit_ptr->alignment[0][length - 1] = '-';
			hit_ptr->alignment[2][length - 1] = sequence2[j - 1];
			hit_ptr->alignment[1][length - 1] = ' ';
			track_array = track[track_array][i][j];
			j--;
		} else if (track_array == 3) {
			length++;
			hit_ptr->alignment[0][length - 1] = sequence1[i - 1];
			hit_ptr->alignment[2][length - 1] = '-';
			hit_ptr->alignment[1][length - 1] = ' ';
			track_array = track[track_array][i][j];
			i--;
		} else {
			hit_ptr->alignment[0][length] = '\0';
			hit_ptr->alignment[1][length] = '\0';
			hit_ptr->alignment[2][length] = '\0';
			break;
		}
	}
	hit_ptr->query_start = i;
	hit_ptr->ref_start = j;
}

int overlap(int start_a, int end_a, int start_b, int end_b) {
	int max_start = start_a > start_b ? start_a : start_b;
	int min_end = end_a < end_b ? end_a : end_b;
	return min_end - max_start + 1;
}

int testfor_overlap(int* good_ones_starts_j, int* good_ones_ends_j, int* good_ones_count, int test_start, int test_end) {
	int index, good_call;
	good_call = 1;
	if (*good_ones_count < 0) {
		good_call = 1;
		*good_ones_count = *good_ones_count + 1;
		return good_call;
	}
	for (index = 0; index <= *good_ones_count; index++) {
		if (overlap(good_ones_starts_j[index],good_ones_ends_j[index],test_start,test_end) > 6) {
			good_call = 0;
			break;
		}
	}
	if (good_call == 1) {*good_ones_count = *good_ones_count + 1;}
	/*
	 printf("ref_s ref_e test_s test_e min_end max_start good_call good_ones_count %d %d %d %d %d %d %d %d\n",
			 ref_start, ref_end, test_start, test_end, min_end, max_start, good_call, *good_ones_count);
	 */
	return good_call;
}
