/* thermo.c */
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

/*
 * This file includes headers and links to the RNAlib library of
 * Ivo Hofackers Vienna Package
 */

#include <config.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "fold_vars.h"
#include "utils.h"
#include "fold.h"
#include "miranda.h"

double vfold(char* sequence) {
	void* struct1;
	double e1;
	struct1 = (char*)space(sizeof(char) * (strlen(sequence) + 1));
	temperature = 30;
	initialize_fold(strlen(sequence));
	e1 = fold(sequence, struct1);
	free_arrays();
	free(struct1);
	return e1;
}

double get_energy(hit_struct* hit) {
	double energy = 0;
	int i = 0;
	int j = 0;
	char foldsequence[5000];
	revstring(hit->alignment[0]);
	revstring(hit->rest[1]);
	revstring(hit->rest[4]);
	foldsequence[0] = '\0';
	for (i = 0; i < strlen(hit->rest[0]); i++) {
		foldsequence[j] = hit->rest[0][i];
		j++;
	}
	for (i = 0; i < strlen(hit->alignment[0]); i++) {
		if (hit->alignment[0][i] != '-') {
			foldsequence[j] = hit->alignment[0][i];
			j++;
		}
	}
	for (i = 0; i < strlen(hit->rest[3]); i++) {
		foldsequence[j] = hit->rest[3][i];
		j++;
	}
	for (i = 0; i < 7; i++) {
		foldsequence[j] = 'X';
		j++;
	}
	for (i = 0; i < strlen(hit->rest[4]); i++) {
		foldsequence[j] = hit->rest[4][i];
		j++;
	}
	for (i = 0; i < strlen(hit->alignment[2]); i++) {
		if (hit->alignment[2][i] != '-') {
			foldsequence[j] = hit->alignment[2][i];
			j++;
		}
	}
	for (i = 0; i < strlen(hit->rest[1]); i++) {
		foldsequence[j] = hit->rest[1][i];
		j++;
	}
	foldsequence[j] = '\0';
	/* printf("FOLD: %s\n", foldsequence);*/
	energy = vfold(foldsequence);
	revstring(hit->alignment[0]);
	revstring(hit->rest[1]);
	revstring(hit->rest[4]);
	return energy;
}
