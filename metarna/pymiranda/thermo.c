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
#include <string.h>
#include "fold_vars.h"
#include "utils.h"
#include "fold.h"
#include "miranda.h"

extern int temperature_override;

double vfold(char* sequence) {
	void* struct1;
	double e1;
	struct1 = (char*)space(sizeof(char) * (strlen(sequence) + 1));
	temperature = temperature_override;
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
	for (i = 0; i < (int)strlen(hit->rest[0]); i++) {
		foldsequence[j] = hit->rest[0][i];
		j++;
	}
	for (i = 0; i < (int)strlen(hit->alignment[0]); i++) {
		if (hit->alignment[0][i] != '-') {
			foldsequence[j] = hit->alignment[0][i];
			j++;
		}
	}
	for (i = 0; i < (int)strlen(hit->rest[3]); i++) {
		foldsequence[j] = hit->rest[3][i];
		j++;
	}
	for (i = 0; i < 7; i++) {
		foldsequence[j] = 'X';
		j++;
	}
	for (i = 0; i < (int)strlen(hit->rest[4]); i++) {
		foldsequence[j] = hit->rest[4][i];
		j++;
	}
	for (i = 0; i < (int)strlen(hit->alignment[2]); i++) {
		if (hit->alignment[2][i] != '-') {
			foldsequence[j] = hit->alignment[2][i];
			j++;
		}
	}
	for (i = 0; i < (int)strlen(hit->rest[1]); i++) {
		foldsequence[j] = hit->rest[1][i];
		j++;
	}
	foldsequence[j] = '\0';
	energy = vfold(foldsequence);
	revstring(hit->alignment[0]);
	revstring(hit->rest[1]);
	revstring(hit->rest[4]);
	return energy;
}
