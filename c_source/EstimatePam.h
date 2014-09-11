/*
 * EstimatePam.h
 *
 *  Created on: Aug 8, 2014
 *      Author: machine
 */

#ifndef ESTIMATEPAM_H_
#define ESTIMATEPAM_H_

#include "Python_extension.h"
#include <stdint.h>

#define MAXPOINTS 7
#define SCALB(x, exp) (x * pow(FLT_RADIX, exp))
#define RoundDM(x,v) (SCALB(rint(SCALB(x,v)),-v))

typedef struct Counters {
	int inds[MAXSEQLEN + 20000], /* encoded pair of symbols, indices in Simi */
	ninds, /* number of entries in inds */
	delF, /* number of indels */
	delI, /* number of amino acids deleted - delF */
	delL[MAXSEQLEN / 10]; /* length of each deletion (ending in 0) */
} Counters;

typedef struct DayMatrix {
	double PamNumber;
	double FixedDel;
	double IncDel;
	/*double Dim;*/
	double* Simi;
} DayMatrix;

extern double NorFre[MAXMUTDIM];

DayMatrix* createDayMatrices(double* gapOpen, double* gapExt,
		double* pamDistances, long long* matrix_pointers, int DMSLen);

void EstimatePam(char* o1, char* o2, int len, DayMatrix* DMS, int DMSLen,
		double* logPAM1, double* result);

void freeDayMatrices(DayMatrix* DMS, int DMSLen);

void CreateOrigDayMatrix(double* log_pam1, double pam, double* new_matrix);

#endif /* ESTIMATEPAM_H_ */
