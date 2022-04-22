/*
 * DynProgr_sse_double.h
 *
 *  Created on: Sep 11, 2014
 *      Author: machine
 */

#ifndef DYNPROGR_SSE_DOUBLE_H_
#define DYNPROGR_SSE_DOUBLE_H_

#include "swps3.h"
#include "matrix.h"
#include "simde/x86/sse2.h"

typedef struct{
	int ls1;
	double* profile;
	double* old_opt;
	double* new_opt;
	double* new_rd;
} ProfileDouble;

double align_double_local(ProfileDouble* profileDouble, const char *s2, int ls2, double gap_open,
		double gap_ext, double threshold, int* max1, int* max2);
ProfileDouble* createProfileDoubleSSE(const char* s1, int ls1, double* matrix);
void free_profile_double_sse(ProfileDouble* profile);

#endif /* DYNPROGR_SSE_DOUBLE_H_ */
