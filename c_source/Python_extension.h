/*
 * Python_extension.h
 *
 *  Created on: Jul 8, 2014
 *      Author: machine
 */

#ifndef PYTHON_EXTENSION_H_
#define PYTHON_EXTENSION_H_

#include "matrix.h"
#include "DynProgr_sse_byte.h"
#include "DynProgr_sse_short.h"

#define MAXSEQLEN	100010
#define MINUSINF (-999999999)
#define MAXMUTDIM       130

#define MMAX(a,b) ((a)>(b)?(a):(b))
/*#define ALIGN_TO(al, typ, ptr) ptr = (typ*)(((unsigned long)ptr)/(al)*(al)+(al));*/

/*typedef union {
	double* ps;
	char* ds;
} SeqCont;*/

extern double coldel[MAXSEQLEN+1], S[MAXSEQLEN+1];
extern int DelFrom[MAXSEQLEN+1];

double c_align_scalar_reference_local(double* matrix, const char *s1, int ls1, const char *s2, int ls2, double gap_open, double gap_ext, double threshold);
double c_align_profile_byte_sse_local(ProfileByte* profile, const char *s2, int ls2, double gap_open, double gap_ext, double threshold);
double c_align_profile_short_sse_local(ProfileShort* profile, const char *s2, int ls2, double gap_open, double gap_ext, double threshold);

ProfileByte* c_create_profile_byte_sse_local(const char* query, int query_len, signed char* matrix);
ProfileShort* c_create_profile_short_sse_local(const char* query, int query_len, signed short* matrix);

void c_free_profile_byte_sse_local(ProfileByte* profile);
void c_free_profile_short_sse_local(ProfileShort* profile);

double c_align_double_global(double* matrix, const char *s1, int ls1, const char *s2, int ls2, double gap_open, double gap_ext);
int c_align_strings(double* matrix, char *s1, int len1, char *s2, int len2, double escore, char *o1, char *o2, double maxerr, double gap_open, double gap_ext);

#endif /* PYTHON_EXTENSION_H_ */
