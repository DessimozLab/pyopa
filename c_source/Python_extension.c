/*
 * Python_extension.c
 *
 *  Created on: Jul 8, 2014
 *      Author: machine
 */

#include <stdlib.h>
#include <stdio.h>

#include "Python_extension.h"
#include "DynProgr_scalar.h"
#include "DynProgr_sse_byte.h"
#include "DynProgr_sse_short.h"

char* denormalize(const char* str, int len) {
	char* ret = (char*) malloc((len + 1)* sizeof(char));
	int i;

	for (i = 0; i < len; ++i) {
		ret[i] = 'A' + str[i];
	}

	ret[len] = '\0';

	return ret;
}

void debug_alignment(char* name, void* profile, const char* s2, int ls2, Options* options) {
	char* dns2 = denormalize(s2, ls2);
	printf("Aligning %s: profile = %lu, s2 = %s, len(s2) = %d, gapOpen = %f, gapExt = %f, threshold = %f\n",
				name, (unsigned long)profile, dns2, ls2, options->gapOpen, options->gapExt, options->threshold);
	free(dns2);
}

void debug_profile(char* name, void* pb, const char* q, int queryLen, void* matrix, int isShort) {
	char* query = denormalize(q, queryLen);

	printf("Successfully created %s profile at %lu from query = %s, len(query) = %d\n", name, (unsigned long)pb, query, queryLen);
	printf("The matrix used for the %s profile: \n", name);
	int i, j;
	for (i = 0; i < 26; ++i) {
		printf("[");
		for (j = 0; j < 26; ++j) {
			double val;
			if(isShort) {
				val = (double) ((SMatrix)matrix)[i * 26 + j];
			} else {
				val = (double) ((BMatrix)matrix)[i * 26 + j];
			}

			printf("%+.2f, ", val);
		}
		printf("]\n");
	}

	free(query);
}

double c_align_profile_byte_sse(ProfileByte* profile, const char *s2, int ls2, double gapOpen, double gapExt, double threshold) {
	Options options;
	options.gapOpen = gapOpen;
	options.gapExt = gapExt;
	options.threshold = threshold;

#ifdef PY_DEBUG
	debug_alignment("BYTE", profile, s2, ls2, &options);
#endif

	return swps3_alignmentByteSSE( profile, s2, ls2, &options);
}
double c_align_profile_short_sse(ProfileShort* profile, const char *s2, int ls2, double gapOpen, double gapExt, double threshold) {
	Options options;
	options.gapOpen = gapOpen;
	options.gapExt = gapExt;
	options.threshold = threshold;

#ifdef PY_DEBUG
	debug_alignment("SHORT", profile, s2, ls2, &options);
#endif

	return swps3_alignmentShortSSE( profile, s2, ls2, &options);
}

ProfileByte* c_create_profile_byte_sse(const char* query, int queryLen, BMatrix matrix) {
	ProfileByte* pb = swps3_createProfileByteSSE(query, queryLen, matrix);
#ifdef PY_DEBUG
	debug_profile("BYTE", pb, query, queryLen, matrix, 0);
#endif
	return pb;
}
ProfileShort* c_create_profile_short_sse(const char* query, int queryLen, SMatrix matrix) {
	ProfileShort* ps = swps3_createProfileShortSSE(query, queryLen, matrix);
#ifdef PY_DEBUG
	debug_profile("SHORT", ps, query, queryLen, matrix, 1);
#endif
	return ps;
}

void c_free_profile_byte_sse(ProfileByte* profile) {
#ifdef PY_DEBUG
	printf("Freeing BYTE profile at %lu\n", (unsigned long)profile);
#endif
	swps3_freeProfileByteSSE(profile);
}
void c_free_profile_short_sse(ProfileShort* profile) {
#ifdef PY_DEBUG
	printf("Freeing SHORT profile at %lu\n", (unsigned long)profile);
#endif
	swps3_freeProfileShortSSE(profile);
}

double c_align_scalar(DMatrix matrix, const char *s1, int ls1, const char *s2, int ls2, double gapOpen, double gapExt, double threshold) {
	Options options;
	options.gapOpen = gapOpen;
	options.gapExt = gapExt;
	options.threshold = threshold;

	return swps3_alignScalar(matrix, s1, ls1, s2, ls2, &options);
}

