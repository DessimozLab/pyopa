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

/*
double python_alignScalar(DMatrix matrix, const char *s1, int ls1, const char *s2, int ls2, double gapOpen, double gapExt, double threshold) {
	Options options;
	options.gapOpen = gapOpen;
	options.gapExt = gapExt;
	options.threshold = threshold;

#ifdef PY_DEBUG
	printf("SCALAR aligning with matrix (gapExt = %.2f, gapOpen = %.2f, threshold = %.2f, matrix loc = %lu, mod 16 = %lu):\n",
			gapExt, gapOpen, threshold, ((unsigned long)matrix), ((unsigned long)matrix) % 16);
	int i, j;
	for (i = 0; i < 26; ++i) {
		printf("[");
		for (j = 0; j < 26; ++j) {
			printf("%+.2f, ", matrix[i * 26 + j]);
		}
		printf("]\n");
	}
#endif

	return swps3_alignScalar(matrix, s1, ls1, s2, ls2, &options);
}

double python_alignByte(BMatrix matrix, const char *s1, int ls1, const char *s2, int ls2, double gapOpen, double gapExt, double threshold) {
	Options options;
	options.gapOpen = gapOpen;
	options.gapExt = gapExt;
	options.threshold = threshold;

#ifdef PY_DEBUG
	printf("BYTE aligning with matrix (gapExt = %.2f, gapOpen = %.2f, threshold = %.2f, matrix loc = %lu, mod 16 = %lu):\n",
			gapExt, gapOpen, threshold, ((unsigned long)matrix), ((unsigned long)matrix) % 16);
	int i, j;
	for (i = 0; i < 26; ++i) {
		printf("[");
		for (j = 0; j < 26; ++j) {
			printf("%+d, ", matrix[i * 26 + j]);
		}
		printf("]\n");
	}
#endif

	ProfileByte  * profileByte = swps3_createProfileByteSSE(s1, ls1, matrix);
	double score = swps3_alignmentByteSSE(profileByte, s2, ls2, &options);
	swps3_freeProfileByteSSE(profileByte);

	return score;
}
double python_alignShort(SMatrix matrix, const char *s1, int ls1, const char *s2, int ls2, double gapOpen, double gapExt, double threshold) {
	Options options;
	options.gapOpen = gapOpen;
	options.gapExt = gapExt;
	options.threshold = threshold;

#ifdef PY_DEBUG
	printf("SHORT aligning with matrix (gapExt = %.2f, gapOpen = %.2f, threshold = %.2f, matrix loc = %lu, mod 16 = %lu):\n",
			gapExt, gapOpen, threshold, ((unsigned long)matrix), ((unsigned long)matrix) % 16);
	int i, j;
	for (i = 0; i < 26; ++i) {
		printf("[");
		for (j = 0; j < 26; ++j) {
			printf("%+d, ", matrix[i * 26 + j]);
		}
		printf("]\n");
	}
#endif

	ProfileShort * profileShort = swps3_createProfileShortSSE(s1, ls1, matrix);
	double score = swps3_alignmentShortSSE( profileShort, s2, ls2, &options);
	swps3_freeProfileShortSSE(profileShort);

	return score;
}*/

double python_alignByteProfileSSE(ProfileByte* profile, const char *s2, int ls2, double gapOpen, double gapExt, double threshold) {
	Options options;
	options.gapOpen = gapOpen;
	options.gapExt = gapExt;
	options.threshold = threshold;

#ifdef PY_DEBUG
	printf("Aligning BYTE: profile = %lu, s2 = %s, len(s2) = %d, gapOpen = %f, gapExt = %f, threshold = %f\n",
			(unsigned long)profile, s2, ls2, options.gapOpen, options.gapExt, options.threshold);
#endif

	return swps3_alignmentByteSSE( profile, s2, ls2, &options);
}
double python_alignShortProfileSSE(ProfileShort* profile, const char *s2, int ls2, double gapOpen, double gapExt, double threshold) {
	Options options;
	options.gapOpen = gapOpen;
	options.gapExt = gapExt;
	options.threshold = threshold;

#ifdef PY_DEBUG
	printf("Aligning SHORT: profile = %lu, s2 = %s, len(s2) = %d, gapOpen = %f, gapExt = %f, threshold = %f\n",
			(unsigned long)profile, s2, ls2, options.gapOpen, options.gapExt, options.threshold);
#endif

	return swps3_alignmentShortSSE( profile, s2, ls2, &options);
}

ProfileByte* python_createByteProfileSSE(const char* query, int queryLen, BMatrix matrix) {
	ProfileByte* pb = swps3_createProfileByteSSE(query, queryLen, matrix);
#ifdef PY_DEBUG
	printf("Successfully created BYTE profile at %lu from query = %s, len(query) = %d\n", (unsigned long)pb, query, queryLen);
	printf("The matrix used for the BYTE profile: \n");
	int i, j;
	for (i = 0; i < 26; ++i) {
		printf("[");
		for (j = 0; j < 26; ++j) {
			printf("%+.2f, ", (double) matrix[i * 26 + j]);
		}
		printf("]\n");
	}
#endif
	return pb;
}
ProfileShort* python_createShortProfileSSE(const char* query, int queryLen, SMatrix matrix) {
	ProfileShort* ps = swps3_createProfileShortSSE(query, queryLen, matrix);
#ifdef PY_DEBUG
	printf("Successfully created SHORT profile at %lu from query = %s, len(query) = %d\n", (unsigned long)ps, query, queryLen);
	printf("The matrix used for the SHORT profile: \n");
	int i, j;
	for (i = 0; i < 26; ++i) {
		printf("[");
		for (j = 0; j < 26; ++j) {
			printf("%+.2f, ", (double) matrix[i * 26 + j]);
		}
		printf("]\n");
	}
#endif
	return ps;
}

void python_freeProfileByteSSE(ProfileByte* profile) {
#ifdef PY_DEBUG
	printf("Freeing BYTE profile at %lu\n", (unsigned long)profile);
#endif
	swps3_freeProfileByteSSE(profile);
}
void python_freeProfileShortSSE(ProfileShort* profile) {
#ifdef PY_DEBUG
	printf("Freeing SHORT profile at %lu\n", (unsigned long)profile);
#endif
	swps3_freeProfileShortSSE(profile);
}

double python_alignScalar(DMatrix matrix, const char *s1, int ls1, const char *s2, int ls2, double gapOpen, double gapExt, double threshold) {
	Options options;
	options.gapOpen = gapOpen;
	options.gapExt = gapExt;
	options.threshold = threshold;

	return swps3_alignScalar(matrix, s1, ls1, s2, ls2, &options);
}

