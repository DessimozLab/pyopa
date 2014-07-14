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


double python_alignScalar(DMatrix matrix, const char *s1, int ls1, const char *s2, int ls2, double gapOpen, double gapExt, double threshold) {
	Options options;
	options.gapOpen = gapOpen;
	options.gapExt = gapExt;
	options.threshold = threshold;

	if(((unsigned long)matrix) % 16 != 0) {
		printf("Matrix has to be aligned to 16 bytes!\n");
		return 0.0;
	}

#ifdef PY_DEBUG
	printf("SCALAR aligning with matrix (gapExt = %.2f, gapOpen = %.2f, threshold = %.2f, matrix loc = %lu, mod 16 = %lu):\n",
			gapOpen, gapExt, threshold, ((unsigned long)matrix), ((unsigned long)matrix) % 16);
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

double python_alignByte(SBMatrix matrix, const char *s1, int ls1, const char *s2, int ls2, double gapOpen, double gapExt, double threshold) {
	Options options;
	options.gapOpen = gapOpen;
	options.gapExt = gapExt;
	options.threshold = threshold;

	if(((unsigned long)matrix) % 16 != 0) {
		printf("Matrix has to be aligned to 16 bytes!\n");
		return 0.0;
	}

#ifdef PY_DEBUG
	printf("BYTE aligning with matrix (gapExt = %.2f, gapOpen = %.2f, threshold = %.2f, matrix loc = %lu, mod 16 = %lu):\n",
			gapOpen, gapExt, threshold, ((unsigned long)matrix), ((unsigned long)matrix) % 16);
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
double python_alignShort(SBMatrix matrix, const char *s1, int ls1, const char *s2, int ls2, double gapOpen, double gapExt, double threshold) {
	Options options;
	options.gapOpen = gapOpen;
	options.gapExt = gapExt;
	options.threshold = threshold;

	if(((unsigned long)matrix) % 16 != 0) {
		printf("Matrix has to be aligned to 16 bytes!\n");
		return 0.0;
	}

#ifdef PY_DEBUG
	printf("SHORT aligning with matrix (gapExt = %.2f, gapOpen = %.2f, threshold = %.2f, matrix loc = %lu, mod 16 = %lu):\n",
			gapOpen, gapExt, threshold, ((unsigned long)matrix), ((unsigned long)matrix) % 16);
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
}





