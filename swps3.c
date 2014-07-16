/** \file swps3.c
 *
 * Main procedure and multi-threading code.
 */
/*
 * Copyright (c) 2007-2008 ETH Zürich, Institute of Computational Science
 *
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following
 * conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

#include "swps3.h"
#include "matrix.h"
#include "fasta.h"
#include "DynProgr_scalar.h"
#if defined(__SSE2__)
#include "DynProgr_sse_byte.h"
#include "DynProgr_sse_short.h"
#endif
#if defined(__ALTIVEC__)
#include "DynProgr_altivec.h"
#endif
#if defined(__PS3)
#include "DynProgr_PPU.h"
#endif
#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <sys/select.h>
#include <unistd.h>
#include <string.h>
#include <assert.h>
#include <float.h>

void normalizeSequence(char* seq, int seqLen) {
	int i;
	for (i = 0; i < seqLen; ++i) {
		seq[i] -= 'A';
	}
}

int main( int argc, char * argv[] ){

	srand(time(NULL));

	int8_t bMatrix[MATRIX_DIM*MATRIX_DIM] __ALIGNED__;
	int16_t sMatrix[MATRIX_DIM*MATRIX_DIM] __ALIGNED__;

	int i, j, k;

	double gapOpen = -20.0;
	double gapExtend = -2.0;

	char* query;
	char* db;

	for(i = 0; i < 1000000; ++i) {

		for (j = 0; j < MATRIX_DIM; ++j) {
			for (k = 0; k < MATRIX_DIM; ++k) {
				int8_t val = (rand() % 256 - 128) / 30;
				sMatrix[k + j * MATRIX_DIM] = val ;
				bMatrix[k + j * MATRIX_DIM] = val;
			}
		}

		int ql = rand() % 5000;
		int dbl = rand() % 5000;

		query = (char*) malloc((ql + 1) * sizeof(char));
		db = (char*) malloc((dbl + 1) * sizeof(char));

		for (k = 0; k < ql; ++k) {
			query[k] = 'A' + rand() % 26;
		}

		for (k = 0; k < dbl; ++k) {
			db[k] = 'A' + rand() % 26;
		}

		query[ql] = '\0';
		db[dbl] = '\0';

		/*printf("Aligning: %s\n to: %s\n", query, db);*/

		int queryLen = strlen(query);
		int dbLen = strlen(db);
		normalizeSequence(query, queryLen);
		normalizeSequence(db, dbLen);
		Options options;
		options.threshold = 120.0;
		options.gapOpen = gapOpen;
		options.gapExt = gapExtend;

		ProfileByte  * profileByte = swps3_createProfileByteSSE( query, queryLen, bMatrix );
		ProfileShort * profileShort = swps3_createProfileShortSSE( query, queryLen, sMatrix );

		double byteScore = swps3_alignmentByteSSE( profileByte, db, dbLen, &options);
		double shortScore = swps3_alignmentShortSSE( profileShort, db, dbLen, &options );

		if(byteScore != shortScore && byteScore != DBL_MAX) {
			printf("ERROR!");
			exit(-1);
		}

		/*printf("SSE score BYTE: %f\n", byteScore);
		printf("SSE score SHORT %f\n", shortScore);*/
		/*printf("Scalar score: %f\n\n", swps3_alignScalar( dmatrix, query[i], queryLen, db[i], dbLen, &options));*/

		swps3_freeProfileByteSSE(profileByte);
		swps3_freeProfileShortSSE(profileShort);
		free(query);
		free(db);

	}

	return 0;
}
