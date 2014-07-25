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
#include <sys/time.h>
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
	/* TODO this main function is not intended to work correctly (the matrices are not even initialized) */
	srand(time(NULL));

	int i, j;

	int8_t bMatrix[MATRIX_DIM*MATRIX_DIM] __ALIGNED__;
	int16_t sMatrix[MATRIX_DIM*MATRIX_DIM] __ALIGNED__;

	const char* loc = "python/test/seq.test";

	FILE *fp;
	fp = fopen(loc,"r");

	int ch, number_of_lines = 0;

	do
	{
	    ch = fgetc(fp);
	    if(ch == '\n')
	    	number_of_lines++;
	} while (ch != EOF);

	char** sequences = (char**) malloc(number_of_lines * sizeof(char*));
	int* lens = (int*) malloc(number_of_lines * sizeof(int));

	for (i = 0; i < number_of_lines; ++i) {
		sequences[i] = (char*) malloc(20000 * sizeof(char*));
	}

	fseek (fp, 9, SEEK_SET);

	i = 0;
	while(fgets(sequences[i], 20000, fp) != NULL)
	{
		/*replace new line with end of string*/
		sequences[i][strlen(sequences[i])-1] = '\0';
		lens[i] = strlen(sequences[i]);
		normalizeSequence(sequences[i], lens[i]);
		++i;
	}

	fclose(fp);

	Options options;
	options.threshold = 85.0;
	options.gapOpen = -20;
	options.gapExt = -2;

	double* res = (double*) malloc((number_of_lines - 1) * number_of_lines / 2 * sizeof(double));

	int resSize = 0;

	float millis = 0.0;
	unsigned long start, end;

	start = clock();

	for (i = 0; i < number_of_lines; ++i) {
		ProfileShort* profile = swps3_createProfileShortSSE(sequences[i], lens[i], sMatrix);
		for (j = i + 1; j < number_of_lines; ++j) {
			res[resSize] = swps3_alignmentShortSSE(profile, sequences[j], lens[j], &options);
			++resSize;
		}

		swps3_freeProfileShortSSE(profile);
	}

	end = clock();

	millis = (end - start) / 1000.0;

	printf("A total of %d alignments (all to all with %d) have been done in %.3fs", resSize, number_of_lines, millis/1000.0);


	/*for (i = 0; i < resSize; ++i) {
		printf("%.3f\n", res[i]);
	}*/

	for (i = 0; i < number_of_lines; ++i) {
		free(sequences[i]);
	}
	free(sequences);
	free(res);

	return 0;
}
