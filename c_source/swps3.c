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
#include "Python_extension.h"
#include "EstimatePam.h"
#include <string.h>

void normalizeSequence(char* seq, int seqLen) {
	int i;
	for (i = 0; i < seqLen; ++i) {
		if (seq[i] != '_') {
			seq[i] -= 'A';
		}
	}
}

void denormalizeSequence(char* str, int len) {
	int i;

	for (i = 0; i < len; ++i) {
		if (str[i] != '_') {
			str[i] = 'A' + str[i];
		}
	}

	str[len] = '\0';
}

double ceil(double val) {
	short tmp = (short) val;
	if (tmp < val) {
		tmp++;
	}
	return tmp;
}

void readDoubleMatrix(DMatrix matrix, const char* fileLoc, Options* options) {
	char line[40000];
	int i, curr_off, offset = 0;

	FILE* f = fopen(fileLoc, "r");

	fgets(line, 40000, f);
	sscanf(line, "%lf", &(options->gapOpen));

	fgets(line, 40000, f);
	sscanf(line, "%lf", &(options->gapExt));

	/* This is now the PanDistance, because I did not want to modify this Options struct but we need it for tests */
	fgets(line, 40000, f);
	sscanf(line, "%lf", &(options->threshold));

	fgets(line, 40000, f);
	for (i = 0; i < MATRIX_DIM * MATRIX_DIM; ++i) {
		sscanf(line + offset, "%lf%n", &matrix[i], &curr_off);
		offset += curr_off;
	}

	fclose(f);
}

void readShortMatrix(SMatrix matrix, const char* fileLoc, Options* options) {
	double scaleFactor = 65535.0 / options->threshold;
	double dMatrix[MATRIX_DIM * MATRIX_DIM];
	int i;

	readDoubleMatrix(dMatrix, fileLoc, options);
	options->gapOpen = (short) ceil(options->gapOpen * scaleFactor);
	options->gapExt = (short) ceil(options->gapExt * scaleFactor);

	for (i = 0; i < MATRIX_DIM * MATRIX_DIM; ++i) {
		matrix[i] = (short) ceil(dMatrix[i] * scaleFactor);
	}
}

int main(int argc, char * argv[]) {

	int16_t sMatrix[MATRIX_DIM * MATRIX_DIM];
	double dMatrix[MATRIX_DIM * MATRIX_DIM];
	/* 1267 instead of 1266 because in darwin indexing starts from 1 and the code is migrated like that */
	double* doubleMatrices[1267];
	double gapOpenCosts[1267];
	double gapExtCosts[1267];
	double pamDistances[1267];
	double logPAM1Matrix[MATRIX_DIM * MATRIX_DIM];
	int i, j;

	Options tmp;
	readDoubleMatrix(logPAM1Matrix,
			"/home/machine/repos/students/2014_Ferenc_Galko_SWPS3_PY/swps3_python_extended/test/data/matrices/C_compatible/logPAM1.dat",
			&tmp);
	for (i = 1; i < 1267; ++i) {
		doubleMatrices[i] = (double*) malloc(
		MATRIX_DIM * MATRIX_DIM * sizeof(double));
		char name[2000];
		sprintf(name,
				"/home/machine/repos/students/2014_Ferenc_Galko_SWPS3_PY/swps3_python_extended/test/data/matrices/C_compatible/%d.dat",
				i);
		readDoubleMatrix(doubleMatrices[i], name, &tmp);
		gapOpenCosts[i] = tmp.gapOpen;
		gapExtCosts[i] = tmp.gapExt;
		/* threshold contains the pamDistance for now */
		pamDistances[i] = tmp.threshold;
	}

	Options options_short;
	Options options_double;

	readShortMatrix(sMatrix,
			"/home/machine/repos/students/2014_Ferenc_Galko_SWPS3_PY/swps3_python_extended/test/data/matrices/C_compatible/1263.dat",
			&options_short);
	readDoubleMatrix(dMatrix,
			"/home/machine/repos/students/2014_Ferenc_Galko_SWPS3_PY/swps3_python_extended/test/data/matrices/C_compatible/98.dat",
			&options_double);

	options_short.threshold = 306.896691;
	options_double.threshold = DBL_MAX;

	double shortFactor = 65535.0 / options_short.threshold;

	/*char query[] =
	 "GANRAKHVKYWRCKWWVASPKNLLFQTVHEMALLVLPEGEWGTMVALARTGFFLLLAFSMGTMSKKFEGNHHWTWVYPFFMELMAGHVAWVLFNLPGEAIVSLRTGYLQRGREKTFVDG";
	 char db[] =
	 "GANRAKHVKYWRTEANPKTCKWWVASPKSNLLFQTVHIKSEGTYLARNSVSATRDTKKVQDLLSRLQTSEYGLRHIFTDARRNETRTGIEMNALLVLPEGEWGTMVALARTGFFLLLAFSMGTMSKKFEGNHHWTWVYPFFMELMAQLHIFNGHVAWVLFNLPGEAIVSLRTGYLQRGREKTFVDG";
	 */
	char o1[MAXSEQLEN], o2[MAXSEQLEN];
	char query[] = "AAA";
	char db[] = "AAA";
	int ql = strlen(query);
	int dbl = strlen(db);

	normalizeSequence(query, ql);
	normalizeSequence(db, dbl);

	ProfileShort* profile = c_create_profile_short_sse_local(query, ql,
			sMatrix);

	double score_short = c_align_profile_short_sse_local(profile, db, dbl,
			options_short.gapOpen, options_short.gapExt,
			options_short.threshold);
	int max1, max2;
	double short_double = c_align_double_local(dMatrix, query, ql, db, dbl,
			options_double.gapOpen, options_double.gapExt,
			options_double.threshold, &max1, &max2);
	double double_global = c_align_double_global(dMatrix, query, ql, db, dbl,
			options_double.gapOpen, options_double.gapExt);
	double scalar_double = c_align_scalar_reference_local(dMatrix, query, ql,
			db, dbl, options_double.gapOpen, options_double.gapExt,
			options_double.threshold);
	/*  */
	int retLen = c_align_strings(dMatrix, query, ql, db, dbl, scalar_double, o1,
			o2, 0.5e-4, options_double.gapOpen, options_double.gapExt);

	denormalizeSequence(o1, retLen);
	denormalizeSequence(o2, retLen);

	score_short /= shortFactor;

	printf("SHORT score: %.5f\n", score_short);
	printf("DOUBLE score: %.15f\n", short_double);
	printf("SCALAR score: %.15f\n", scalar_double);
	printf("double_global: %.15f\n", double_global);
	printf("Concrete alignment (%d):\nSeq1: %s\nSeq2: %s\n", retLen, o1, o2);

	normalizeSequence(o1, retLen);
	normalizeSequence(o2, retLen);
	DayMatrix* DMS = createDayMatrices(gapOpenCosts, gapExtCosts, pamDistances,
			doubleMatrices, 1266);
	double result[3];

	char t1[100] = "KPL__ELA";
	char t2[100] = "SAVLCEDN";
	int tLen = strlen(t1);
	normalizeSequence(t1, strlen(t1));
	normalizeSequence(t2, strlen(t2));
	EstimatePam(t1, t2, tLen, DMS, 1266, logPAM1Matrix, result);
	printf("EstimatePam:\n%f\n%f\n%f\n", result[0], result[1], result[2]);

	freeDayMatrices(DMS, 1266);
	swps3_freeProfileShortSSE(profile);

	return 0;
}
