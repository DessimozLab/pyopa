/*
 * Python_extension.c
 *
 *  Created on: Jul 8, 2014
 *      Author: machine
 */

#include <stdlib.h>
#include <stdio.h>

#include "Python_extension.h"
#include "swps3.h"


double python_alignScalar(double* matrix, const char *s1, int ls1, const char *s2, int ls2, double gapOpen, double gapExt, double threshold) {
	Options options;
	options.gapOpen = gapOpen;
	options.gapExt = gapExt;
	options.threshold = threshold;

	return swps3_alignScalar(matrix, s1, ls1, s2, ls2, &options);
}

double python_alignScalarConvertMatrix(double** matrix, const char* columns, int columnLen, const char *s1, int ls1, const char *s2, int ls2, double gapOpen, double gapExt, double threshold) {
	double* converted = convertMatrix(matrix, columns, columnLen);
	double ret = python_alignScalar(matrix, s1, ls1, s2, ls2, gapOpen, gapExt, threshold);

	free(converted);

	return ret;
}

double** convertMatrix(double** matrix, const char* columns, int columnLen) {
	int conLen = columnLen * columnLen;
	int i, j;
	double** converted = (double*) malloc(sizeof(double) * conLen);

	memset(converted, 0.0, conLen);

	for (i = 0; i < columnLen; ++i) {
		for (j = 0; j < columnLen; ++j) {
			/* sbmatrix[ (int)leftChar ][ (int)topChar[ i ] ] --> from matrix.c / swps3_readSBMatrix*/
			converted[(int)(columns[i] - 'A')][(int)(columns[j] - 'A')] = matrix[i][j];
		}
	}

	return converted;
}
