/*
 * Python_extension.h
 *
 *  Created on: Jul 8, 2014
 *      Author: machine
 */

#ifndef PYTHON_EXTENSION_H_
#define PYTHON_EXTENSION_H_

#include "matrix.h"
#include "DynProgr_scalar.h"

struct {
	double* matrix;

} AlignmentEnvironment;

/* Accessible from Python */
double python_alignScalar(double* matrix, const char *s1, int ls1, const char *s2, int ls2, double gapOpen, double gapExt, double threshold);
double python_alignScalarConvertMatrix(double** matrix, const char* columns, int columnLen, const char *s1, int ls1, const char *s2, int ls2, double gapOpen, double gapExt, double threshold);

/* Not accessible from Python */
double** convertMatrix(double** matrix, const char* columns, int columnLen);

#endif /* PYTHON_EXTENSION_H_ */
