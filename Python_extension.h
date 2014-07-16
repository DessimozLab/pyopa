/*
 * Python_extension.h
 *
 *  Created on: Jul 8, 2014
 *      Author: machine
 */

#ifndef PYTHON_EXTENSION_H_
#define PYTHON_EXTENSION_H_

#include "matrix.h"

/* Accessible from Python */
double python_alignScalar(DMatrix matrix, const char *s1, int ls1, const char *s2, int ls2, double gapOpen, double gapExt, double threshold);
double python_alignByte(BMatrix matrix, const char *s1, int ls1, const char *s2, int ls2, double gapOpen, double gapExt, double threshold);
double python_alignShort(SMatrix matrix, const char *s1, int ls1, const char *s2, int ls2, double gapOpen, double gapExt, double threshold);


#endif /* PYTHON_EXTENSION_H_ */
