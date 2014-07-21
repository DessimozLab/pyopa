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

/* Accessible from Python */
/*double python_alignByte(BMatrix matrix, const char *s1, int ls1, const char *s2, int ls2, double gapOpen, double gapExt, double threshold);
double python_alignShort(SMatrix matrix, const char *s1, int ls1, const char *s2, int ls2, double gapOpen, double gapExt, double threshold);*/

double python_alignScalar(DMatrix matrix, const char *s1, int ls1, const char *s2, int ls2, double gapOpen, double gapExt, double threshold);
double python_alignByteProfileSSE(ProfileByte* profile, const char *s2, int ls2, double gapOpen, double gapExt, double threshold);
double python_alignShortProfileSSE(ProfileShort* profile, const char *s2, int ls2, double gapOpen, double gapExt, double threshold);

ProfileByte* python_createByteProfileSSE(const char* query, int queryLen, BMatrix matrix);
ProfileShort* python_createShortProfileSSE(const char* query, int queryLen, SMatrix matrix);

void python_freeProfileByteSSE(ProfileByte* profile);
void python_freeProfileShortSSE(ProfileShort* profile);

#endif /* PYTHON_EXTENSION_H_ */
