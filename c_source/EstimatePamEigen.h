#ifndef ESTIMATEPAMEIGEN_H_
#define ESTIMATEPAMEIGEN_H_

#include "EstimatePam.h"

/* logPAM1 and the data derived from it that EstimatePamEigen needs. It only
 depends on logPAM1, so it is computed once (see MutipleAlEnv in pyopa.pyx).
 logPAM1 = V diag(lam) W with W = V^-1, restricted to the nlam letters for
 which logPAM1 has non-zero rows/columns; V and W are zero elsewhere. */
typedef struct PamModel {
	int nlam;
	double lam[MATRIX_DIM];
	double V[MATRIX_DIM * MATRIX_DIM]; /* V[a * MATRIX_DIM + k] */
	double W[MATRIX_DIM * MATRIX_DIM]; /* W[k * MATRIX_DIM + b] */
	double L1[MATRIX_DIM * MATRIX_DIM]; /* logPAM1 */
	double L2[MATRIX_DIM * MATRIX_DIM]; /* logPAM1^2 */
	double L3[MATRIX_DIM * MATRIX_DIM]; /* logPAM1^3 */
} PamModel;

/* Same results as EstimatePam, but evaluates exp(pam * logPAM1) and its
 derivatives from the eigendecomposition, only for the pairs of symbols
 occurring in the alignment. */
void EstimatePamEigen(const char* o1, const char* o2, int len, DayMatrix* DMS,
		int DMSLen, const PamModel* model, double* result);

#endif /* ESTIMATEPAMEIGEN_H_ */
