#include <math.h>
#include "EstimatePamEigen.h"

/*
 # EstimatePam (see EstimatePam.c) with two changes that make it faster:
 #  - the aligned pairs of symbols are counted, so scores and likelihoods
 #    are sums over the (at most MATRIX_DIM^2) distinct pairs instead of
 #    over all positions of the alignment
 #  - exp(p*logPAM1) and its derivatives are evaluated from the
 #    eigendecomposition of logPAM1, only for the pairs that occur,
 #    instead of computing matrix exponentials and products in every
 #    Newton iteration
 */

#define NPAIRS (MATRIX_DIM * MATRIX_DIM)

typedef struct PairCounts {
	int delF, /* number of indels */
	delI, /* number of amino acids deleted - delF */
	nmatch, /* number of aligned pairs of symbols */
	npairs, /* number of distinct pairs */
	idx[NPAIRS]; /* distinct pairs, as index into a Dayhoff matrix */
	double cnt[NPAIRS]; /* number of occurrences of each pair */
} PairCounts;

static void CountPairs(const char *o1, const char *o2, int len, PairCounts *pc) {
	int counts[NPAIRS] = { 0 };
	int i, ind;

	pc->delF = pc->delI = pc->nmatch = pc->npairs = 0;
	for (ind = i = 0; i < len; i++) {
		if (o1[i] == '_') {
			if (o2[i] == '_')
				continue; /* do not count del==del */
			if (ind == 1)
				pc->delI++;
			else {
				pc->delF++;
				ind = 1;
			}
		} else if (o2[i] == '_') {
			if (ind == 2)
				pc->delI++;
			else {
				pc->delF++;
				ind = 2;
			}
		} else {
			if (o1[i] != 'X' && o2[i] != 'X') {
				counts[o1[i] * MATRIX_DIM + o2[i]]++;
				pc->nmatch++;
			}
			ind = 0;
		}
	}
	for (i = 0; i < NPAIRS; i++)
		if (counts[i] > 0) {
			pc->idx[pc->npairs] = i;
			pc->cnt[pc->npairs++] = counts[i];
		}
}

static double ComputeScorePairs(int idms, const PairCounts *pc, const DayMatrix *DMS) {
	int k;
	const double *Simi = DMS[idms].Simi;
	double cj = pc->delF * DMS[idms].FixedDel + pc->delI * DMS[idms].IncDel;

	for (k = 0; k < pc->npairs; k++)
		cj += pc->cnt[k] * Simi[pc->idx[k]];
	return cj;
}

/* Multi-point multiple maximum searching a la Brent, as FindMaxPamR */
static int FindMaxPamPairs(int x[], double f[], int n, double *res,
		const PairCounts *pc, const DayMatrix *DMS) {
	int j, j1, j2, x2[3 * MAXPOINTS + 1];
	double cj, c1, c2, f2[3 * MAXPOINTS + 1];

	restart:
	/* end condition */
	for (j = 1; j < n && x[j] - x[j - 1] <= 1; j++)
		;
	if (j >= n) {
		for (j1 = j = 0; j < n; j++)
			if (f[j] > f[j1])
				j1 = j;
		*res = f[j1];
		return (x[j1]);
	}

	/* can insert new points */
	if (n < MAXPOINTS) {
		for (j1 = j = 1; j < n; j++)
			if (x[j] - x[j - 1] > x[j1] - x[j1 - 1])
				j1 = j;
		j = (x[j1] + x[j1 - 1]) / 2;
		cj = ComputeScorePairs(j, pc, DMS);
		for (j2 = n - 1; j2 >= j1; j2--) {
			f[j2 + 1] = f[j2];
			x[j2 + 1] = x[j2];
		}
		x[j1] = j;
		f[j1] = cj;
		n++;
		goto restart;
	}

	/* if a valley is found, split into two searches */
	for (j = 1; j < n - 1; j++)
		if (f[j] < f[j - 1] && f[j] < f[j + 1]) {
			for (j2 = j; j2 < n; j2++) {
				x2[j2 - j] = x[j2];
				f2[j2 - j] = f[j2];
			}
			j2 = FindMaxPamPairs(x2, f2, n - j, &c2, pc, DMS);
			j1 = FindMaxPamPairs(x, f, j + 1, &c1, pc, DMS);
			if (c1 > c2) {
				*res = c1;
				return (j1);
			} else {
				*res = c2;
				return (j2);
			}
		}

	/* discard end farthest from maximum */
	for (j1 = j = 0; j < n; j++)
		if (f[j] > f[j1])
			j1 = j;
	if (x[j1] - x[0] > x[n - 1] - x[j1])
		for (j = 1; j < n; j++) {
			f[j - 1] = f[j];
			x[j - 1] = x[j];
		}
	n--;
	goto restart;
}

/* PamVariance by maximum likelihood, as ComputeVarianceML */
static double VarianceML(const PairCounts *pc, int dmsnumb, const DayMatrix *DMS,
		int DMSLen, const PamModel *m) {
	double e[MATRIX_DIM], le[MATRIX_DIM], l2e[MATRIX_DIM];
	double FD1, FD1d, p, logL1, logL2, incr, x1, Mp, logMp, log2Mp, vw, rmpam;
	const double *Va, *Wb;
	int j, k, l;

	for (j = DMSLen; j > 3 && DMS[j - 1].FixedDel >= -1; j--)
		;
	FD1 = log(10.0) * (DMS[1].FixedDel - DMS[j - 1].FixedDel)
			/ log(DMS[1].PamNumber / DMS[j - 1].PamNumber);
	FD1d = pc->delF * FD1 / 10;
	p = DMS[dmsnumb].PamNumber;

	for (;;) {
		/* exp(p*logPAM1), logPAM1*exp(p*logPAM1) and logPAM1^2*exp(p*logPAM1)
		 are V diag(f(lam)) W with f = exp(p*lam), lam*exp(p*lam), lam^2*exp(p*lam) */
		for (l = 0; l < m->nlam; l++) {
			e[l] = exp(p * m->lam[l]);
			le[l] = m->lam[l] * e[l];
			l2e[l] = m->lam[l] * le[l];
		}
		logL1 = FD1d / p;
		logL2 = -FD1d / (p * p);
		for (k = 0; k < pc->npairs; k++) {
			Va = m->V + (pc->idx[k] / MATRIX_DIM) * MATRIX_DIM;
			Wb = m->W + pc->idx[k] % MATRIX_DIM;
			Mp = logMp = log2Mp = 0;
			for (l = 0; l < m->nlam; l++) {
				vw = Va[l] * Wb[l * MATRIX_DIM];
				Mp += vw * e[l];
				logMp += vw * le[l];
				log2Mp += vw * l2e[l];
			}
			/* symbols not in the model (e.g. X), ignore entry */
			if (Mp == 0 && logMp == 0 && log2Mp == 0)
				continue;
			x1 = logMp / Mp;
			logL1 += pc->cnt[k] * x1;
			logL2 += pc->cnt[k] * (log2Mp / Mp - x1 * x1);
		}

		if (logL2 >= 0)
			break;
		incr = -logL1 / logL2;
		/* do not go beyond the limit of the DMS array */
		if (p + incr > DMS[DMSLen - 1].PamNumber)
			incr = DMS[DMSLen - 1].PamNumber - p;
		x1 = 0.04 * (p + 25);
		if (fabs(incr) > x1)
			incr = (incr < 0) ? -x1 : x1;
		p += incr;
		if (fabs(incr) < 1e-8 * p)
			break;
	}
	if (logL2 >= 0 && p >= DMS[DMSLen - 1].PamNumber) {
		/* at the end of the DMS, estimate the variance by the idealized
		 mutations model */
		p = DMS[DMSLen - 1].PamNumber;
		rmpam = pow(1 - 0.20 / 19, -p);
		logL2 = -pc->nmatch * 19 * log(1 - 0.20 / 19) * log(1 - 0.20 / 19)
				/ ((rmpam - 1) * (rmpam + 19));
	}
	return -1 / logL2;
}

/* PamVariance for very small pam, as ComputeSmallVariance */
static double VarianceSmall(const PairCounts *pc, int dmsnumb, const DayMatrix *DMS,
		const PamModel *m) {
	double a1, a2, b0, c, t1, t2, t3, t4, t5, var;
	const double *L1 = m->L1, *L2 = m->L2, *L3 = m->L3;
	int ij, k, n = MATRIX_DIM;

	a1 = a2 = 0;
	b0 = pc->delF * 0.7434;
	for (k = 0; k < pc->npairs; k++) {
		ij = pc->idx[k];
		c = pc->cnt[k];
		if (L1[ij] == 0)
			continue;
		if (ij % (n + 1) == 0) {
			a1 += c * L1[ij];
			a2 += c * 0.5 * (L2[ij] - L1[ij] * L1[ij]);
		} else {
			b0 += c;
			a1 += c * 0.5 * L2[ij] / L1[ij];
			a2 += c * (L3[ij] / 6 - L2[ij] * L2[ij] / 8 / L1[ij]) / L1[ij];
		}
	}

	t5 = b0 * b0;
	t3 = a1 * a1;
	t2 = 2.0 * t3 + a2 * t5;
	t4 = b0 + 1.0;
	t1 = (3.0 * b0 + 2.0) * a2 + t2;
	var = (4.0 * t3 * t3
			+ ((32.0 + 4.0 * t5 + 24.0 * b0) * t3
					+ (28.0 * b0 + 12.0 + (23.0 + 8.0 * b0 + t5) * t5) * a2)
					* a2) / t3 * t4 / (t1 * t1);
	if (a1 == 0 || t1 == 0)
		var = DMS[dmsnumb].PamNumber * DMS[dmsnumb].PamNumber;
	return var;
}

void EstimatePamEigen(const char* o1, const char* o2, int len, DayMatrix* DMS,
		int DMSLen, const PamModel* model, double* result) {
	PairCounts pc;
	double f[3 * MAXPOINTS + 1], sim;
	int k, x[3 * MAXPOINTS + 1], dmsnumb;

	/* Increase the length to be compatible with the core */
	DMSLen += 1;

	CountPairs(o1, o2, len, &pc);
	for (k = 0; k <= 3 * MAXPOINTS; k++) {
		x[k] = k * (DMSLen - 2) / (3 * MAXPOINTS) + 1;
		f[k] = ComputeScorePairs(x[k], &pc, DMS);
	}
	dmsnumb = FindMaxPamPairs(x, f, 3 * MAXPOINTS + 1, &sim, &pc, DMS);

	result[0] = ComputeScorePairs(dmsnumb, &pc, DMS);
	result[1] = DMS[dmsnumb].PamNumber;
	if (dmsnumb <= 1)
		result[2] = VarianceSmall(&pc, dmsnumb, DMS, model);
	else
		result[2] = VarianceML(&pc, dmsnumb, DMS, DMSLen, model);
}
