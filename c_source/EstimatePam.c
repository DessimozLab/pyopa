#include <math.h>
#include <stdio.h>
#include "EstimatePam.h"
#include <float.h>
#include <string.h>

/*
 #  compute the best Pam index for a given
 #  alignment.  This gives the maximum likelihood estimator
 #  for the pam number or the number of evolutionary steps
 #  separating two sequences
 #				Gaston H. Gonnet (Oct 1990)
 #  (algorithm improved an rewritten in C) (Jan 1991)
 #  (code for multiple maxima) (Dec 2001)
 #  (rearranged to compute indices and to work for any dimension Mar/2005)
 #  (special code for very small pam, Aug/2005)
 #  (ML Distance and Variance, Mar/2006)
 #							*/

double NorFre[MAXMUTDIM];
int DelFrom[MAXSEQLEN+1];

#define ABS(a)	((a)>=0?(a):-(a))

#ifdef __i386
#define ADD_ERR(a,b,c) { volatile double tmpa, tmpb, tmpc, tmpe, res; \
        tmpa = a;  tmpb = b; tmpc = c; res = a; \
        res += tmpb; \
        if( ABS(tmpa) > ABS(tmpb) ) \
             { tmpe = res-tmpa;  tmpc += tmpb-tmpe; } \
        else { tmpe = res-tmpb;  tmpc += tmpa-tmpe; } \
        a = res; c = tmpc; }
#else
#define ADD_ERR(a,b,c) { double tmpa, tmpb, tmpe; \
        tmpa = a;  tmpb = b;\
        a += tmpb; \
        if( ABS(tmpa) > ABS(tmpb) ) \
             { tmpe = a-tmpa;  c += tmpb-tmpe; } \
        else { tmpe = a-tmpb;  c += tmpa-tmpe; } }
#endif

static void mmul(double *AA, double *BB, double *CC, int n)
/* multiply AAxBB -> CC (all are n x n) */

{
	int in, j, k, kn, n2;
	double cij, cij2, *AAin, *BBj;
	if (n > 30 && n <= 10000) {
		/*mmul_inorder(AA, BB, CC, n);*/
		printf("Matrix multiplication error?");
		return;
	}
	n2 = n * n;
	for (in = 0; in < n2; in += n) {
		AAin = AA + in;
		for (j = 0; j < n; j++) {
#ifdef HASLONGDOUBLE
			long double ldr= {0};
			BBj = BB+j;
			for( kn=k=0; k<n; k++,kn+=n ) ldr = ldr + AAin[k] * BBj[kn];
			CheckFinite( CC[in+j] = ldr );
#else
			cij = cij2 = 0;
			BBj = BB + j;
			for (kn = k = 0; k < n; k++, kn += n)
				ADD_ERR(cij, AAin[k] * BBj[kn], cij2);
			CC[in + j] = cij + cij2;
#endif
		}
	}
}

static void mexp(double *x2k, double *t, double *c, double *r, int n)
/* compute exp(matrix)  input (destroyed) is in x2k and
 output is placed in r */
{
	int f, i, in, j, k, m, mn, n2, nonz;
	double cij, norm, norm1, p, *tj, twok, *x2kin;

	n2 = n * n;
	/* Compute norm */
	norm = nonz = 0;
	for (i = 0; i < n; i++) {
		norm1 = 0;
		for (j = 0; j < n; j++)
			if (x2k[i * n + j] != 0) {
				nonz++;
				norm1 += ABS(x2k[i * n + j]);
			}
		if (norm1 > norm)
			norm = norm1;
		/*if( isnan(norm) ) userror("matrix too ill-conditioned");*/
	}
	if (n > 30 && nonz <= 9 * n) {
		/*mexp_sparse( x2k, t, c, r, n, nonz, norm );*/
		printf("mexp error?");
		return;
	}

	/* Compute k, number of powerings to leave norm<0.01 */
	/*k = ilogb(200 * norm);*/
	k = (int) (log(fabs(200 * norm)) / log(FLT_RADIX));
	if (k < 0)
		k = 0;
	/*if( k >= 1023 ) userror("cannot compute exp of a matrix");*/

	/*twok = scalb(1.0, k);
	 norm = scalb(norm, -k);*/
	twok = 1.0 * pow(FLT_RADIX, k);
	norm = norm * pow(FLT_RADIX, -k);

	/* Compute exp(x/2^k)-I

	 exp(x/2^k) - I - x/2^k = r
	 (x/2^k)^f / f! = t
	 */
	for (m = 0; m < n2; m++) {
		t[m] = x2k[m] /= twok;
		r[m] = 0;
	}
	norm = 1;
	for (f = 2; norm > f * 0.01 * DBL_EPSILON; f++) {
		norm = 0;
		for (in = 0; in < n2; in += n) {
			x2kin = x2k + in;
			for (j = 0; j < n; j++) {
				tj = t + j;
				cij = 0;
				for (mn = m = 0; m < n; m++, mn += n)
					cij += x2kin[m] * tj[mn];
				p = c[in + j] = cij / f;
				if (p < 0)
					p = -p;
				if (p > norm)
					norm = p;
			}
		}
		for (m = 0; m < n2; m++)
			r[m] += (t[m] = c[m]);
	}
	/* add x/2^k */
	for (m = 0; m < n2; m++)
		r[m] += x2k[m];

	/* Square k times */
	while (k-- > 0) {
		mmul(r, r, c, n);
		for (m = 0; m < n2; m++)
			r[m] = 2 * r[m] + c[m];
		/* if the diagonal becomes significant, then add the
		 identity and continue normal powering */
		for (i = 0; i < n; i++)
			if ( ABS(r[i*n+i]) > 0.5)
				break;
		if (i < n)
			break;
	}

	/* Add identity */
	for (i = 0; i < n; i++)
		r[i * n + i] += 1;

	/* Square k times (whatever is left) */
	while (k-- > 0) {
		mmul(r, r, c, n);
		for (m = 0; m < n2; m++)
			r[m] = c[m];
	}
}

/* a deletion in both sequences (like in an MSA, forced by some other
 insertion) has to be completely eliminated			*/
double ComputeScore(int idms, Counters *cnt, DayMatrix *DMS) {
	int i;
	double cj, *Simp;
	/*if (ISNULL(DMS[idms]->DelCost))*/
	cj = cnt->delF * DMS[idms].FixedDel + cnt->delI * DMS[idms].IncDel;
	/*else
	 for (cj = i = 0; cnt->delL[i] > 0; i++)
	 cj += ((cnt->delL[i] <= DELCOSTARRAYSIZE) ?
	 D(DMS[idms]->DelCostArray, cnt->delL[i]) :
	 DBL(
	 eval(
	 New4(FUNCTION, DMS[idms]->DelCost,
	 Newint(cnt->delL[i]),
	 NewNumber(DMS[idms]->PamNumber)))));*/
	Simp = DMS[idms].Simi;
	for (i = cnt->ninds - 1; i >= 0; i--)
		cj += Simp[cnt->inds[i]];
	return (cj);
}

/* On a failure of finding the maximum, consider enlarging MAXPOINTS */
int FindMaxPamR(int x[], double f[], int n, double *res, Counters *cnt,
		DayMatrix *DMS);

static void FindMaxPam(char *o1, char *o2, int len, DayMatrix *DMS, int DMSLen,
		Counters *cnt, int *dmsnumb, double *sim) {
	double f[3 * MAXPOINTS + 1];
	int i, j, ind, k, x[3 * MAXPOINTS + 1], l;

	/*if (len >= MAXSEQLEN + 20000)
	 userror("molecular sequence too long");*/
	cnt->delF = cnt->delI = cnt->ninds = 0;
	for (ind = i = j = l = 0; i < len; i++) {
		if (o1[i] == '_') {
			if (o2[i] == '_')
				continue; /* do not count del==del */
			if (ind == 2) {
				cnt->delL[j++] = l + 1;
				l = 0;
			}
			if (ind == 1) {
				cnt->delI++;
				l++;
			} else {
				cnt->delF++;
				ind = 1;
			}
		} else if (o2[i] == '_') {
			if (ind == 1) {
				cnt->delL[j++] = l + 1;
				l = 0;
			}
			if (ind == 2) {
				cnt->delI++;
				l++;
			} else {
				cnt->delF++;
				ind = 2;
			}
		} else {
			if (/*DMS[1]->Mapping != 0 ||*/(o1[i] != 'X' && o2[i] != 'X'))
				cnt->inds[cnt->ninds++] = /*MapSymbol(o1[i], DMS[1])*/o1[i]
						* MATRIX_DIM + /*MapSymbol(o2[i], DMS[1])*/o2[i];
			if (ind != 0) {
				cnt->delL[j++] = l + 1;
				ind = l = 0;
			}
		}
		/*if (j >= MAXSEQLEN / 10)
		 userror("too many deletions for EstimatePam");*/
	}
	if (ind != 0)
		cnt->delL[j++] = l + 1;
	cnt->delL[j] = 0;

	for (k = 0; k <= 3 * MAXPOINTS; k++) {
		x[k] = k * (DMSLen - 2) / (3 * MAXPOINTS) + 1;
		f[k] = ComputeScore(x[k], cnt, DMS);
	}
	*dmsnumb = FindMaxPamR(x, f, 3 * MAXPOINTS + 1, sim, cnt, DMS);
}

/* Multi-point multiple maximum searching a la Brent */
int FindMaxPamR(int x[], double f[], int n, double *res, Counters *cnt,
		DayMatrix *DMS) {
	int j, j1, j2, x2[3 * MAXPOINTS + 1];
	double cj, c1, c2, f2[3 * MAXPOINTS + 1];

	restart:
	/*DBG( if( n<2 ) notimpl("FindMaxPamR");
	 for (j = 1; j < n; j++)
	 if (x[j] < x[j - 1])
	 notimpl("FindMaxPamR2");
	 )*/

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
		cj = ComputeScore(j, cnt, DMS);
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
			j2 = FindMaxPamR(x2, f2, n - j, &c2, cnt, DMS);
			j1 = FindMaxPamR(x, f, j + 1, &c1, cnt, DMS);
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

/* Input: cnt, dmsnumb, DMS.
 Output: Sim, PamNumber, PamVariance, MLPamDistance (global) */
static void ComputeVarianceML(Counters *cnt, int dmsnumb, DayMatrix *DMS,
		int DMSLen, double *Sim, double *PamNumber, double *PamVariance,
		double* logPAM1) {
	double /*logPAM1[MAXMUTDIM * MAXMUTDIM],*/Mp[MAXMUTDIM * MAXMUTDIM],
			logMp[MAXMUTDIM * MAXMUTDIM], log2Mp[MAXMUTDIM * MAXMUTDIM],
			t1[MAXMUTDIM * MAXMUTDIM], t2[MAXMUTDIM * MAXMUTDIM], t3[MAXMUTDIM
					* MAXMUTDIM], t4[MAXMUTDIM * MAXMUTDIM];
	double /*cj,*/FD1d, logL1, logL2, p, incr, x1;
	int f21, i, iter, j, n, n2;
	/*ALGEB s;

	 s = DMS[1]->logPAM1;*/
	n = /*LENGTH(s[1]) - 1*/26;
	/*if (!((n == 20 && DMS[1]->dim == 21) || n == DMS[1]->dim))
	 userror("logPAM1 and DMS dimensions do not agree");
	 if (LENGTH(DMS) < 6)
	 userror("not enough Dayhoff matrices to compute PAM-variance");*/

	/* is this always 1 in darwin?? */
	f21 = /*n == 20 && DMS[1]->dim == 21*/0;
	/*ListToArray(s, logPAM1, n, n);*/
	/* is this always 400 in darwin? */
	n2 = n * n;

	/*  model FixedDel = FD0 + FD1*log10(pam) + ... == Score
	 Score = 10*log10(Prob) = 10/ln(10)*ln(Prob)

	 FixedDel[1] = FD0+FD1*log[10](pam[1]);
	 FixedDel[n] = FD0+FD1*log[10](pam[n]);
	 sol := solve({%%,%},{FD0,FD1});
	 diff( ln(10)/10 * (FD0 + FD1*log[10](pam)), pam );
	 normal( subs(sol,%) );
	 */
	for (j = DMSLen; j > 3 && DMS[j - 1].FixedDel >= -1; j--)
		;
	{
		double FD1/*, v1, v2*/;
		FD1 = log(10.0) * (DMS[1].FixedDel - DMS[j - 1].FixedDel)
				/ log(DMS[1].PamNumber / DMS[j - 1].PamNumber);
		/*v1 = FD1 * (log10(DMS[j / 2].PamNumber) - log10(DMS[1].PamNumber));*/
		/*v2 = DMS[j / 2].FixedDel - DMS[1].FixedDel;*/
		/*if (fabs(v1 - v2) > 1e-8 * fabs(v2))*/
		/* be a bit more tolerant due to rounding of Day entries */
		/*userror("deletion cost is not logarithmic on pam distance");*/
		FD1d = cnt->delF * FD1 / 10;
	}
	p = DMS[dmsnumb].PamNumber;
	/* printf( "DMS[%d]->PamNumber=%.12g, FD1d=%.12g\n", dmsnumb, p, FD1d ); */

	for (i = 0; i < n2; i++)
		t1[i] = p * logPAM1[i];
	mexp(t1, t2, t3, Mp, n);
	mmul(logPAM1, Mp, logMp, n);
	mmul(logPAM1, logMp, log2Mp, n);

	for (iter = 0; 1; iter++) {
		/*double tot1, tot2, tot3;*/
		/*DBG( for( i=0; i < n;
		 i++
		 ) {
		 tot1 = tot2 = tot3 = 0;
		 for( j=0; j<n; j++ ) {
		 tot1 += Mp[j*n+i];
		 tot2 += logMp[j*n+i];
		 tot3 += log2Mp[j*n+i];
		 }
		 if( fabs(tot1-1) > 1e-11 || fabs(tot2) > 1e-12 ||
		 fabs(tot3) > 1e-12 ) {
		 printf( "tot1=%g, tot2=%g, tot3=%g\n", tot1-1, tot2, tot3 );
		 userror("matrices computed incorrectly, failed assertion");
		 }
		 }
		 )*/
		logL1 = FD1d / p;
		logL2 = -FD1d / (p * p);
		for (i = cnt->ninds - 1; i >= 0; i--) {
			j = cnt->inds[i];
			if (f21)
				j -= j / 21;
			/* if the matrix is exactly diagonal, ignore entry (possibly an X) */
			if (Mp[j] == 0 && logMp[j] == 0 && log2Mp[j] == 0)
				continue;
			x1 = logMp[j] / Mp[j];
			logL1 += x1;
			logL2 += log2Mp[j] / Mp[j] - x1 * x1;
		}

		if (logL2 >= 0)
			break;
		incr = -logL1 / logL2;
		/* printf( "p:=%.10g: incr:=%.10g: logL1:=%.10g: logL2:=%.10g:\n",
		 p, incr, logL1, logL2 ); */
		/* for cases where the maximum is at +infinity, make sure
		 that we do not go beyond the limit of the DMS array.
		 Cases like this may happen when there are too few matches
		 and deletions, and the deletion gain is less than the
		 loss due to the matches, and then pam -> infinity */
		if (p + incr > DMS[DMSLen - 1].PamNumber) {
			incr = DMS[DMSLen - 1].PamNumber - p;
		}
		x1 = 0.04 * (p + 25);
		if (fabs(incr) > x1)
			incr = (incr < 0) ? -x1 : x1;
		p += incr;
		if (fabs(incr) < 1e-8 * p)
			break;
		/*if (iter >= 10)
		 userror(
		 "failure to compute PAM by ML (likely to be a problem with DMS)");*/
		for (i = 0; i < n2; i++)
			t1[i] = incr * logPAM1[i];
		mexp(t1, t2, t3, t4, n);
		mmul(Mp, t4, t1, n);
		memcpy(Mp, t1, n2 * sizeof(double));
		mmul(logMp, t4, t1, n);
		memcpy(logMp, t1, n2 * sizeof(double));
		mmul(log2Mp, t4, t1, n);
		memcpy(log2Mp, t1, n2 * sizeof(double));
	}
	if (logL2 >= 0) {
		double rmpam;
		if (p >= DMS[DMSLen - 1].PamNumber) {
			p = DMS[DMSLen - 1].PamNumber;
			/* if at the end of the DMS, estimate the variance by the
			 idealized mutations model, it is the best we can do */
			rmpam = pow(1 - 0.20 / 19, -p);
			logL2 = -cnt->ninds * 19 * log(1 - 0.20 / 19) * log(1 - 0.20 / 19)
					/ ((rmpam - 1) * (rmpam + 19));
		} /*else
		 userror("ComputeVariance by ML found a minimum not a maximum");*/
	}
	/*if ((dmsnumb > 1 && p < pamDistances[dmsnumb - 1])
	 || (dmsnumb < DMSLen - 1 && p > pamDistances[dmsnumb + 1])) {*/
	/* printf( "p=%.12g, dmsnumb=%d\n", p, dmsnumb );
	 if( dmsnumb > 1 ) printf( "DMS[dmsnumb-1]->PamNumber=%.12g, LENGTH(DMS)=%ld\n", DMS[dmsnumb-1]->PamNumber, LENGTH(DMS) );
	 if( dmsnumb < LENGTH(DMS)-1 ) printf( "DMS[dmsnumb+1]->PamNumber=%.12g\n", DMS[dmsnumb+1]->PamNumber ); */
	/*userror(
	 "ComputeVariance by ML, internal conflict, probably bad DMS, or DayMatrix cannot be derived from rate matrix (logPAM1), or stop codons may be present in the sequences");
	 }*/
	/*assign(naminstall("MLPamDistance"), NewNumber(p), 0);*/
	*PamVariance = -1 / logL2;

	*Sim = ComputeScore(dmsnumb, cnt, DMS);
	*PamNumber = DMS[dmsnumb].PamNumber;

	/* check the neighbours to see if they are better */
	/*for (j = (dmsnumb > 3 ? dmsnumb - 2 : 1); j < DMSLen && j <= dmsnumb + 2;
	 j++)
	 if (j != dmsnumb) {
	 cj = ComputeScore(j, cnt, DMS);
	 if (cj > *Sim)
	 userror("maximum likelihood Distance not correct");
	 }*/
}

static void ComputeSmallVariance(Counters *cnt, int dmsnumb, DayMatrix *DMS,
		double *Sim, double *PamNumber, double *PamVariance, double* logPAM1) {
	double a1, a2, b0, t1, t2, t3, t4, t5, ave, var;
	double /*L1[MAXMUTDIM * MAXMUTDIM],*/L2[MAXMUTDIM * MAXMUTDIM], L3[MAXMUTDIM
			* MAXMUTDIM];
	int i, ij, n;
	double* L1 = logPAM1;
	/*ALGEB p;*/

	/* Create L1=logPAM1, L2=L1^2 and L3=L1^3.  */
	/*p = DMS[1]->logPAM1;*/
	n = 26;
	/*if (!((n == 20 && DMS[1]->dim == 21)
	 || n == DMS[1]->dim))
	 userror(
	 "logPAM1 and DMS dimensions do not agree");*/
	/*ListToArray(p, L1, n, n);*/
	mmul(L1, L1, L2, n);
	mmul(L1, L2, L3, n);

	/*
	 >         convert( series( Dij, p, 3 ), polynom );
	 /                   2\
                       L1ij        L2ij p   |    L3ij       L2ij |  2
	 ln(p) + ln(----) + 1/2 ------ + |1/6 ---- - 1/8 -----| p
	 f[i]         L1ij    |    L1ij           2|
	 \               L1ij /

	 >         convert( series( Dii, p, 3 ), polynom );
	 1                                   2   2
	 ln(----) + L1ii p + (1/2 L2ii - 1/2 L1ii ) p
	 f[i]
	 */
	a1 = a2 = 0;
	b0 = cnt->delF * 0.7434;
	for (i = cnt->ninds - 1; i >= 0; i--) {
		ij = cnt->inds[i];
		/* ij is an index for the DayMatrix array, correct if 20/21 */
		/*if (n == 20 && DMS[1]->dim == 21)
		 ij -= ij / 21;*/
		if (L1[ij] == 0)
			continue;
		if (ij % (n + 1) == 0) {
			a1 += L1[ij];
			a2 += 0.5 * (L2[ij] - L1[ij] * L1[ij]);
		} else {
			b0 += 1;
			a1 += 0.5 * L2[ij] / L1[ij];
			a2 += (L3[ij] / 6 - L2[ij] * L2[ij] / 8 / L1[ij]) / L1[ij];
		}
	}

	t5 = b0 * b0;
	t3 = a1 * a1;
	t2 = 2.0 * t3 + a2 * t5;
	t4 = b0 + 1.0;
	t1 = (3.0 * b0 + 2.0) * a2 + t2;
	ave = -t4 * ((5.0 * b0 + 6.0) * a2 + t2) / a1 / t1;
	*PamVariance = var = (4.0 * t3 * t3
			+ ((32.0 + 4.0 * t5 + 24.0 * b0) * t3
					+ (28.0 * b0 + 12.0 + (23.0 + 8.0 * b0 + t5) * t5) * a2)
					* a2) / t3 * t4 / (t1 * t1);
	if (a1 == 0 || t1 == 0) {
		ave = DMS[dmsnumb].PamNumber;
		*PamVariance = ave * ave;
	}
	/*assign(naminstall("ExpectedPamDistance"), NewNumber(ave), 0);*/
	/*assign(naminstall("MLPamDistance"), Newint(0), 0);*/

	*PamNumber = DMS[dmsnumb].PamNumber;
	*Sim = ComputeScore(dmsnumb, cnt, DMS);
}

void EstimatePam(char* o1, char* o2, int len, DayMatrix* DMS, int DMSLen,
		double* logPAM1, double* result)

{
	Counters cnt;
	int dmsnumb;
	double Sim, PamNumber, PamVariance;

	/* Increase the length to be compatible with the core */
	DMSLen += 1;

	/*if (length(t[1]) != length(t[2]))
	 userror("s1 and s2 must have the same length");*/
	/* return an arbitrary distance/variance for an empty string */
	/*if (length(t[1]) == 0)
	 return (New2(LIST,
	 New4(EXPSEQ, Newint(0), Newint(250),
	 Newint(62500))));*/
	/*DMS = (DayMatrix**) (t[3][1]);*/
	FindMaxPam(o1, o2, len, DMS, DMSLen, &cnt, &dmsnumb, &Sim);
	if (dmsnumb <= 1) {
		ComputeSmallVariance(&cnt, dmsnumb, DMS, &Sim, &PamNumber, &PamVariance,
				logPAM1);
	} else {
		ComputeVarianceML(&cnt, dmsnumb, DMS, DMSLen, &Sim, &PamNumber,
				&PamVariance, logPAM1);
	}

	result[0] = Sim;
	result[1] = PamNumber;
	result[2] = PamVariance;

}

DayMatrix* createDayMatrices(double* gapOpen, double* gapExt,
		double* pamDistances, long long* matrix_pointers, int DMSLen) {
	int i, j;
#ifdef PY_DEBUG
	int k, l;
#endif
	/* indexing starts from one */
	DayMatrix* ret = (DayMatrix*) malloc((DMSLen + 1) * sizeof(DayMatrix));

	for (i = 1; i <= DMSLen; ++i) {
#ifdef PY_DEBUG
		printf("Creating day matrix:\n");
		printf("FixedDel = %f\nIncDel = %f\nPamNumber = %f\n", gapOpen[i], gapExt[i], pamDistances[i]);
		printf("Simi at %lx:\n", (long)matrix_pointers[i]);
		for(k = 0; k < MATRIX_DIM; ++k) {
			printf("\t[");
			for(l = 0; l < MATRIX_DIM; ++l) {
				printf("%.3f, ", ((double*)(matrix_pointers[i]))[k * MATRIX_DIM + l]);
			}
			printf("]\n");
		}
#endif
		ret[i].FixedDel = gapOpen[i];
		ret[i].IncDel = gapExt[i];
		ret[i].PamNumber = pamDistances[i];
		ret[i].Simi = (double*) malloc(
		MATRIX_DIM * MATRIX_DIM * sizeof(double));
		for (j = 0; j < MATRIX_DIM * MATRIX_DIM; ++j) {
			ret[i].Simi[j] = ((double*) (matrix_pointers[i]))[j];
		}
	}

	return ret;
}

void freeDayMatrices(DayMatrix* DMS, int DMSLen) {
	int i;
	for (i = 1; i <= DMSLen; ++i) {
		free(DMS[i].Simi);
	}
	free(DMS);
}

void CreateOrigDayMatrix(double* log_pam1, double PamNum, double* new_matrix)

/*  Computations following the paper:

 A Model of Evolutionary Change in Proteins by
 M.O. Dayhoff, R.M. Schwartz and B.C. Orcutt.

 and our own limiting computation	Gaston H. Gonnet (Sep 1990) */

{
	/*double AA[MAXMUTDIM][MAXMUTDIM], norm, RelMut[MAXMUTDIM], RowsA[MAXMUTDIM];*/
	double /*a, b, lambda,*/M[MAXMUTDIM * MAXMUTDIM], MP[MAXMUTDIM * MAXMUTDIM],
			t1[MAXMUTDIM * MAXMUTDIM], t2[MAXMUTDIM * MAXMUTDIM];
	double /*k, */pam1, pam2, tot, /*v, */maxs, t;
	/*ALGEB p, res, logpam;*/
	int i, j, n;
	/*double *dm;*/
	int d1, d2;

	d1 = d2 = 26;

	/*if (ID(Mutations) != LIST)
	 userror2(Mutations, "invalid argument");*/
	n = /*LENGTH(Mutations[1]) - 1;*/d1;
	/*if (numeric(PamNum))*/
	pam1 = pam2 = PamNum;
	/*else {
	 pam1 = DBL(A(PamNum[1]));
	 pam2 = DBL(A(PamNum[2]));
	 }*/
	/*if (pam1 <= 0 || pam2 > 10000 || pam1 > pam2 || (pam1 < pam2 && pam1 != 1))
	 userror("incorrect Pam number");
	 if (n < 1 || n > MAXMUTDIM)
	 userror2(Mutations, "dimension out of range");*/

	/*if (LastMutMatrix == Mutations && LastAAFreq == AaCounts
	 && LastLogPAM1 == valuenv_ININ(ININlogPAM1)) {*/
	/* Dayhoff was called earlier with the same Mutations/AaCounts */
	/*ListToArray(LastLogPAM1, M, n, n);
	 logpam = LastLogPAM1;
	 goto resume;
	 }*/

	/*if (AaCounts == 0 *//* signal that a logPAM1 matrix is passed *//*) {*/
	/*ListToArray(Mutations, M, n, n);
	 LastAAFreq = 0;
	 logpam = Mutations;
	 goto resume;
	 }

	 for (i = 0; i < n; i++) {
	 if (ID(Mutations) != LIST)
	 userror("CreateOrigDayMatrix: invalid mutations matrix");
	 p = A(Mutations[1][i + 1]);
	 if (ID(p) == LISTNUM)
	 for (j = 0; j <= i; j++) {
	 AA[i][j] = AA[j][i] = D(p, j + 1);
	 if (AA[i][j] < 0)
	 userror("negative entries in Mutations matrix");
	 }
	 else if (ID(p) == LIST)
	 for (j = 0; j <= i; j++) {
	 AA[i][j] = AA[j][i] = DBL(A(p[1][j + 1]));
	 if (AA[i][j] < 0)
	 userror("negative entries in Mutations matrix");
	 }
	 else
	 userror("CreateOrigDayMatrix: invalid mutations matrix");
	 }
	 norm = 0;

	 for (i = 0; i < n; i++) {
	 RowsA[i] = 0.0;
	 for (j = 0; j < n; j++)
	 RowsA[i] += AA[i][j];
	 if (ID(AaCounts) == LIST)
	 NorFre[i] = DBL(A(AaCounts[1][i + 1]));
	 else
	 NorFre[i] = D(AaCounts, i + 1);
	 if (NorFre[i] <= 0)
	 userror("non-positive frequecies");
	 norm += NorFre[i];
	 }

	 for (i = 0; i < n; i++) {
	 NorFre[i] /= norm;
	 RelMut[i] = RowsA[i] / NorFre[i];
	 }*/

	/* PAM 1 matrix; solve the equation for lambda */
	/*a = b = 0.0;
	 for (i = 0; i < n; i++) {
	 a += NorFre[i];
	 b -= NorFre[i] * RelMut[i];
	 }*/
	/* 100*(a + lambda*b) = 99 */
	/*lambda = (0.99 - a) / b;

	 for (i = 0; i < n; i++) {
	 for (j = 0; j < n; j++)
	 M[i * n + j] = lambda * RelMut[j] * AA[i][j] / RowsA[j];
	 M[i * (n + 1)] = -lambda * RelMut[i];
	 }
	 assign(ININlogPAM1, NewMatrix(M, n, n), 0);
	 LastMutMatrix = Mutations;
	 LastAAFreq = AaCounts;
	 logpam = LastLogPAM1 = valuenv_ININ(ININlogPAM1);

	 resume:*/
	/* the following variables must be set: pam1, M, pam2, NorFre */
	/* if AaCounts==0, then NorFre has to be set after computing
	 the exp(M) */
	for (i = 0; i < d1 * d2; ++i) {
		M[i] = log_pam1[i];
	}

	if (pam1 != 1) {
		for (i = 0; i < n * n; i++) {
			M[i] *= pam1;
		}
	}

	mexp(M, t1, t2, MP, n);
	/*if (AaCounts == 0 *//* signal that a logPAM1 matrix is passed *//*) {*/
	for (tot = i = 0; i < n; i++) {
		if (MP[i] == 0 || MP[i * n] == 0)
			NorFre[i] = 0;
		else
			tot += NorFre[i] = MP[i * n] / MP[i];
	}
	for (i = 0; i < n; i++)
		NorFre[i] /= tot;
	/*}*/
	/*if (pam1 < pam2) {
	 res = Newl1(I(pam2 - pam1 + 2), EXPSEQ);
	 for (i = 0; i < n * n; i++)
	 M[i] = MP[i];
	 } else*/
	/*res = 0;*/

	/* To insure that the algorithms for dynamic programming and
	 backwards dynamic programming work without the need of
	 checking for approximate values, we will work with doubles
	 which are <integer>/2^v.

	 The largest possible positive score (negative scores do
	 not matter too much) should be representable as an exact
	 sum of all of its components.  I.e

	 MAXSEQLEN * MaxSim * 2^v < 2^DBL_MANT_DIG (usually 53)
	 or
	 2^v < 2 / (DBL_EPSILON * MAXSEQLEN * MaxSim) */

	/*for (k = pam1; k <= pam2; k += 1) {*/
	/* Compute Dayhoff Matrix */
	maxs = 0.0;
	/* if there is a negative or 0 entry in the mutation matrix,
	 this means that there were not enough observations for
	 this mutation.  We assign it an arbitrarily low value,
	 representing a probability of 1e-5 */
	for (i = 0; i < d1; i++) {
		for (j = i; j < d1; j++) {
			t = new_matrix[d2 * i + j] = new_matrix[d2 * j + i] =
					NorFre[i] <= 0 || MP[i * d1 + j] <= 0 ?
							/*-50.0*/0 : 10 * log10(MP[i * d1 + j] / NorFre[i]);
			if (t > maxs) {
				maxs = t;
			}
		}
	}

	/*v = ilogb(1 / (DBL_EPSILON * MAXSEQLEN * (maxs < 1 ? 1 : maxs)));*/

	/* TODO fix this RoundDM macro? */
	/*for (i = 0; i < d1; i++) {
	 for (j = 0; j < d1; j++) {
	 new_matrix[d2 * i + j] = RoundDM(new_matrix[d2 * i + j], v);
	 }
	 }*/
	/*if (n == 20)
	 dm->type = naminstall("Peptide");
	 else if (n == 64)
	 dm->type = naminstall("Nucleotide");
	 else
	 dm->type = ININunknown;

	 eval(A(dm));
	 if (pam1 == pam2)*/
	/*return (A(dm));*/
	/*res[I(k - pam1 + 1)] = A2(dm);

	 if (k < pam2) {*/
	/* MP = MP*M */
	/*mmul(MP, M, t1, n);
	 for (i = 0; i < n * n; i++)
	 MP[i] = t1[i];
	 }
	 }

	 return (New2(LIST, res));*/
}
