/*
 * Python_extension.c
 *
 *  Created on: Jul 8, 2014
 *      Author: machine
 */

#include <stdlib.h>
#include <stdio.h>
#include <float.h>

#include "Python_extension.h"
#include "DynProgr_scalar.h"
#include "DynProgr_sse_byte.h"
#include "DynProgr_sse_short.h"
#include "DynProgr_sse_double.h"

double coldel[MAXSEQLEN + 1], S[MAXSEQLEN + 1];

char* denormalize(const char* str, int len) {
	char* ret = (char*) malloc((len + 1) * sizeof(char));
	int i;

	for (i = 0; i < len; ++i) {
		ret[i] = 'A' + str[i];
	}

	ret[len] = '\0';

	return ret;
}

void debug_alignment(char* name, void* profile, const char* s2, int ls2,
		Options* options) {
	char* dns2 = denormalize(s2, ls2);
	printf(
			"Aligning %s: profile = %lu, s2 = %s, len(s2) = %d, gapOpen = %f, gapExt = %f, threshold = %f\n",
			name, (unsigned long) profile, dns2, ls2, options->gapOpen,
			options->gapExt, options->threshold);
	free(dns2);
}

void debug_profile(char* name, void* pb, const char* q, int queryLen,
		void* matrix, int isShort) {
	char* query = denormalize(q, queryLen);
	int i, j;

	printf(
			"Successfully created %s profile at %lu from query = %s, len(query) = %d\n",
			name, (unsigned long) pb, query, queryLen);
	printf("The matrix used for the %s profile: \n", name);

	for (i = 0; i < 26; ++i) {
		printf("[");
		for (j = 0; j < 26; ++j) {
			double val;
			if (isShort) {
				val = (double) ((SMatrix) matrix)[i * 26 + j];
			} else {
				val = (double) ((BMatrix) matrix)[i * 26 + j];
			}

			printf("%+.2f, ", val);
		}
		printf("]\n");
	}

	free(query);
}

double c_align_profile_byte_sse_local(ProfileByte* profile, const char *s2,
		int ls2, double gapOpen, double gapExt, double threshold) {
	Options options;
	options.gapOpen = gapOpen;
	options.gapExt = gapExt;
	options.threshold = threshold;

#ifdef PY_DEBUG
	debug_alignment("BYTE", profile, s2, ls2, &options);
#endif

	return swps3_alignmentByteSSE(profile, s2, ls2, &options);
}
double c_align_profile_short_sse_local(ProfileShort* profile, const char *s2,
		int ls2, double gapOpen, double gapExt, double threshold) {
	Options options;
	options.gapOpen = gapOpen;
	options.gapExt = gapExt;
	options.threshold = threshold;

#ifdef PY_DEBUG
	debug_alignment("SHORT", profile, s2, ls2, &options);
#endif

	return swps3_alignmentShortSSE(profile, s2, ls2, &options);
}

ProfileByte* c_create_profile_byte_sse_local(const char* query, int queryLen,
		BMatrix matrix) {
	ProfileByte* pb = swps3_createProfileByteSSE(query, queryLen, matrix);
#ifdef PY_DEBUG
	debug_profile("BYTE", pb, query, queryLen, matrix, 0);
#endif
	return pb;
}
ProfileShort* c_create_profile_short_sse_local(const char* query, int queryLen,
		SMatrix matrix) {
	ProfileShort* ps = swps3_createProfileShortSSE(query, queryLen, matrix);
#ifdef PY_DEBUG
	debug_profile("SHORT", ps, query, queryLen, matrix, 1);
#endif
	return ps;
}

void c_free_profile_byte_sse_local(ProfileByte* profile) {
#ifdef PY_DEBUG
	printf("Freeing BYTE profile at %lu\n", (unsigned long)profile);
#endif
	swps3_freeProfileByteSSE(profile);
}
void c_free_profile_short_sse_local(ProfileShort* profile) {
#ifdef PY_DEBUG
	printf("Freeing SHORT profile at %lu\n", (unsigned long)profile);
#endif
	swps3_freeProfileShortSSE(profile);
}

double c_align_scalar_reference_local(DMatrix matrix, const char *s1, int ls1,
		const char *s2, int ls2, double gapOpen, double gapExt,
		double threshold) {
	Options options;
	options.gapOpen = gapOpen;
	options.gapExt = gapExt;
	options.threshold = threshold;

	return swps3_alignScalar(matrix, s1, ls1, s2, ls2, &options);
}

int c_align_strings(double* matrix, char *s1, int len1, char *s2, int len2,
		double escore, char *o1, char *o2, double maxerr, double gap_open,
		double gap_ext) {
	int i, i1, j;
	int *DelFrom1 = NULL;
	double *S1 = NULL, *coldel1 = NULL;
	double Slen2i, maxs, t;
	char *tprs1 = NULL, *tprs2 = NULL, *prs1, *prs2;
	double tot;

	/* make len1 >= len2 */
	if (len1 < len2) {
		prs1 = s1;
		s1 = s2;
		s2 = prs1;
		i = len1;
		len1 = len2;
		len2 = i;
		prs1 = o1;
		o1 = o2;
		o2 = prs1;
	}
	/*ASSERT( len2 >= 0, Newint(len2) );*/

	/* s2 is null */
	if (len2 == 0) {
		/*if( escore != DeletionCost(len1, gap_open, gap_ext) )
		 userror("score incompatible with aligned strings");*/
		for (i = 0; i < len1; i++) {
			o1[i] = s1[i];
			o2[i] = '_';
		}
		return (len1);
	}

	/* 1 against 1 (this case is needed for recursion termination) */
	if (len1 == 1) {
		if (2 * gap_open >= escore - maxerr) {
			o1[0] = s1[0];
			o2[0] = '_';
			o1[1] = '_';
			o2[1] = s2[0];
			return (2);
		}
		/*if( (NoSelf && s1==s2) || escore > DMScore(s1[0],s2[0],matrix) )
		 userror("Alignment/Match has incorrect data");*/
		o1[0] = s1[0];
		o2[0] = s2[0];
		return (1);
	}

	/* equal length, try to see if there is a match without deletions */
	if (len1 == len2 && s1 != s2) {
		for (i = tot = 0; i < len1; i++)
			tot += matrix[s1[i] * MATRIX_DIM + s2[i]];
		if (tot >= escore) {
			for (i = 0; i < len1; i++) {
				o1[i] = s1[i];
				o2[i] = s2[i];
			}
			return (len1);
		}
	}

	/* len1 >= len2 >= 1, try to see if there is an all deleted match */
	/* allow some tolerance for a Match read with score error */
	if (2 * gap_open + (len1 + len2 - 2) * gap_ext >= escore - maxerr) {
		for (i = 0; i < len1; i++) {
			o1[i] = s1[i];
			o2[i] = '_';
		}
		for (i = 0; i < len2; i++) {
			o1[len1 + i] = '_';
			o2[len1 + i] = s2[i];
		}
		return (len1 + len2);
	}

	/* reverse strings */
	if (s1 - s2 >= 0 && s2 + len2 - s1 > 0) {
		j = MMAX(s1 - s2 + len1, len2);
		tprs1 = (char*) malloc(j * sizeof(char));
		for (i = 0; i < j; i++)
			tprs1[j - 1 - i] = s2[i];
		prs1 = tprs1 + (j + s2 - s1 - len1);
		prs2 = tprs1 + (j - len2);
	} else if (s2 - s1 >= 0 && s1 + len1 - s2 > 0) {
		j = MMAX(s2 - s1 + len2, len1);
		tprs1 = (char*) malloc(j * sizeof(char));
		for (i = 0; i < j; i++)
			tprs1[j - 1 - i] = s1[i];
		prs2 = tprs1 + (j + s1 - s2 - len2);
		prs1 = tprs1 + (j - len1);
	} else {
	    tprs1 = (char*) malloc(len1 * sizeof(char));
	    tprs2 = (char*) malloc(len2 * sizeof(char));
		for (i = 0; i < len2; i++)
			tprs2[len2 - 1 - i] = s2[i];
		for (i = 0; i < len1; i++)
			tprs1[len1 - 1 - i] = s1[i];
		prs1 = tprs1;
		prs2 = tprs2;
	}

	/* divides s1 in half */
	i1 = len1 / 2;
	c_align_double_global(matrix, s1, i1, s2, len2, gap_open, gap_ext);
	S1 = (double*) malloc((len2+1) * sizeof(double));
	coldel1 = (double*) malloc((len2 + 1) * sizeof(double));
	DelFrom1 = (int*) malloc((len2 + 1) * sizeof(int));
	for (i = 0; i <= len2; i++) {
		S1[i] = S[i]; coldel1[i] = coldel[i]; DelFrom1[i] = DelFrom[i];
	}
	c_align_double_global(matrix, prs1, len1 - i1, prs2, len2, gap_open,
			gap_ext);

	/* Find the best sum of scores from the two halves */
	maxs = -DBL_MAX;
	for (j = 0; j <= len2; j++) {
		t = S1[j] + S[len2 - j];
		if (t > maxs) {
			maxs = t;
			i = j;
		}
		t = coldel1[j] + coldel[len2 - j] - gap_open + gap_ext;
		if (t > maxs) {
			maxs = t;
			i = j;
		}
	}

	/* splitting on a match */
	if (maxs == S1[i] + S[len2 - i]) {
		Slen2i = S[len2 - i];
		j = c_align_strings(matrix, s1, i1, s2, i, S1[i], o1, o2, 0.0, gap_open,
				gap_ext);
		j += c_align_strings(matrix, s1 + i1, len1 - i1, s2 + i, len2 - i,
				Slen2i, o1 + j, o2 + j, 0.0, gap_open, gap_ext);
		free(tprs1); free(tprs2); free(S1); free(DelFrom1); free(coldel1);
		return (j);
	}

	/* splitting on a vertical deletion */
	{
		int i3, i4, len;
		i3 = DelFrom1[i] - 1;
		i4 = len1 - DelFrom[len2 - i] + 2;
		Slen2i = coldel[len2 - i];
		len = c_align_strings(matrix, s1, i3, s2, i,
				coldel1[i] - gap_open - gap_ext * (i1 - i3 - 1), o1, o2, 0.0,
				gap_open, gap_ext);
		for (j = i3 + 1; j < i4; j++) {
			o1[len] = s1[j - 1];
			o2[len++] = '_';
		}
		len += c_align_strings(matrix, s1 + i4 - 1, len1 - i4 + 1, s2 + i,
				len2 - i, Slen2i - gap_open - gap_ext * (i4 - i1 - 2), o1 + len,
				o2 + len, 0.0, gap_open, gap_ext);
		free(tprs1); free(tprs2); free(S1); free(DelFrom1); free(coldel1);
		return (len);
	}
	return 0;
}

double c_align_double_global(double* matrix, const char *s1, int ls1,
		const char *s2, int ls2, double gap_open, double gap_ext) {
	int i, j, k;
	int *AToInts2 = NULL;
	double DelFixed, DelIncr, *Score_s1 = NULL;
	double t, t2/*, MaxScore*/, rowdel, Sj1;
	/*double vScore[MAXMUTDIM];*/
	int NoSelf = 0;
	/* This totcells was a system variable and I have no idea what it is used for */
	double totcells = 0.0;

	DelFixed = gap_open;
	DelIncr = gap_ext;

	/*MaxScore = MINUSINF;*/
	AToInts2 = (int*) malloc((ls2 + 1) * sizeof(int));
	S[0] = coldel[0] = 0;
	for (j = 1; j <= ls2; j++) {
		/*if (s2[j - 1] == '_')
		 userror("underscores cannot be used in sequence alignment");*/
		AToInts2[j] = /*MapSymbol(s2[j - 1], DM)*/s2[j - 1];
		coldel[j] = MINUSINF;
		S[j] = /*(mode == CFE || mode == Local) ? 0 : */S[j - 1]
				+ (j == 1 ? DelFixed : DelIncr);
	}

	DelFrom[0] = 1;
	for (i = 1; i <= ls1; i++) {

		Sj1 = S[0];
		coldel[0] += i == 1 ? DelFixed : DelIncr;
		S[0] = /*(mode == CFE || mode == Local) ? 0 :*/coldel[0];
		rowdel = MINUSINF;

		/*if (ProbSeqCnt == 0) {*/
		/* setup Score_s1 */
		/*if (s1[i - 1] == '_')
		 userror("underscores cannot be used in sequence alignment");*/
		k = /*MapSymbol(s1[i - 1], DM)*/s1[i - 1];
		Score_s1 = matrix + k * MATRIX_DIM;
		/*} else if (ProbSeqCnt == 1) {
		 ProbScore(s1.ps + (i - 1) * NrDims, M, AF, NrDims, vScore);
		 Score_s1 = vScore;
		 }*/

		/* Complete code for the inner loop of dynamic programming.
		 Any changes should be made to this code and generate the
		 others by explicitly eliminating the tests tested on the
		 outside.  The purpose is to maximize speed of the most
		 common case.		March 27, 2005 */
		for (j = 1; j <= ls2; j++) {
			/* current row is i (1..ls1), column is j (1..ls2) */

			coldel[j] += DelIncr;
			t = (t2 = S[j]) + DelFixed;
			if (coldel[j] < t) {
				coldel[j] = t;
				DelFrom[j] = i;
			}

			rowdel += DelIncr;
			t = S[j - 1] + DelFixed;
			if (rowdel < t)
				rowdel = t;

			/* TODO: check that det vs prob can never be a Self */
			t = NoSelf /*&& ProbSeqCnt == 0*/&& s1 + i == s2 + j ?
					MINUSINF : Sj1 + Score_s1[AToInts2[j]];
			if (t < rowdel)
				t = rowdel;
			if (t < coldel[j])
				t = coldel[j];

			/*if (mode == Local) {
			 if (t < 0)
			 t = 0;
			 if (t > MaxScore) {
			 MaxScore = t;
			 *Max1 = i;
			 *Max2 = j;
			 if (MaxScore >= goal) {
			 totcells += j;
			 return (MaxScore);
			 }
			 }
			 } else if (mode == Shake && t > MaxScore) {
			 MaxScore = t;
			 *Max1 = i;
			 *Max2 = j;
			 if (MaxScore >= goal) {
			 totcells += j;
			 return (MaxScore);
			 }
			 }*/
			S[j] = t;
			Sj1 = t2;
		}
		totcells += ls2;

		/*if ((mode == CFE || mode == CFEright) && S[ls2] > MaxScore) {
		 MaxScore = S[ls2];
		 *Max1 = i;
		 *Max2 = ls2;
		 if (MaxScore >= goal)
		 return (MaxScore);
		 }*/

	}

	/*if (mode == CFE || mode == CFEright)
	 for (j = 1; j <= ls2; j++)
	 if (S[j] > MaxScore) {
	 MaxScore = S[j];
	 *Max1 = ls1;
	 *Max2 = j;
	 }

	 if (mode == Global) {
	 *Max1 = ls1;
	 *Max2 = ls2;*/
	free(AToInts2);
	return (S[ls2]);
	/*}
	 return (MaxScore);*/
}
