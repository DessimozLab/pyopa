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

/* According to http://stackoverflow.com/questions/1919183/how-to-allocate-and-free-aligned-memory-in-c */

void *aligned_malloc(int align, int size) {
	char *mem = malloc(size + align + sizeof(char*));
	char **ptr = (char**) ((long) (mem + align + sizeof(char*)) & ~(align - 1));
	ptr[-1] = mem;
	return ptr;
}

void aligned_free(void *ptr) {
	free(((void**) ptr)[-1]);
}

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
	int DelFrom1[MAXSEQLEN + 1], i, i1, j/*, M1*/;
	double S1[MAXSEQLEN + 1], coldel1[MAXSEQLEN + 1];
	double Slen2i, maxs, t;
	char rs1[2 * MAXSEQLEN], rs2[MAXSEQLEN], *prs1, *prs2;
	double tot;
	/*int Global = 2;*/
	/* this NoSelf could be a parameter */
	/*int NoSelf = 0;*/

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
			tot += /*DMScore(s1[i],s2[i],matrix)*/matrix[s1[i] * MATRIX_DIM
					+ s2[i]];
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
		for (i = 0; i < j; i++)
			rs1[j - 1 - i] = s2[i];
		prs1 = rs1 + (j + s2 - s1 - len1);
		prs2 = rs1 + (j - len2);
	} else if (s2 - s1 >= 0 && s1 + len1 - s2 > 0) {
		j = MMAX(s2 - s1 + len2, len1);
		for (i = 0; i < j; i++)
			rs1[j - 1 - i] = s1[i];
		prs2 = rs1 + (j + s1 - s2 - len2);
		prs1 = rs1 + (j - len1);
	} else {
		for (i = 0; i < len2; i++)
			rs2[len2 - 1 - i] = s2[i];
		for (i = 0; i < len1; i++)
			rs1[len1 - 1 - i] = s1[i];
		prs1 = rs1;
		prs2 = rs2;
	}

	/* divides s1 in half */
	i1 = len1 / 2;
	/*seq1.ds = s1;
	 seq2.ds = s2;*/
	c_align_double_global(matrix, s1, i1, s2, len2, gap_open, gap_ext);
	for (i = 0; i <= len2; i++) {
		S1[i] = S[i];
		coldel1[i] = coldel[i];
		DelFrom1[i] = DelFrom[i];
	}
	/*seq1.ds = prs1;
	 seq2.ds = prs2;*/
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
	/*if( !(maxerr != 0 || maxs==escore) ) {
	 printf( "DM->bits=%d, gap_open=%.18g, gap_ext=%.18g\n",
	 DM->bits, gap_open, gap_ext );
	 printf( "maxerr=%.18g, maxs=%.18g, escore=%.18g\n",
	 maxerr, maxs, escore );
	 }
	 ASSERT2( maxerr != 0 || maxs==escore, NewNumber(maxs), NewNumber(escore) );
	 if( ENVprintlevel > 0 && maxs < escore-maxerr )
	 fprintf( w_unit, "Warning: DynProgStrings could not reach "
	 "score %.12g, reached %.12g instead\n", escore, maxs );*/

	/* splitting on a match */
	if (maxs == S1[i] + S[len2 - i]) {
		Slen2i = S[len2 - i];
		j = c_align_strings(matrix, s1, i1, s2, i, S1[i], o1, o2, 0.0, gap_open,
				gap_ext);
		j += c_align_strings(matrix, s1 + i1, len1 - i1, s2 + i, len2 - i,
				Slen2i, o1 + j, o2 + j, 0.0, gap_open, gap_ext);
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
		return (len);
	}
	return 0;
}

double c_align_double_local(double* matrix, const char *s1, int ls1,
		const char *s2, int ls2, double gap_open, double gap_ext,
		double threshold, int* max1, int* max2) {
	unsigned int i, j, k;
	int AToInts2[MAXSEQLEN + 1];

	double DelFixed = gap_open, DelIncr = gap_ext, *Score_s1;
	double Tcd, t, MaxScore, Sj, Sj1, Tj, Tj1, Trd;

	/* the segment length */
	unsigned int segLength = (ls1 + 1) / 2;

	/* This totcells was a system variable and I have no idea what it is used for */
	double totcells;

	/* allocate the arrays needed locally and align them to 16 by hand */
	/*the score profile*/
	double * profile = (double*) aligned_malloc(16,
			(2 * segLength * MATRIX_DIM + 2) * sizeof(double));

	/* the old optimal score */
	double * old_opt = (double*) aligned_malloc(16,
			(2 * segLength + 2) * sizeof(double));

	/* the new optimal score */
	double * new_opt = (double*) aligned_malloc(16,
			(2 * segLength + 2) * sizeof(double));

	double * swap;

	/* the new row deletion score */
	double * new_rd = (double*) aligned_malloc(16,
			(2 * segLength + 2) * sizeof(double));

	double * current_profile = profile;
	double new_cd1 = 0;
	double new_cd = 0;
	double MaxScore_old = -1;

	double gapOpen;

	/* m128_DelIncr = [DelIncr,DelIncr] */
	__m128d m128_DelIncr = _mm_load1_pd(&DelIncr);
	/* m128_DelFixed = [DelFixed,DelFixed] */
	__m128d m128_DelFixed = _mm_load1_pd(&DelFixed);

	/* m128_zero = [0,0] */
	__m128d m128_zero = _mm_setzero_pd();
	__m128d m128_temp1;
	double * temp1 = (double*) &m128_temp1; /* m128_MaxScore = [0,0] */
	__m128d m128_MaxScore = _mm_setzero_pd();

	double * tempMax = (double*) &m128_MaxScore;
	__m128d m128_profile = _mm_load_pd(profile); /* the score profile */
	__m128d m128_new_opt; /* the new optimal score */
	__m128d m128_new_rd; /* the new row deletion score */
	__m128d m128_new_cd = _mm_setzero_pd();

#ifdef PY_DEBUG
	char* ns1 = denormalize(s1, ls1);
	char* ns2 = denormalize(s2, ls2);

	printf("c_align_double: s1 = %s, len(s1) = %d, s2 = %s, len(s2) = %d, gap_open = %f, gap_ext = %f, threshold = %f\n",
			ns1, ls1, ns2, ls2, gap_open, gap_ext, threshold);

	free(ns1);
	free(ns2);
#endif

	MaxScore = 0;
	S[0] = coldel[0] = 0;
	for (j = 1; j <= ls2; j++) {
		/*if( s2[j-1]=='_' )
		 userror("underscores cannot be used in sequence alignment");*/
		AToInts2[j] = /*MapSymbol(s2[j-1],DM)*/s2[j - 1];
		coldel[j] = MINUSINF;
		S[j] = 0;
	}

	/* this version of the code does two lines per loop iteration */

	/* This version of the code implements the idea presented in
	 *
	 ************************************************************
	 * Six-fold speed-up of Smith-Waterman sequence database searches using
	 * parallel processing on common microprocessors
	 * Torbj��rn Rognes 1, and Erling Seeberg 1
	 * Bioinformatics Vol. 16 no. 8 2000
	 ************************************************************
	 *
	 * This version implements all the features present in NewDynProgrCode2
	 * and is using doubles
	 * e.g. early exit if the 'goal' was reached, updating totcells, computing
	 * the location where the maximum is achieved.
	 */

	/* initializing the profile */
	if (ls1 > 2) {
		for (i = 0; i < MATRIX_DIM; i++) {
			for (j = 0; j < 2 * segLength; j += 2) {
				/*if (s1[j / 2] == '_')
				 userror("underscores cannot be used in sequence alignment");*/

				if (j / 2 + segLength < ls1) {
					/*if (s1[j / 2 + segLength] == '_')
					 userror(
					 "underscores cannot be used in sequence alignment");*/

					m128_profile = _mm_set_pd((matrix + (/*MapSymbol(s1[j / 2 + segLength],
					 DM)*/s1[j / 2 + segLength]) * MATRIX_DIM)[i],
							(matrix
									+ (/*MapSymbol(s1[j / 2], DM)*/s1[j / 2])
											* MATRIX_DIM)[i]);
					_mm_store_pd(profile + i * 2 * segLength + j, m128_profile);
				} else {
					m128_profile = _mm_set_pd(0,
							(matrix
									+ /*MapSymbol(s1[j / 2], DM)*/s1[j / 2]
											* MATRIX_DIM)[i]);
					_mm_store_pd(profile + i * 2 * segLength + j, m128_profile);
				}
			}
		}

		/* initialize the other arrays used for the dynProg code */
		for (i = 0; i < 2 * segLength; i += 2) {
			_mm_store_pd(old_opt + i, _mm_setzero_pd());
			_mm_store_pd(new_opt + i, _mm_setzero_pd());
			_mm_store_pd(new_rd + i, _mm_setzero_pd());
		}

		for (j = 0; j < ls2; j++) {
			/* set the column deletion score to zero, has to be fixed later on */
			m128_new_cd = _mm_setzero_pd();

			/* shift old_opt */
			/* compute the temporal max score depending on the previous column */
			m128_new_opt = _mm_set_pd(new_opt[2 * segLength - 2], 0);

			/* compute the current profile, depending on the character in s2 */
			/*if (s2[j] == '_')
			 userror("underscores cannot be used in sequence alignment");*/
			current_profile = profile
					+ /*MapSymbol(s2[j], DM)*/s2[j] * segLength * 2;

			/* swap the old optimal score with the new one */
			swap = new_opt;
			new_opt = old_opt;
			old_opt = swap;

			for (i = 0; i < 2 * segLength; i += 2) {
				m128_profile = _mm_load_pd(current_profile + i);
				m128_new_opt = _mm_add_pd(m128_new_opt, m128_profile);

				m128_MaxScore = _mm_max_pd(m128_MaxScore, m128_new_opt);

				m128_new_rd = _mm_load_pd(new_rd + i);
				m128_new_opt = _mm_max_pd(m128_new_opt, m128_new_rd);

				m128_new_opt = _mm_max_pd(m128_new_opt, m128_new_cd);

				m128_new_opt = _mm_max_pd(m128_new_opt, m128_zero);

				_mm_store_pd(new_opt + i, m128_new_opt);

				m128_new_opt = _mm_add_pd(m128_new_opt, m128_DelFixed);

				m128_new_rd = _mm_add_pd(m128_new_rd, m128_DelIncr);

				m128_new_rd = _mm_max_pd(m128_new_rd, m128_new_opt);
				_mm_store_pd(new_rd + i, m128_new_rd);

				m128_new_cd = _mm_add_pd(m128_new_cd, m128_DelIncr);

				m128_new_cd = _mm_max_pd(m128_new_cd, m128_new_opt);

				m128_new_opt = _mm_load_pd(old_opt + i);
			}

			/* set totcells */
			totcells += ls1;

			/* check for a changed MaxScore, if so, find the location */
			k = 0;
			MaxScore_old = MaxScore;
			if (tempMax[0] > MaxScore) {
				MaxScore = tempMax[0];
				k = 1;
			}
			if (tempMax[1] > MaxScore) {
				MaxScore = tempMax[1];
				k = 2;
			}

			/* if the MaxScore changed, go and search for it in the segment */
			if (k >= 1) {
				for (i = 0; i < segLength * 2; i += 2) {
					if (new_opt[i + k - 1] > MaxScore_old) {
						MaxScore_old = new_opt[i + k - 1];
						*max1 = (k - 1) * segLength + i / 2 + 1; /*+1 since the other versions use 1<=i<=ls1 and I use 0<=i<ls1*/
						*max2 = j + 1;
					}
				}
			}

			/* if the goal was reached, exit */
			if (MaxScore > threshold) {
				/*WSReturn();*/
				aligned_free(profile);
				aligned_free(old_opt);
				aligned_free(new_opt);
				return ((double) MaxScore);
			}

			/* cleaning up with missed arrows (this part of the code seems to be
			 having (almost) no effect on the runtime,
			 therefore I didnot bother vectorizing it) */

			_mm_store_pd(temp1, m128_new_cd);
			new_cd1 = temp1[0];
			gapOpen = new_opt[segLength * 2 - 2] + DelFixed;

			/* filling S and coldel for backtracking */
			coldel[j + 1] = temp1[1];

			for (i = 0; i < segLength * 2 && new_cd1 >= gapOpen; i += 2) {
				if (new_opt[i + 1] < new_cd1)
					new_opt[i + 1] = new_cd1;
				new_cd1 += DelIncr;
				gapOpen = new_opt[i + 1] + DelFixed;
				if (new_rd[i + 1] < gapOpen)
					new_rd[i + 1] = gapOpen;
				if (i == 2 * segLength - 2) {
					if (ls2 % 2 == 1) {
						/* depending on whether cd changed at the last position we */
						/* have to adapt coldel[j+1] */
						coldel[j + 1] = new_cd1 - DelIncr;
					} else if (new_cd > gapOpen) {
						coldel[j + 1] = new_cd1;
					}
				}
			}
			/* filling S and coldel for backtracking */
			S[j + 1] = new_opt[ls1 - 1 - ls1 % 2];
		}
	} else {

		for (i = 1; i <= ls1; i++) {
			Tj1 = Sj1 = 0;
			Trd = MINUSINF;

			/* setup Score_s1 */
			/*if (s1[i - 1] == '_')
			 userror("underscores cannot be used in sequence alignment");*/
			k = /*MapSymbol(s1[i - 1], DM)*/s1[i - 1];
			Score_s1 = matrix + k * MATRIX_DIM;

			for (j = 1; j <= ls2; j++) {
				Sj = S[j];
				Tcd = coldel[j] + DelIncr;

				if (Tcd < (t = Sj + DelFixed))
					Tcd = t;
				Tj = Sj1 + Score_s1[AToInts2[j]];

				Trd += DelIncr;
				if (Trd < (t = Tj1 + DelFixed))
					Trd = t;

				if (Tj < Tcd)
					Tj = Tcd;
				if (Tj < Trd)
					Tj = Trd;
				if (Tj < 0)
					Tj = 0;
				if (Tj > MaxScore) {
					MaxScore = Tj;
					*max1 = i;
					*max2 = j;
					if (MaxScore >= threshold) {
						/*WSReturn();*/
						totcells += j;
						aligned_free(profile);
						aligned_free(old_opt);
						aligned_free(new_opt);
						return (MaxScore);
					}
				}

				coldel[j] = Tcd;
				S[j] = Tj1 = Tj;
				Sj1 = Sj;
			}
			totcells += ls2;
		}
	}

	/*WSReturn();*/
	aligned_free(profile);
	aligned_free(old_opt);
	aligned_free(new_opt);
#ifdef PY_DEBUG
	printf("c_align_double: return %f\n", (double) MaxScore);
#endif
	return ((double) MaxScore);

}

double c_align_double_global(double* matrix, const char *s1, int ls1,
		const char *s2, int ls2, double gap_open, double gap_ext) {
	unsigned int i, j, k;
	int AToInts2[MAXSEQLEN + 1];
	double DelFixed, DelIncr, *Score_s1 = NULL;
	double t, t2, MaxScore, rowdel, Sj1;
	double vScore[MAXMUTDIM];
	int NoSelf = 0;
	/* This totcells was a system variable and I have no idea what it is used for */
	double totcells;

	DelFixed = gap_open;
	DelIncr = gap_ext;

	MaxScore = MINUSINF;
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
	return (S[ls2]);
	/*}
	 return (MaxScore);*/
}
