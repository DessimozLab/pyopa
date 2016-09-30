/*
 * DynProgr_sse_double.c
 *
 *  Created on: Sep 11, 2014
 *      Author: machine
 */

/* According to http://stackoverflow.com/questions/1919183/how-to-allocate-and-free-aligned-memory-in-c */

#include "Python_extension.h"
#include "DynProgr_sse_double.h"

void *aligned_malloc(int align, int size) {
	char *mem = malloc(size + align + sizeof(char*));
	char **ptr = (char**) ((long) (mem + align + sizeof(char*)) & ~(align - 1));
	ptr[-1] = mem;
	return ptr;
}

void aligned_free(void *ptr) {
	free(((void**) ptr)[-1]);
}

void free_profile_double_sse(ProfileDouble* profile) {
	aligned_free(profile->profile);
	aligned_free(profile->new_opt);
	aligned_free(profile->old_opt);
	aligned_free(profile->new_rd);
	free(profile);
}

ProfileDouble* createProfileDoubleSSE(const char* s1, int ls1, double* matrix) {
	ProfileDouble* profileDouble = (ProfileDouble*) malloc(
			sizeof(ProfileDouble));

	int i, j;

	/* the segment length */
	int segLength = (ls1 + 1) / 2;

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

	/* the new row deletion score */
	double * new_rd = (double*) aligned_malloc(16,
			(2 * segLength + 2) * sizeof(double));

	__m128d m128_profile = _mm_load_pd(profile); /* the score profile */

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

	profileDouble->ls1 = ls1;
	profileDouble->profile = profile;
	profileDouble->new_opt = new_opt;
	profileDouble->old_opt = old_opt;
	profileDouble->new_rd = new_rd;

	return profileDouble;
}

double align_double_local(ProfileDouble* profileDouble, const char *s2, int ls2,
		double gap_open, double gap_ext, double threshold, int* max1, int* max2) {
	int i, j, k;
	/*int AToInts2[MAXSEQLEN + 1];*/

	int segLength = (profileDouble->ls1 + 1) / 2;
	double DelFixed = gap_open, DelIncr = gap_ext/*, *Score_s1*/;
	double /*Tcd, t, */MaxScore/*, Sj, Sj1, Tj, Tj1, Trd*/;

	/* This totcells was a system variable and I have no idea what it is used for */
	double totcells = 0;
	double* swap;
	double * current_profile = profileDouble->profile;
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
	__m128d m128_profile = _mm_load_pd(profileDouble->profile); /* the score profile */
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
		/*AToInts2[j] = MapSymbol(s2[j-1],DM)s2[j - 1];*/
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

	/*if (profileDouble->ls1 > 2) {*/
	/* initialize the other arrays used for the dynProg code */
	for (i = 0; i < 2 * segLength; i += 2) {
		_mm_store_pd(profileDouble->old_opt + i, _mm_setzero_pd());
		_mm_store_pd(profileDouble->new_opt + i, _mm_setzero_pd());
		_mm_store_pd(profileDouble->new_rd + i, _mm_setzero_pd());
	}

	for (j = 0; j < ls2; j++) {
		/* set the column deletion score to zero, has to be fixed later on */
		m128_new_cd = _mm_setzero_pd();

		/* shift old_opt */
		/* compute the temporal max score depending on the previous column */
		m128_new_opt = _mm_set_pd(profileDouble->new_opt[2 * segLength - 2], 0);

		/* compute the current profile, depending on the character in s2 */
		/*if (s2[j] == '_')
		 userror("underscores cannot be used in sequence alignment");*/
		current_profile = profileDouble->profile
				+ /*MapSymbol(s2[j], DM)*/s2[j] * segLength * 2;

		/* swap the old optimal score with the new one */
		swap = profileDouble->new_opt;
		profileDouble->new_opt = profileDouble->old_opt;
		profileDouble->old_opt = swap;

		for (i = 0; i < 2 * segLength; i += 2) {
			m128_profile = _mm_load_pd(current_profile + i);
			m128_new_opt = _mm_add_pd(m128_new_opt, m128_profile);

			m128_MaxScore = _mm_max_pd(m128_MaxScore, m128_new_opt);

			m128_new_rd = _mm_load_pd(profileDouble->new_rd + i);
			m128_new_opt = _mm_max_pd(m128_new_opt, m128_new_rd);

			m128_new_opt = _mm_max_pd(m128_new_opt, m128_new_cd);

			m128_new_opt = _mm_max_pd(m128_new_opt, m128_zero);

			_mm_store_pd(profileDouble->new_opt + i, m128_new_opt);

			m128_new_opt = _mm_add_pd(m128_new_opt, m128_DelFixed);

			m128_new_rd = _mm_add_pd(m128_new_rd, m128_DelIncr);

			m128_new_rd = _mm_max_pd(m128_new_rd, m128_new_opt);
			_mm_store_pd(profileDouble->new_rd + i, m128_new_rd);

			m128_new_cd = _mm_add_pd(m128_new_cd, m128_DelIncr);

			m128_new_cd = _mm_max_pd(m128_new_cd, m128_new_opt);

			m128_new_opt = _mm_load_pd(profileDouble->old_opt + i);
		}

		/* set totcells */
		totcells += profileDouble->ls1;

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
				if (profileDouble->new_opt[i + k - 1] > MaxScore_old) {
					MaxScore_old = profileDouble->new_opt[i + k - 1];
					*max1 = (k - 1) * segLength + i / 2 + 1; /*+1 since the other versions use 1<=i<=ls1 and I use 0<=i<ls1*/
					*max2 = j + 1;
				}
			}
		}

		/* if the goal was reached, exit */
		if (MaxScore > threshold) {
			/*WSReturn();*/
			return ((double) MaxScore);
		}

		/* cleaning up with missed arrows (this part of the code seems to be
		 having (almost) no effect on the runtime,
		 therefore I didnot bother vectorizing it) */

		_mm_store_pd(temp1, m128_new_cd);
		new_cd1 = temp1[0];
		gapOpen = profileDouble->new_opt[segLength * 2 - 2] + DelFixed;

		/* filling S and coldel for backtracking */
		coldel[j + 1] = temp1[1];

		for (i = 0; i < segLength * 2 && new_cd1 >= gapOpen; i += 2) {
			if (profileDouble->new_opt[i + 1] < new_cd1)
				profileDouble->new_opt[i + 1] = new_cd1;
			new_cd1 += DelIncr;
			gapOpen = profileDouble->new_opt[i + 1] + DelFixed;
			if (profileDouble->new_rd[i + 1] < gapOpen)
				profileDouble->new_rd[i + 1] = gapOpen;
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
		S[j + 1] = profileDouble->new_opt[profileDouble->ls1 - 1
				- profileDouble->ls1 % 2];
	}
	/*} else {*/

	/*for (i = 1; i <= profileDouble->ls1; i++) {
	 Tj1 = Sj1 = 0;
	 Trd = MINUSINF;*/

	/* setup Score_s1 */
	/*if (s1[i - 1] == '_')
	 userror("underscores cannot be used in sequence alignment");*/
	/*k = s1[i - 1];
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
	 totcells += j;
	 return (MaxScore);
	 }
	 }

	 coldel[j] = Tcd;
	 S[j] = Tj1 = Tj;
	 Sj1 = Sj;
	 }
	 totcells += ls2;
	 }*/
	/*printf("NotImplemented!");
	 }*/

	/*WSReturn();*/
#ifdef PY_DEBUG
	printf("c_align_double: return %f\n", (double) MaxScore);
#endif
	return ((double) MaxScore);

}
