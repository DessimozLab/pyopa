/** \file DynProgr_sse_byte.c
 *
 * Profile generation and alignment for packed byte vectors on SSE2.
 */
/*
 * Copyright (c) 2007-2008 ETH Zürich, Institute of Computational Science
 *
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following
 * conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

#include "DynProgr_sse_byte.h"
#include "Page_size.h"
#include "debug.h"
#include <stdio.h>
#include <float.h>

#define PAGE_ALIGN(x) (((size_t)(x)+getPageSize()-1)&~(getPageSize()-1))
/**
 *  Creates a profile with unsigned 8 bit integers
 */

/*
void print128_num(__m128i * var) {
    uint32_t *val = (uint32_t*) &var;
    printf("Numerical: %d %d %d %d \n", 
           val[0], val[1], val[2], val[3]);
}
*/

//TODO: remove query (empty string) - not used anywhere in the swps3_createProfileByteSSE except for debug
EXPORT ProfileByte * swps3_createProfileByteSSE(const char * query,
		int queryLen, BMatrix matrix) {
	int segLen = (queryLen + 15) / 16;
//	printf("Query Length: %d\n", queryLen);
	int i, j, k;
	int bias = 0;
	uint8_t * pprofile;
	ProfileByte * profile = malloc(
			sizeof(ProfileByte) + segLen * (MATRIX_DIM + 3) * sizeof(__m128i )
					+ 64 + 2 * getPageSize());

	profile->loadOpt = (__m128i *) ((size_t)(profile->data + 15) & ~(0xf));
//	printf("%s", "Loadopt:\n");
//	print128_num(profile->loadOpt);

	profile->storeOpt = profile->loadOpt + segLen;
//      	printf("%s", "Storeopt:\n");
//        print128_num(profile->storeOpt);


	profile->rD = profile->storeOpt + segLen;
	profile->profile = (__m128i *) PAGE_ALIGN(profile->rD + segLen);

	/* Init the profile */
	profile->len = queryLen;
	/* Init the byte profile */
	for (i = 0; i < MATRIX_DIM; i++)
		for (j = 0; j < MATRIX_DIM; j++)
			if (bias < -matrix[i * MATRIX_DIM + j])
				bias = -matrix[i * MATRIX_DIM + j];
	pprofile = (uint8_t*) profile->profile;

	for (i = 0; i < MATRIX_DIM; i++)
		for (j = 0; j < segLen; j++)
			for (k = 0; k < 16; k++)
				if (j + k * segLen < queryLen)
					*(pprofile++) = matrix[query[j + k * segLen] * MATRIX_DIM
							+ i] + bias;
				else
					*(pprofile++) = bias;
	profile->bias = bias;

#ifdef DEBUG
	for(i=0; i<queryLen; ++i) debug("\t%c",query[i]+'A');
	debug("\n");
#endif
	return profile;
}

EXPORT double swps3_alignmentByteSSE_lin(ProfileByte * query, const char * db,
		int dbLen, Options * options) {

	/**********************************************************************
	 * This version of the code implements the idea presented in
	 *
	 ***********************************************************************
	 * Striped Smith-Waterman speeds database searches six times over other
	 * SIMD implementations
	 *
	 * Michael Farrar, Bioinformatics, 23(2), pp. 156-161, 2007
	 **********************************************************************/

	int i, j;
	unsigned char MaxScore = 0;
	int segLength = (query->len + 15) / 16; /* the segment length */

	__m128i * loadOpt = query->loadOpt;
	__m128i * storeOpt = query->storeOpt;
	__m128i * current_profile;
	__m128i * swap;

	__m128i vMinimums = _mm_set1_epi32(0);

	__m128i vDelFixed = _mm_set1_epi8(-options->gapOpen);
	__m128i vBias = _mm_set1_epi8(query->bias);

	__m128i vMaxScore = vMinimums; /* vMaxScore = [0,0] */

	__m128i vStoreOpt; /* the new optimal score */
	__m128i vRD; /* the new row deletion score */
	__m128i vCD = vMinimums; /* the column deletion score */
	__m128i zero = vMinimums; /* the column deletion score */
	__m128i vTmp;
#ifdef DEBUG
	int ii,jj;
#endif

	/* initialize the other arrays used for the dynProg code */
	/*********************************************************/
	for (i = 0; LIKELY(i < segLength); i++) {
		_mm_store_si128(loadOpt + i, zero);
		_mm_store_si128(storeOpt + i, zero);
	}

	/* looping through all the columns */
	/***********************************/

	for (j = 0; LIKELY(j < dbLen); j++) {

		/* compute the opt and cd score depending on the previous column
		 *******************************************************************
		 * set the column deletion score to zero, has to be fixed later on */
		vCD = zero;

		/* set the opt score to the elements computed in the previous column*/
		/* set the low of storeOpt to MaxS[j]                               */
		vStoreOpt = _mm_load_si128(storeOpt + segLength - 1);
		vStoreOpt = _mm_slli_si128(vStoreOpt, 1);

		/* compute the current profile, depending on the character in s2 */
		/*****************************************************************/
		current_profile = query->profile + db[j] * segLength;

		/* swap the old optimal score with the new one */
		/***********************************************/
		swap = storeOpt;
		storeOpt = loadOpt;
		loadOpt = swap;

		/* main loop computing the max, precomputing etc. */
		/**************************************************/
		for (i = 0; LIKELY(i < segLength); i++) {
			vTmp = _mm_load_si128(loadOpt + i);
			vRD = _mm_subs_epu8(vTmp, vDelFixed);

			/* add the profile the prev. opt */
			vStoreOpt = _mm_adds_epu8(vStoreOpt, *(current_profile + i));
			vStoreOpt = _mm_subs_epu8(vStoreOpt, vBias);

			/* update the maxscore found so far (gaps only decrease score) */
			vMaxScore = _mm_max_epu8(vMaxScore, vStoreOpt);

			/* compute the correct opt score of the cell */
			vStoreOpt = _mm_max_epu8(vStoreOpt, vRD);
			vStoreOpt = _mm_max_epu8(vStoreOpt, vCD);

			/* store the opt score of the cell */
			_mm_store_si128(storeOpt + i, vStoreOpt);

			/* precompute cd for next iteration */
			vCD = _mm_subs_epu8(vStoreOpt, vDelFixed);

			/* load precomputed opt for next iteration */
			vStoreOpt = vTmp;
		}

		for (i = 0; LIKELY(i < 16); ++i) {
			int k;
			/* compute the gap extend penalty for the current cell */
			vCD = _mm_slli_si128(vCD, 1);

			for (k = 0; LIKELY(k < segLength); ++k) {
				/* compute the current optimal value of the cell */
				vTmp = _mm_load_si128(storeOpt + k);
				vStoreOpt = _mm_max_epu8(vTmp, vCD);
				_mm_store_si128(storeOpt + k, vStoreOpt);

				/* break if vStoreOpt unchanged */
				if (UNLIKELY(
						_mm_movemask_epi8(_mm_cmpeq_epi8(vTmp, vStoreOpt))
								== 0xFFFF))
					goto shortcut;

				/* precompute the scores for the next cell */
				vCD = _mm_subs_epu8(vStoreOpt, vDelFixed);
			}
		}
		shortcut:

#ifdef DEBUG
		debug("%c\t",db[j]);
		for(ii=0; ii<16;++ii) {
			for(jj=0; jj<segLength;++jj) {
				if(ii*segLength+jj < query->len)
				debug("%d\t",(int)((unsigned char*)storeOpt)[ii+jj*16]);
			}
		}
		debug("\n");
#else
		i = 0;
#endif
	}

	vMaxScore = _mm_max_epu8(vMaxScore, _mm_srli_si128(vMaxScore, 8));
	vMaxScore = _mm_max_epu8(vMaxScore, _mm_srli_si128(vMaxScore, 4));
	vMaxScore = _mm_max_epu8(vMaxScore, _mm_srli_si128(vMaxScore, 2));
	vMaxScore = _mm_max_epu8(vMaxScore, _mm_srli_si128(vMaxScore, 1));
	MaxScore = (unsigned char) _mm_extract_epi16(vMaxScore, 0);
	if ((int) MaxScore + (int) query->bias >= 255)
		return DBL_MAX;
	return ((double) MaxScore);
}

EXPORT double swps3_alignmentByteSSE(ProfileByte * query, const char * db,
		int dbLen, Options * options) {

	/**********************************************************************
	 * This version of the code implements the idea presented in
	 *
	 ***********************************************************************
	 * Striped Smith-Waterman speeds database searches six times over other
	 * SIMD implementations
	 *
	 * Michael Farrar, Bioinformatics, 23(2), pp. 156-161, 2007
	 **********************************************************************/

	int i, j;
	unsigned char MaxScore = 0;
	int segLength = (query->len + 15) / 16; /* the segment length */

	__m128i * loadOpt = query->loadOpt;
	__m128i * storeOpt = query->storeOpt;
	__m128i * rD = query->rD;
	__m128i * current_profile;
	__m128i * swap;

	__m128i vMinimums = _mm_set1_epi32(0);

	__m128i vDelIncr = _mm_set1_epi8(-options->gapExt);
	__m128i vDelFixed = _mm_set1_epi8(-options->gapOpen);
	__m128i vBias = _mm_set1_epi8(query->bias);

	__m128i vMaxScore = vMinimums; /* vMaxScore = [0,0] */

	__m128i vStoreOpt; /* the new optimal score */
	__m128i vRD; /* the new row deletion score */
	__m128i vCD = vMinimums; /* the column deletion score */
	__m128i zero = vMinimums; /* the column deletion score */
	__m128i vTmp;
#ifdef DEBUG
	int ii,jj;
#endif

	if (options->gapExt <= options->gapOpen) {
		return swps3_alignmentByteSSE_lin(query, db, dbLen, options);
	}

	/* initialize the other arrays used for the dynProg code */
	/*********************************************************/
	for (i = 0; LIKELY(i < segLength); i++) {
		_mm_store_si128(loadOpt + i, zero);
		_mm_store_si128(storeOpt + i, zero);
		_mm_store_si128(rD + i, zero);
	}

	/* looping through all the columns */
	/***********************************/

	for (j = 0; LIKELY(j < dbLen); j++) {
		/* compute the opt and cd score depending on the previous column
		 *******************************************************************
		 * set the column deletion score to zero, has to be fixed later on */
		vCD = zero;

		/* set the opt score to the elements computed in the previous column*/
		vStoreOpt = _mm_load_si128(storeOpt + segLength - 1);
		vStoreOpt = _mm_slli_si128(vStoreOpt, 1);

		/* compute the current profile, depending on the character in s2 */
		/*****************************************************************/
		current_profile = query->profile + db[j] * segLength;

		/* swap the old optimal score with the new one */
		/***********************************************/
		swap = storeOpt;
		storeOpt = loadOpt;
		loadOpt = swap;

		/* main loop computing the max, precomputing etc. */
		/**************************************************/
		for (i = 0; LIKELY(i < segLength); i++) {
			vRD = _mm_load_si128(rD + i);
			vRD = _mm_subs_epu8(vRD, vDelIncr);
			vTmp = _mm_load_si128(loadOpt + i);
			vTmp = _mm_subs_epu8(vTmp, vDelFixed);
			vRD = _mm_max_epu8(vRD, vTmp);
			_mm_store_si128(rD + i, vRD);

			/* add the profile the prev. opt */
			vStoreOpt = _mm_adds_epu8(vStoreOpt, *(current_profile + i));
			vStoreOpt = _mm_subs_epu8(vStoreOpt, vBias);

			/* update the maxscore found so far (gaps only decrease score) */
			vMaxScore = _mm_max_epu8(vMaxScore, vStoreOpt);

			/* compute the correct opt score of the cell */
			vStoreOpt = _mm_max_epu8(vStoreOpt, vRD);
			vStoreOpt = _mm_max_epu8(vStoreOpt, vCD);

			/* store the opt score of the cell */
			_mm_store_si128(storeOpt + i, vStoreOpt);

			/* precompute cd for next iteration */
			vStoreOpt = _mm_subs_epu8(vStoreOpt, vDelFixed);
			vCD = _mm_subs_epu8(vCD, vDelIncr);
			vCD = _mm_max_epu8(vCD, vStoreOpt);

			/* load precomputed opt for next iteration */
			vStoreOpt = _mm_load_si128(loadOpt + i);
		}

		for (i = 0; LIKELY(i < 16); ++i) {
			int k;
			/* compute the gap extend penalty for the current cell */
			vCD = _mm_slli_si128(vCD, 1);

			for (k = 0; LIKELY(k < segLength); ++k) {
				/* compute the current optimal value of the cell */
				vStoreOpt = _mm_load_si128(storeOpt + k);
				vStoreOpt = _mm_max_epu8(vStoreOpt, vCD);

				_mm_store_si128(storeOpt + k, vStoreOpt);

				/* precompute the scores for the next cell */
				vStoreOpt = _mm_subs_epu8(vStoreOpt, vDelFixed);
				vCD = _mm_subs_epu8(vCD, vDelIncr);

				if (UNLIKELY(
						_mm_movemask_epi8(
								_mm_cmpeq_epi8(_mm_subs_epu8(vCD, vStoreOpt),
										zero)) == 0xFFFF))
					goto shortcut;
			}
		}
		shortcut:

#ifdef DEBUG
		debug("%c\t",db[j]);
		for(ii=0; ii<16;++ii) {
			for(jj=0; jj<segLength;++jj) {
				if(ii*segLength+jj < query->len)
				debug("%d\t",(int)((unsigned char*)storeOpt)[ii+jj*16]);
			}
		}
		debug("\n");
#else
		i = 0;
#endif
	}

	vMaxScore = _mm_max_epu8(vMaxScore, _mm_srli_si128(vMaxScore, 8));
	vMaxScore = _mm_max_epu8(vMaxScore, _mm_srli_si128(vMaxScore, 4));
	vMaxScore = _mm_max_epu8(vMaxScore, _mm_srli_si128(vMaxScore, 2));
	vMaxScore = _mm_max_epu8(vMaxScore, _mm_srli_si128(vMaxScore, 1));
	MaxScore = (unsigned char) _mm_extract_epi16(vMaxScore, 0);
	if ((int) MaxScore + (int) query->bias >= 255)
		return DBL_MAX;
	return ((double) MaxScore);
}

EXPORT void swps3_freeProfileByteSSE(ProfileByte * profile) {
	free(profile);
}

