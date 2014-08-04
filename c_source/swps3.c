/** \file swps3.c
 *
 * Main procedure and multi-threading code.
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

#include "swps3.h"
#include "matrix.h"
#include "fasta.h"
#include "DynProgr_scalar.h"
#if defined(__SSE2__)
#include "DynProgr_sse_byte.h"
#include "DynProgr_sse_short.h"
#endif
#if defined(__ALTIVEC__)
#include "DynProgr_altivec.h"
#endif
#if defined(__PS3)
#include "DynProgr_PPU.h"
#endif
#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <sys/select.h>
#include <sys/time.h>
#include <unistd.h>
#include <string.h>
#include <assert.h>
#include <float.h>
#include "Python_extension.h"

void normalizeSequence(char* seq, int seqLen) {
	int i;
	for (i = 0; i < seqLen; ++i) {
		seq[i] -= 'A';
	}
}

double ceil(double val) {
	short tmp = (short) val;
	if(tmp < val) {
		tmp++;
	}
	return tmp;
}

void readDoubleMatrix(DMatrix matrix, const char* fileLoc, Options* options) {
	char line[40000];
	int i, curr_off, offset = 0;

	FILE* f = fopen(fileLoc, "r");

	fgets(line , 40000, f);
	sscanf(line, "%lf", &(options->gapOpen));

	fgets(line , 40000, f);
	sscanf(line, "%lf", &(options->gapExt));

	fgets(line , 40000, f);
	for (i = 0; i < MATRIX_DIM*MATRIX_DIM; ++i) {
		sscanf(line + offset, "%lf%n", &matrix[i], &curr_off);
		offset += curr_off;
	}

	fclose(f);
}

void readShortMatrix(SMatrix matrix, const char* fileLoc, Options* options) {
	double scaleFactor = 65535.0 / options->threshold;
	double dMatrix[MATRIX_DIM*MATRIX_DIM];
	int i;

	readDoubleMatrix(dMatrix, fileLoc, options);
	options->gapOpen = (short) ceil(options->gapOpen * scaleFactor);
	options->gapExt = (short) ceil(options->gapExt * scaleFactor);

	for (i = 0; i < MATRIX_DIM*MATRIX_DIM; ++i) {
		matrix[i] = (short) ceil(dMatrix[i] * scaleFactor);
	}
}

int main( int argc, char * argv[] ){

	int16_t sMatrix[MATRIX_DIM*MATRIX_DIM] __ALIGNED__;
	int i, j;

	Options options;
	options.threshold = 306.896691;

	readShortMatrix(sMatrix, "/home/machine/repos/students/2014_Ferenc_Galko_SWPS3_PY/swps3_python_extended/test/data/matrices/C_compatible/1263.dat", &options);

	double shortFactor = 65535.0 / options.threshold;

	char query[] = "PISRIDNNKILGNTGIISVTIGVIIFKDLHAKVLLGLHGFWKFIYYYDGLDVVLTVLRDLKGTTNTSDICKHKMSIAEGQFDSAAIQGCEKWILRIEIVTDLILTRSAVLCEDNSFDSRMPGFLLIAIIWANSYEIGCSYSQAASTLIARGYASFDAAVRSRIIQPTIKMEGNNDAQEPLKQNVVWQSQSFCVCFGRLPGTAESYSTCDFLLTGTELPHTFQTIRCDIFGTDKPFTNFDGAGRFAYPNFNPFGGALRNLSIGEVNHITIDHASKIEPSVKGNLITYYILEKKGFFPDGCLLSLLVDPLFLLSSPSEIKTVNLYSAKTRCTSSNAEMPIVVSIGKEGANEYTLIHLSFYVPAWRAGEYRLCSALEFTQFENSYWAHYIVTDIAADLAETQANASNGDRQEKQVGTRLMVLKAKGLTEPTASQASEPRENFFPEGKLRLSQSAAVAVGMLITVVDMTAKYGCYGETNFVVRAQLSTLLQGFPGILKIHVSVAEKCIVGIICATLKKGYLHQTAVETLPRPYRYKANARANKYQELCRLLRKATPDGKLQLFIPLAVVIWFQPTVPHEARRTNLCFVKVKLLPRVERDLCDSVQIYEWFYTQPWSIRTARPGVDPSGPEASEQWKWLDFPDFCIVLACKVQMIVHIHYKDWPIICHPFDGRKQVVDDKTMYEYQVACLLVEGLDRPERKQGQFKWRFVQYGKQNPLLPQPALVGGLGIASRFEPPTHIQADQDLTPSDTIMKTSAEAPRGDVIQYVGKEGLPTDNTWIVVQVFFDTPGQWVGAARAMPTKYLS";
	char db[] = "PISRVENNKILANTGHISVTIGCIILKELHGPPHGLSSTTHAKVLRGLHGFWKFIYYLDGLDVVMTPLRNTKGVTNTSDHCKHKMAIAEGRFDSAVIQGCEKWILRIEIVWDSILTRPAVLCLDKSFDSRMPGFLLIRIIWANSYEISFSYSNETASTLIANGYATFDAATRSRIIQPTIKMEGNNDAQKPLNQNVVWEKQTFCVCFGRLPGTAESFSTCDFLLTGEELPNSLQTIRCDLFGTDKPFTNFDGAGRFATPSFNEFGGALRNLSVGSVNHITIEHASEIEPTVKGNWVTYYVLEKKGFFPTACLLSILTDPAYLLTSPSEYRVINLYSPRTRCTSSNAEMPIVVAIGKEGAEDYTLIHLSFYVPAWRAGEYRLCSSLEFTEFSNSYWAHYIVTDIQAKRAETQANASNGKRQEKQKGTRLMVLKAKGATEPTATDQADEPRENFFPEGSIRLSQAAAVAVGHLFTVCDMTARYGCYGETNFVVRAQWSRLTQGFPGILRIGVSVAEKCIVGIICAALKPKGLLHQTAVERLPLPYRYKANARANDYQELCKLLRKSTPDGKLQLFIGPAAVITFQPHEARRTNLCFVKVKLLPRVERDACDKILIYTWFYAEPWSIRTGRPASGPEASEQYKWLDFPDFAIVLACKVQMVVHIHYQDWPINCHPFDGRKQVMDDKTMYQYQVACVLVEGLDHPQRVQGEFKWKMIQYGKTDPLLPQPSLVGGLGIASRFEPPTHIQADQDLTPTDSIMRTSAEAPRGDIIQMVGKDVFFDTPGQWVGAGRALRTKYLK";
	int ql = strlen(query);
	int dbl = strlen(db);

	normalizeSequence(query, ql);
	normalizeSequence(db, dbl);

	ProfileShort* profile = c_create_profile_short_sse(query, ql, sMatrix);

	double score = c_align_profile_short_sse(profile, db, dbl, options.gapOpen, options.gapExt, options.threshold);

	score /= shortFactor;

	printf("SHORT score: %.5f", score);

	swps3_freeProfileShortSSE(profile);

	return 0;
}
