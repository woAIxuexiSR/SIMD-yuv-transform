#ifndef _TRANSFORM_H_
#define _TRANSFORM_H_

#include <time.h>
#include <cstring>
#include <string>

#include <mmintrin.h>
#include <emmintrin.h>
#include <immintrin.h>

enum ISAtype
{
	BASIC, MMX, SSE2, AVX
};

union MMXunion
{
	unsigned short sh[4];
	__m64 m;
};

union SSE2union
{
	unsigned short sh[8];
	__m128i m;
};

union AVXunion
{
	unsigned short sh[16];
	__m256i m;
};

const int pic_w = 1920;
const int pic_h = 1080;

double transform(int A);

#endif
