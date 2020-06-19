#include "transform.h"

ISAtype type;
int alpha;

unsigned char Y[pic_w * pic_h];
unsigned char U[pic_w * pic_h / 4];
unsigned char V[pic_w * pic_h / 4];

//unsigned char R[pic_w * pic_h];
//unsigned char G[pic_w * pic_h];
//unsigned char B[pic_w * pic_h];

unsigned char newY[pic_w * pic_h];
unsigned char newU[pic_w * pic_h / 4];
unsigned char newV[pic_w * pic_h / 4];

void basic_transform()
{
	unsigned char R, G, B;

	memset(newY, 0, sizeof(newY));
	memset(newU, 0, sizeof(newU));
	memset(newV, 0, sizeof(newV));

	for(int i = 0; i < pic_h; ++i)
	{
		for(int j = 0; j < pic_w; ++j)
		{
			double y = Y[i * pic_w + j], u = U[(i/2) * (pic_w/2) + (j/2)], v = V[(i/2) * (pic_w/2) + (j/2)];
			double r = 1.164383 * (y - 16) + 1.596027 * (v - 128);
			double b = 1.164383 * (y - 16) + 2.017232 * (u - 128);
			double g = 1.164383 * (y - 16) - 0.391762 * (u - 128) - 0.812968 * (v - 128);
			
			R = r < 0 ? 0 : (r > 255 ? 255 : r);
			G = g < 0 ? 0 : (g > 255 ? 255 : g);
			B = b < 0 ? 0 : (b > 255 ? 255 : b);

			R = alpha / 256.0 * R;
			G = alpha / 256.0 * G;
			B = alpha / 256.0 * B;

			y = 0.256788 * R + 0.504129 * G + 0.097906 * B + 16;
			u = -0.148223 * R - 0.290993 * G + 0.439216 * B + 128;
			v = 0.439216 * R - 0.367788 * G - 0.071427 * B + 128;

			newY[i * pic_w + j] = y;
			newU[(i/2) * (pic_w/2) + (j/2)] += u / 4.0;
			newV[(i/2) * (pic_w/2) + (j/2)] += v / 4.0;
		}
	}

}

void mmx_transform()
{
	memset(newY, 0, sizeof(newY));
	memset(newU, 0, sizeof(newU));
	memset(newV, 0, sizeof(newV));
	
	for(int i = 0; i < pic_h; ++i)
	{
		for(int j = 0; j < pic_w; j += 4)
		{
			MMXunion y, u, v;
			for(int k = 0; k < 4; ++k)
			{
				y.sh[k] = Y[i * pic_w + j + k];
				u.sh[k] = U[(i/2) * (pic_w/2) + (j + k)/2];
				v.sh[k] = V[(i/2) * (pic_w/2) + (j + k)/2];
			}
			
			__m64 i16 = _mm_set1_pi16(16), i128 = _mm_set1_pi16(128);
			y.m = _mm_sub_pi16(y.m, i16);
			u.m = _mm_sub_pi16(u.m, i128);
			v.m = _mm_sub_pi16(v.m, i128);

			__m64 param1, param2, param3, param4, param5;
			param1 = _mm_set1_pi16((short)(1.164383 * 256));
			param2 = _mm_set1_pi16((short)(1.596027 * 256));
			param3 = _mm_set1_pi16((short)(2.017232 * 256));
			param4 = _mm_set1_pi16((short)(0.391762 * 256));
			param5 = _mm_set1_pi16((short)(0.812968 * 256));

			__m64 r, g, b;
			r = g = b = _mm_mullo_pi16(y.m, param1);
			r = _mm_adds_pi16(r, _mm_mullo_pi16(v.m, param2));
			b = _mm_adds_pi16(b, _mm_mullo_pi16(u.m, param3));
			g = _mm_subs_pi16(g, _mm_add_pi16(_mm_mullo_pi16(u.m, param4), _mm_mullo_pi16(v.m, param5)));

			r = _mm_srli_pi16(r, 8);
			g = _mm_srli_pi16(g, 8);
			b = _mm_srli_pi16(b, 8);

			param1 = _mm_set1_pi16((short)alpha);
			r = _mm_srli_pi16(_mm_mullo_pi16(r, param1), 8);
			g = _mm_srli_pi16(_mm_mullo_pi16(g, param1), 8);
			b = _mm_srli_pi16(_mm_mullo_pi16(b, param1), 8);

			param1 = _mm_set1_pi16((short)(0.256788 * 256));
			param2 = _mm_set1_pi16((short)(0.504129 * 256));
			param3 = _mm_set1_pi16((short)(0.097906 * 256));
			y.m = _mm_add_pi16(_mm_mullo_pi16(r, param1), _mm_mullo_pi16(g, param2));
			y.m = _mm_add_pi16(y.m, _mm_mullo_pi16(b, param3));
			y.m = _mm_srai_pi16(y.m, 8);
			y.m = _mm_add_pi16(y.m, i16);

			param1 = _mm_set1_pi16((short)(0.148223 * 256));
			param2 = _mm_set1_pi16((short)(0.290993 * 256));
			param3 = _mm_set1_pi16((short)(0.439216 * 256));
			u.m = _mm_sub_pi16(_mm_mullo_pi16(b, param3), _mm_mullo_pi16(r, param1));
			u.m = _mm_sub_pi16(u.m, _mm_mullo_pi16(g, param2));
			u.m = _mm_srai_pi16(u.m, 8);
			u.m = _mm_add_pi16(u.m, i128);

			param1 = _mm_set1_pi16((short)(0.439216 * 256));
			param2 = _mm_set1_pi16((short)(0.367788 * 256));
			param3 = _mm_set1_pi16((short)(0.071427 * 256));
			v.m = _mm_sub_pi16(_mm_mullo_pi16(r, param1), _mm_mullo_pi16(g, param2));
			v.m = _mm_sub_pi16(v.m, _mm_mullo_pi16(b, param3));
			v.m = _mm_srai_pi16(v.m, 8);
			v.m = _mm_add_pi16(v.m, i128);
			
			for(int k = 0; k < 4; ++k)
			{
				newY[i * pic_w + j + k] = y.sh[k];
				newU[(i/2) * (pic_w/2) + (j + k)/2] += u.sh[k] / 4.0;
				newV[(i/2) * (pic_w/2) + (j + k)/2] += v.sh[k] / 4.0;
			}
		}
	}
}

void sse2_transform()
{
	memset(newY, 0, sizeof(newY));
	memset(newU, 0, sizeof(newU));
	memset(newV, 0, sizeof(newV));
	
	for(int i = 0; i < pic_h; ++i)
	{
		for(int j = 0; j < pic_w; j += 8)
		{
			SSE2union y, u, v;
			for(int k = 0; k < 8; ++k)
			{
				y.sh[k] = Y[i * pic_w + j + k];
				u.sh[k] = U[(i/2) * (pic_w/2) + (j + k)/2];
				v.sh[k] = V[(i/2) * (pic_w/2) + (j + k)/2];
			}

			__m128i i16 = _mm_set1_epi16(16), i128 = _mm_set1_epi16(128);
			y.m = _mm_sub_epi16(y.m, i16);
			u.m = _mm_sub_epi16(u.m, i128);
			v.m = _mm_sub_epi16(v.m, i128);

			__m128i param1, param2, param3, param4, param5;
			param1 = _mm_set1_epi16((short)(1.164383 * 256));
			param2 = _mm_set1_epi16((short)(1.596027 * 256));
			param3 = _mm_set1_epi16((short)(2.017232 * 256));
			param4 = _mm_set1_epi16((short)(0.391762 * 256));
			param5 = _mm_set1_epi16((short)(0.812968 * 256));

			__m128i r, g, b;
			r = g = b = _mm_mullo_epi16(y.m, param1);
			r = _mm_adds_epi16(r, _mm_mullo_epi16(v.m, param2));
			b = _mm_adds_epi16(b, _mm_mullo_epi16(u.m, param3));
			g = _mm_subs_epi16(g, _mm_add_epi16(_mm_mullo_epi16(u.m, param4), _mm_mullo_epi16(v.m, param5)));

			r = _mm_srli_epi16(r, 8);
			g = _mm_srli_epi16(g, 8);
			b = _mm_srli_epi16(b, 8);

			param1 = _mm_set1_epi16((short)alpha);
			r = _mm_srli_epi16(_mm_mullo_epi16(r, param1), 8);
			g = _mm_srli_epi16(_mm_mullo_epi16(g, param1), 8);
			b = _mm_srli_epi16(_mm_mullo_epi16(b, param1), 8);

			param1 = _mm_set1_epi16((short)(0.256788 * 256));
			param2 = _mm_set1_epi16((short)(0.504129 * 256));
			param3 = _mm_set1_epi16((short)(0.097906 * 256));
			y.m = _mm_add_epi16(_mm_mullo_epi16(r, param1), _mm_mullo_epi16(g, param2));
			y.m = _mm_add_epi16(y.m, _mm_mullo_epi16(b, param3));
			y.m = _mm_srai_epi16(y.m, 8);
			y.m = _mm_add_epi16(y.m, i16);

			param1 = _mm_set1_epi16((short)(0.148223 * 256));
			param2 = _mm_set1_epi16((short)(0.290993 * 256));
			param3 = _mm_set1_epi16((short)(0.439216 * 256));
			u.m = _mm_sub_epi16(_mm_mullo_epi16(b, param3), _mm_mullo_epi16(r, param1));
			u.m = _mm_sub_epi16(u.m, _mm_mullo_epi16(g, param2));
			u.m = _mm_srai_epi16(u.m, 8);
			u.m = _mm_add_epi16(u.m, i128);

			param1 = _mm_set1_epi16((short)(0.439216 * 256));
			param2 = _mm_set1_epi16((short)(0.367788 * 256));
			param3 = _mm_set1_epi16((short)(0.071427 * 256));
			v.m = _mm_sub_epi16(_mm_mullo_epi16(r, param1), _mm_mullo_epi16(g, param2));
			v.m = _mm_sub_epi16(v.m, _mm_mullo_epi16(b, param3));
			v.m = _mm_srai_epi16(v.m, 8);
			v.m = _mm_add_epi16(v.m, i128);
			
			for(int k = 0; k < 8; ++k)
			{
				newY[i * pic_w + j + k] = y.sh[k];
				newU[(i/2) * (pic_w/2) + (j + k)/2] += u.sh[k] / 4.0;
				newV[(i/2) * (pic_w/2) + (j + k)/2] += v.sh[k] / 4.0;
			}
		}
	}

}

void avx_transform()
{
	memset(newY, 0, sizeof(newY));
	memset(newU, 0, sizeof(newU));
	memset(newV, 0, sizeof(newV));
	
	for(int i = 0; i < pic_h; ++i)
	{
		for(int j = 0; j < pic_w; j += 16)
		{
			AVXunion y, u, v;
			for(int k = 0; k < 16; ++k)
			{
				y.sh[k] = Y[i * pic_w + j + k];
				u.sh[k] = U[(i/2) * (pic_w/2) + (j + k)/2];
				v.sh[k] = V[(i/2) * (pic_w/2) + (j + k)/2];
			}

			__m256i i16 = _mm256_set1_epi16(16), i128 = _mm256_set1_epi16(128);
			y.m = _mm256_sub_epi16(y.m, i16);
			u.m = _mm256_sub_epi16(u.m, i128);
			v.m = _mm256_sub_epi16(v.m, i128);

			__m256i param1, param2, param3, param4, param5;
			param1 = _mm256_set1_epi16((short)(1.164383 * 256));
			param2 = _mm256_set1_epi16((short)(1.596027 * 256));
			param3 = _mm256_set1_epi16((short)(2.017232 * 256));
			param4 = _mm256_set1_epi16((short)(0.391762 * 256));
			param5 = _mm256_set1_epi16((short)(0.812968 * 256));

			__m256i r, g, b;
			r = g = b = _mm256_mullo_epi16(y.m, param1);
			r = _mm256_adds_epi16(r, _mm256_mullo_epi16(v.m, param2));
			b = _mm256_adds_epi16(b, _mm256_mullo_epi16(u.m, param3));
			g = _mm256_subs_epi16(g, _mm256_add_epi16(_mm256_mullo_epi16(u.m, param4), _mm256_mullo_epi16(v.m, param5)));

			r = _mm256_srli_epi16(r, 8);
			g = _mm256_srli_epi16(g, 8);
			b = _mm256_srli_epi16(b, 8);

			param1 = _mm256_set1_epi16((short)alpha);
			r = _mm256_srli_epi16(_mm256_mullo_epi16(r, param1), 8);
			g = _mm256_srli_epi16(_mm256_mullo_epi16(g, param1), 8);
			b = _mm256_srli_epi16(_mm256_mullo_epi16(b, param1), 8);

			param1 = _mm256_set1_epi16((short)(0.256788 * 256));
			param2 = _mm256_set1_epi16((short)(0.504129 * 256));
			param3 = _mm256_set1_epi16((short)(0.097906 * 256));
			y.m = _mm256_add_epi16(_mm256_mullo_epi16(r, param1), _mm256_mullo_epi16(g, param2));
			y.m = _mm256_adds_epi16(y.m, _mm256_mullo_epi16(b, param3));
			y.m = _mm256_srai_epi16(y.m, 8);
			y.m = _mm256_add_epi16(y.m, i16);

			param1 = _mm256_set1_epi16((short)(0.148223 * 256));
			param2 = _mm256_set1_epi16((short)(0.290993 * 256));
			param3 = _mm256_set1_epi16((short)(0.439216 * 256));
			u.m = _mm256_sub_epi16(_mm256_mullo_epi16(b, param3), _mm256_mullo_epi16(r, param1));
			u.m = _mm256_sub_epi16(u.m, _mm256_mullo_epi16(g, param2));
			u.m = _mm256_srai_epi16(u.m, 8);
			u.m = _mm256_add_epi16(u.m, i128);

			param1 = _mm256_set1_epi16((short)(0.439216 * 256));
			param2 = _mm256_set1_epi16((short)(0.367788 * 256));
			param3 = _mm256_set1_epi16((short)(0.071427 * 256));
			v.m = _mm256_sub_epi16(_mm256_mullo_epi16(r, param1), _mm256_mullo_epi16(g, param2));
			v.m = _mm256_sub_epi16(v.m, _mm256_mullo_epi16(b, param3));
			v.m = _mm256_srai_epi16(v.m, 8);
			v.m = _mm256_add_epi16(v.m, i128);
			
			for(int k = 0; k < 16; ++k)
			{
				newY[i * pic_w + j + k] = y.sh[k];
				newU[(i/2) * (pic_w/2) + (j + k)/2] += u.sh[k] / 4.0;
				newV[(i/2) * (pic_w/2) + (j + k)/2] += v.sh[k] / 4.0;
			}
		}
	}
}

double transform(int A)
{
	clock_t start, end;
	start = clock();

	alpha = A;
	switch(type)
	{
	case BASIC: basic_transform(); break;
	case MMX: mmx_transform(); break;
	case SSE2: sse2_transform(); break;
	case AVX: avx_transform(); break;
	}

	end = clock();

	std::string outpath = "./alpha/" + std::to_string(alpha) + ".yuv";
	FILE* fp = fopen(outpath.c_str(), "wb");
	fwrite(newY, 1, pic_w * pic_h, fp);
	fwrite(newU, 1, pic_w * pic_h / 4, fp);
	fwrite(newV, 1, pic_w * pic_h / 4, fp);
	fclose(fp);

	return (double)(end - start) / CLOCKS_PER_SEC;
}

