#include <iostream>
#include <string>
#include <cstring>
#include <fstream>
#include "transform.h"

extern ISAtype type;
char* filepath;

extern unsigned char Y[pic_w * pic_h];
extern unsigned char U[pic_w * pic_h / 4];
extern unsigned char V[pic_w * pic_h / 4];

int main(int argc, char* argv[])
{
	if(argc != 3)
	{
		printf("you should use \"./trans [yuvfile] -ISAtype\" to run\n");
		exit(1);
	}

	filepath = argv[1];

	if(!strcmp(argv[2], "-BASIC")) type = BASIC;
	else if(!strcmp(argv[2], "-MMX")) type = MMX;
	else if(!strcmp(argv[2], "-SSE2")) type = SSE2;
	else if(!strcmp(argv[2], "-AVX")) type = AVX;
	else 
	{
		printf("ISAtype must be BASIC, MMX, SSE2, AVX\n");
		exit(1);
	}

	FILE* fp = fopen(filepath, "rb");
	if(!fp)
	{
		printf("can't open file %s\n", filepath);
		exit(1);
	}

	int count;
	count = fread(Y, 1, pic_w * pic_h, fp);
	count = fread(U, 1, pic_w * pic_h / 4, fp);
	count = fread(V, 1, pic_w * pic_h / 4, fp);
	fclose(fp);

	double time = 0;
	for(int i = 1; i < 255; i += 3)
	{
		time = transform(i);
	}

	printf(" total 85 pictures transform time: %lfs\n", time);
	printf("                 time per picture: %lfs\n", time / 85.0);

	/*
	FILE* test = fopen("./hh.yuv", "wb");
	fwrite(Y, 1, pic_w * pic_h, test);
	fwrite(U, 1, pic_w * pic_h / 4, test);
	fwrite(V, 1, pic_w * pic_h / 4, test);
	fclose(fp);
	*/
	return 0;
}
