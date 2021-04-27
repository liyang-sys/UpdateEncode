#define _CRT_SECURE_NO_WARNINGS
#include "readYUV.h"
#include <stdio.h>
#include <stdlib.h>


unsigned char ***readYUV8(char* filesourse, int h, int w, int f1, int Nframe, int fmt)
{
	FILE  *fp;
	int hh = h >> 1;
	int hw = w >> 1;
	int nby = 1;
	fp = fopen(filesourse, "rb");
	if (fp == NULL) {
		printf("打开文件失败！\n");
		return NULL;
	}

	unsigned char*** ImgData = (unsigned char ***)malloc(sizeof(unsigned char **)* Nframe);
	for (int i = 0; i < Nframe; i++)
	{
		ImgData[i] = (unsigned char**)malloc(sizeof(unsigned char*)* h);
		for (int j = 0; j < h; j++)
		{
			ImgData[i][j] = (unsigned char*)malloc(sizeof(unsigned char)* w);
		}
	}

	for (int i = f1; i < f1 + Nframe; i++) {
		if (fmt == 444) {
			fseek(fp, (i - 1)*w*h * 3 * nby, 0);
			for (int j = 0; j < h; j++)
			{
				fread(&ImgData[i - f1][j][0], sizeof(unsigned char), w, fp);
			}
		}
		else if (fmt == 422) {
			fseek(fp, (i - 1)*w*h * 2 * nby, 0);
			for (int j = 0; j < h; j++)
			{
				fread(&ImgData[i - f1][j][0], sizeof(unsigned char), w, fp);
			}
		}
		else { //fmt=420
			fseek(fp, (i - 1)*w*h * 1.5 * nby, 0);
			for (int j = 0; j < h; j++)
			{
				fread(&ImgData[i - f1][j][0], sizeof(unsigned char), w, fp);
			}
		}
	}


	return ImgData;
}

unsigned short ***readYUV10(char* filesourse, int h, int w, int f1, int Nframe, int fmt)
{
	FILE  *fp;
	int hh = h >> 1;
	int hw = w >> 1;
	int nby = 2;
	unsigned short *u, *v;

	fp = fopen(filesourse, "rb");
	if (fp == NULL) {
		printf("打开文件失败！\n");
		return NULL;
	}

	unsigned short*** ImgData = (unsigned short ***)malloc(sizeof(unsigned short **)* Nframe);
	for (int i = 0; i < Nframe; i++)
	{
		ImgData[i] = (unsigned short**)malloc(sizeof(unsigned short*)* h);
		for (int j = 0; j < h; j++)
		{
			ImgData[i][j] = (unsigned short*)malloc(sizeof(unsigned short)* w);
		}
	}

	if (fmt == 444) {
		u = (unsigned short*)malloc(sizeof(unsigned short)*w);
		v = (unsigned short*)malloc(sizeof(unsigned short)*w);
		for (int i = f1; i < f1 + Nframe; i++) {
			fseek(fp, (i - 1)*w*h * 3 * nby, 0);
			for (int j = 0; j < h; j++)
			{
				fread(&ImgData[i - f1][j][0], sizeof(unsigned short), w, fp);
			}
			for (int j = 0; j < h; j++)
			{
				fread(u, sizeof(unsigned short), w, fp);
			}
			for (int j = 0; j < h; j++)
			{
				fread(v, sizeof(unsigned short), w, fp);
			}
		}
	}
	else if (fmt == 422) {
		u = (unsigned short*)malloc(sizeof(unsigned short)*w);
		v = (unsigned short*)malloc(sizeof(unsigned short)*w);
		for (int i = f1; i < f1 + Nframe; i++) {
			for (int j = 0; j < h; j++)
			{
				fread(&ImgData[i - f1][j][0], sizeof(unsigned short), w, fp);
			}
			for (int j = 0; j < hh; j++)
			{
				fread(u, sizeof(unsigned short), w, fp);
			}
			for (int j = 0; j < hh; j++)
			{
				fread(v, sizeof(unsigned short), w, fp);
			}
		}
	}
	else { //fmt=420
		u = (unsigned short*)malloc(sizeof(unsigned short)*hw);
		v = (unsigned short*)malloc(sizeof(unsigned short)*hw);
		for (int i = f1; i < f1 + Nframe; i++) {
			for (int j = 0; j < h; j++)
			{
				fread(&ImgData[i - f1][j][0], sizeof(unsigned short), w, fp);
			}
			for (int j = 0; j < hh; j++)
			{
				fread(u, sizeof(unsigned short), hw, fp);
			}
			for (int j = 0; j < hh; j++)
			{
				fread(v, sizeof(unsigned short), hw, fp);
			}
		}
	}

	free(u);
	free(v);

	return ImgData;
}

unsigned short ***readY4M(char* filesourse, int h, int w, int f1, int Nframe) //暂时默认fmt=420
{
	FILE  *fp;
	int hh = h >> 1;
	int hw = w >> 1;
	unsigned short *u, *v;
	short temp[3];
	int offset = 50;

	fp = fopen(filesourse, "rb");
	if (fp == NULL) {
		printf("打开文件失败！\n");
		return NULL;
	}

	unsigned short*** ImgData = (unsigned short ***)malloc(sizeof(unsigned short **)* Nframe);
	for (int i = 0; i < Nframe; i++)
	{
		ImgData[i] = (unsigned short**)malloc(sizeof(unsigned short*)* h);
		for (int j = 0; j < h; j++)
		{
			ImgData[i][j] = (unsigned short*)malloc(sizeof(unsigned short)* w);
		}
	}

	u = (unsigned short*)malloc(sizeof(unsigned short)*hw);
	v = (unsigned short*)malloc(sizeof(unsigned short)*hw);
	fseek(fp, offset, SEEK_SET);
	for (int i = f1; i < f1 + Nframe; i++) {
		for (int j = 0; j < h; j++)
		{
			fread(&ImgData[i - f1][j][0], sizeof(unsigned short), w, fp);
		}
		for (int j = 0; j < hh; j++)
		{
			fread(u, sizeof(unsigned short), hw, fp);
		}
		for (int j = 0; j < hh; j++)
		{
			fread(v, sizeof(unsigned short), hw, fp);
		}
		//fseek(fp, 6, SEEK_CUR);
		fread(temp, 2, 3, fp);
	}


	return ImgData;
}