#include "quant_rest.h"
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>

union Fabs {
	float **qcf;
	int **temp;
};

int quanEVEN_DC(int a, int b, int c, int lg, int pos, int ptv, int harr, int *m, double delta)
{
	extern unsigned char **sn;
	extern int **nc;
	extern int len;
	int i = 0, j, lth, qf, max = 0, offset1, offset0;
	extern float ***PTVData;
	extern union Fabs f1;
	lth = pos + lg;
	offset0 = 1 << ptv;
	offset1 = 4 << harr;

	memset(sn[0], 0, sizeof(uchar)*(lg));
	memset(nc[0], 0, sizeof(int)*(lg));

	for (j = pos; j < lth; j++) {
		f1.qcf[0][j] = PTVData[a][b][c];
		if (f1.qcf[0][j]) {
			sn[0][j] = (((f1.temp[0][j] & 0x80000000) >> 31) & 0x01);//取符号
			sn[0][j]++; //1表示正，2表示负，0表示零
		}
		f1.temp[0][j] &= 0x7fffffff; //将符号位置零求绝对值
		qf = (int)(f1.qcf[0][j] * delta);
		if (qf > max) {
			max = qf;
		}
		nc[0][qf]++;//统计nc
		switch (m[i])
		{
		case 3:
			c += offset0;
			break;
		case -3:
			c -= offset0;
			break;
		case 2:
			b += offset0;
			break;
		case -2:
			b -= offset0;
			break;
		case 1:
			a += offset1;
			break;
		case -1:
			a -= offset1;
			break;
		default:
			break;
		}i++;
	}
	max++;
	return max;
}

int quanEVEN8(int *max, int a, int b, int lg, int pos, int harr, int *m, double delta0) //b,a分别为起点的行和深,量化8个子带
{
	extern unsigned char **sn;
	int i = 0, j, qf, c = 0, offset;
	extern int **nc;
	extern int len;
	extern float ***VideoData;
	extern union Fabs f1;
	lg += pos;
	offset = 4 << harr;

	for (int i = 0; i < 8; i++) {
		memset(sn[i], 0, sizeof(uchar)*(lg));
		memset(nc[i], 0, sizeof(int)*(lg));
	}

	for (j = pos; j < lg; j++) {
		/*第1个子带*/
		f1.qcf[0][j] = VideoData[a][b][c];
		if (f1.qcf[0][j]) {
			sn[0][j] = (((f1.temp[0][j] & 0x80000000) >> 31) & 0x01);//取符号
			sn[0][j]++; //1表示正，2表示负，0表示零
		}
		f1.temp[0][j] &= 0x7fffffff;		//通过将符号位置零来求绝对值
		qf = (int)(f1.qcf[0][j] * delta0);	//量化
		if (qf > max[0]) {//求出量化后系数的最大值
			max[0] = qf;
		}
		nc[0][qf]++; //求nc

		/*第2个子带*/
		f1.qcf[1][j] = VideoData[a][b][c + 1];
		if (f1.qcf[1][j]) {
			sn[1][j] = (((f1.temp[1][j] & 0x80000000) >> 31) & 0x01);//取符号
			sn[1][j]++; //1表示正，2表示负，0表示零
		}
		f1.temp[1][j] &= 0x7fffffff;		//通过将符号位置零来求绝对值
		qf = (int)(f1.qcf[1][j] * delta0);	//量化
		if (qf > max[1]) {//求出量化后系数的最大值
			max[1] = qf;
		}
		nc[1][qf]++; //求nc

		/*3*/
		f1.qcf[2][j] = VideoData[a][b][c + 2];
		if (f1.qcf[2][j]) {
			sn[2][j] = (((f1.temp[2][j] & 0x80000000) >> 31) & 0x01);//取符号
			sn[2][j]++; //1表示正，2表示负，0表示零
		}
		f1.temp[2][j] &= 0x7fffffff;		//通过将符号位置零来求绝对值
		qf = (int)(f1.qcf[2][j] * delta0);	//量化
		if (qf > max[2]) {//求出量化后系数的最大值
			max[2] = qf;
		}
		nc[2][qf]++; //求nc

		/*4*/
		f1.qcf[3][j] = VideoData[a][b][c + 3];
		if (f1.qcf[3][j]) {
			sn[3][j] = (((f1.temp[3][j] & 0x80000000) >> 31) & 0x01);//取符号
			sn[3][j]++; //1表示正，2表示负，0表示零
		}
		f1.temp[3][j] &= 0x7fffffff;		//通过将符号位置零来求绝对值
		qf = (int)(f1.qcf[3][j] * delta0);	//量化
		if (qf > max[3]) {//求出量化后系数的最大值
			max[3] = qf;
		}
		nc[3][qf]++; //求nc

		/*5*/
		f1.qcf[4][j] = VideoData[a][b][c + 4];
		if (f1.qcf[4][j]) {
			sn[4][j] = (((f1.temp[4][j] & 0x80000000) >> 31) & 0x01);//取符号
			sn[4][j]++; //1表示正，2表示负，0表示零
		}
		f1.temp[4][j] &= 0x7fffffff;		//通过将符号位置零来求绝对值
		qf = (int)(f1.qcf[4][j] * delta0);	//量化
		if (qf > max[4]) {//求出量化后系数的最大值
			max[4] = qf;
		}
		nc[4][qf]++; //求nc

		/*6*/
		f1.qcf[5][j] = VideoData[a][b][c + 5];
		if (f1.qcf[5][j]) {
			sn[5][j] = (((f1.temp[5][j] & 0x80000000) >> 31) & 0x01);//取符号
			sn[5][j]++; //1表示正，2表示负，0表示零
		}
		f1.temp[5][j] &= 0x7fffffff;		//通过将符号位置零来求绝对值
		qf = (int)(f1.qcf[5][j] * delta0);	//量化
		if (qf > max[5]) {//求出量化后系数的最大值
			max[5] = qf;
		}
		nc[5][qf]++; //求nc

		/*7*/
		f1.qcf[6][j] = VideoData[a][b][c + 6];
		if (f1.qcf[6][j]) {
			sn[6][j] = (((f1.temp[6][j] & 0x80000000) >> 31) & 0x01);//取符号
			sn[6][j]++; //1表示正，2表示负，0表示零
		}
		f1.temp[6][j] &= 0x7fffffff;		//通过将符号位置零来求绝对值
		qf = (int)(f1.qcf[6][j] * delta0);	//量化
		if (qf > max[6]) {//求出量化后系数的最大值
			max[6] = qf;
		}
		nc[6][qf]++; //求nc

		/*8*/
		f1.qcf[7][j] = VideoData[a][b][c + 7];
		if (f1.qcf[7][j]) {
			sn[7][j] = (((f1.temp[7][j] & 0x80000000) >> 31) & 0x01);//取符号
			sn[7][j]++; //1表示正，2表示负，0表示零
		}
		f1.temp[7][j] &= 0x7fffffff;		//通过将符号位置零来求绝对值
		qf = (int)(f1.qcf[7][j] * delta0);	//量化
		if (qf > max[7]) {//求出量化后系数的最大值
			max[7] = qf;
		}
		nc[7][qf]++; //求nc

		/*保存sn*/
		switch (m[i])
		{
		case 1:
			a += offset;
			break;
		case -1:
			a -= offset;
			break;
		case 2:
			b += 8;
			break;
		case -2:
			b -= 8;
			break;
		case 3:
			c += 8;
			break;
		case -3:
			c -= 8;
			break;
		default:
			break;
		}i++;
	}
	for (int i = 0; i < 8; i++)
		max[i]++;

	return lg;
}

int quanEVEN7(int *max, int a, int lg, int pos, int harr, int *m, double delta0) //a为起点的深
{
	extern unsigned char **sn;
	extern int **nc;
	extern int len;
	int i = 0, j, qf, c = 1, b = 0, offset;
	extern float ***VideoData;
	extern union Fabs f1;
	lg += pos;
	offset = 4 << harr;

	for (int i = 1; i < 8; i++) {
		memset(sn[i], 0, sizeof(uchar)*(lg));
		memset(nc[i], 0, sizeof(int)*(lg));
	}

	for (j = pos; j < lg; j++) {
		/*1*/
		f1.qcf[1][j] = VideoData[a][b][c];
		if (f1.qcf[1][j]) {
			sn[1][j] = (((f1.temp[1][j] & 0x80000000) >> 31) & 0x01);//取符号
			sn[1][j]++; //1表示正，2表示负，0表示零
		}
		f1.temp[1][j] &= 0x7fffffff;		//通过将符号位置零来求绝对值
		qf = (int)(f1.qcf[1][j] * delta0);	//量化
		if (qf > max[0]) {//求出量化后系数的最大值
			max[0] = qf;
		}
		nc[1][qf]++; //求nc

		/*2*/
		f1.qcf[2][j] = VideoData[a][b][c + 1];
		if (f1.qcf[2][j]) {
			sn[2][j] = (((f1.temp[2][j] & 0x80000000) >> 31) & 0x01);//取符号
			sn[2][j]++; //1表示正，2表示负，0表示零
		}
		f1.temp[2][j] &= 0x7fffffff;		//通过将符号位置零来求绝对值
		qf = (int)(f1.qcf[2][j] * delta0);	//量化
		if (qf > max[1]) {//求出量化后系数的最大值
			max[1] = qf;
		}
		nc[2][qf]++; //求nc

		/*3*/
		f1.qcf[3][j] = VideoData[a][b][c + 2];
		if (f1.qcf[3][j]) {
			sn[3][j] = (((f1.temp[3][j] & 0x80000000) >> 31) & 0x01);//取符号
			sn[3][j]++; //1表示正，2表示负，0表示零
		}
		f1.temp[3][j] &= 0x7fffffff;		//通过将符号位置零来求绝对值
		qf = (int)(f1.qcf[3][j] * delta0);	//量化
		if (qf > max[2]) {//求出量化后系数的最大值
			max[2] = qf;
		}
		nc[3][qf]++; //求nc

		/*4*/
		f1.qcf[4][j] = VideoData[a][b][c + 3];
		if (f1.qcf[4][j]) {
			sn[4][j] = (((f1.temp[4][j] & 0x80000000) >> 31) & 0x01);//取符号
			sn[4][j]++; //1表示正，2表示负，0表示零
		}
		f1.temp[4][j] &= 0x7fffffff;		//通过将符号位置零来求绝对值
		qf = (int)(f1.qcf[4][j] * delta0);	//量化
		if (qf > max[3]) {//求出量化后系数的最大值
			max[3] = qf;
		}
		nc[4][qf]++; //求nc

		/*5*/
		f1.qcf[5][j] = VideoData[a][b][c + 4];
		if (f1.qcf[5][j]) {
			sn[5][j] = (((f1.temp[5][j] & 0x80000000) >> 31) & 0x01);//取符号
			sn[5][j]++; //1表示正，2表示负，0表示零
		}
		f1.temp[5][j] &= 0x7fffffff;		//通过将符号位置零来求绝对值
		qf = (int)(f1.qcf[5][j] * delta0);	//量化
		if (qf > max[4]) {//求出量化后系数的最大值
			max[4] = qf;
		}
		nc[5][qf]++; //求nc

		/*6*/
		f1.qcf[6][j] = VideoData[a][b][c + 5];
		if (f1.qcf[6][j]) {
			sn[6][j] = (((f1.temp[6][j] & 0x80000000) >> 31) & 0x01);//取符号
			sn[6][j]++; //1表示正，2表示负，0表示零
		}
		f1.temp[6][j] &= 0x7fffffff;		//通过将符号位置零来求绝对值
		qf = (int)(f1.qcf[6][j] * delta0);	//量化
		if (qf > max[5]) {//求出量化后系数的最大值
			max[5] = qf;
		}
		nc[6][qf]++; //求nc

		/*7*/
		f1.qcf[7][j] = VideoData[a][b][c + 6];
		if (f1.qcf[7][j]) {
			sn[7][j] = (((f1.temp[7][j] & 0x80000000) >> 31) & 0x01);//取符号
			sn[7][j]++; //1表示正，2表示负，0表示零
		}
		f1.temp[7][j] &= 0x7fffffff;		//通过将符号位置零来求绝对值
		qf = (int)(f1.qcf[7][j] * delta0);	//量化
		if (qf > max[6]) {//求出量化后系数的最大值
			max[6] = qf;
		}
		nc[7][qf]++; //求nc

		/*保存sn*/
		switch (m[i])
		{
		case 3:
			c += 8;
			break;
		case -3:
			c -= 8;
			break;
		case 2:
			b += 8;
			break;
		case -2:
			b -= 8;
			break;
		case 1:
			a += offset;
			break;
		case -1:
			a -= offset;
			break;
		default:
			break;
		}i++;
	}
	for (int i = 0; i < 7; i++)
		max[i]++;

	return lg;
}

int quanEVEN2(int len, int pos, float *absA, double delta)  //只进行量化，不求nc，runs
{
	int i, num1 = 0;
	int ctr[2], intpart;
	double sum1 = 0.0;
	double temp;
	len += pos;
	for (i = pos; i < len; i++) {
		temp = (double)*(absA + i);
		if (temp) {
			intpart = (int)(temp / delta);
			if ((temp / delta) - (double)intpart) //如果不能被整除
				*(absA + i) = intpart + 1; //+1 达到向上取整效果
			else
				*(absA + i) = intpart;
		}
		else
			*(absA + i) = 1;

		if (*(absA + i) == 1) {
			num1++;
			sum1 += temp;
		}
	}
	if (num1 > 0) {
		ctr[0] = (int)(((sum1 / delta) / num1) * 512 - 191.5);
	}
	else ctr[0] = 63;
	return ctr[0];
}

QT quanTHD(uchar *sgn, int len, int pos, float *fabsA, double delta)
{
	QT qtTHD;
	int i, j = 0, *temp;
	extern double *absA2;
	extern double *absA;
	qtTHD.runs = (int*)calloc(len, sizeof(int));
	if (!qtTHD.runs) {
		printf("创建动态数组runs失败！\n");
		exit(1);
	}
	len += pos;
	for (i = pos; i < len; i++)
	{
		qtTHD.runs[j]++;		//求runs
		absA2[i] = fabsA[i] * fabsA[i];
		absA[i] = fabsA[i];
		fabsA[i] = (int)(fabsA[i] * delta + 0.45); //量化
		if (fabsA[i])
			j++;
		else
			sgn[i] = 0;
	}
	if (!fabsA[len - 1])
		j++;
	temp = (int*)realloc(qtTHD.runs, j * sizeof(int));
	if (temp == NULL) {
		printf("重新分配内存失败！\n");
		exit(1);
	}
	qtTHD.runs = temp;
	qtTHD.lg = j;
	return qtTHD;
}

int quanTHD2(float *absA, double delta)  //只进行量化
{
	extern int len;
	int i;
	for (i = 0; i < len; i++)
	{
		absA[i] = (int)(absA[i] * delta + 0.48); //量化
	}
	return 0;
}

void quant3(int *max, int H, int W, int a, int b, int c, int lg, int harr, int mul, int *m, int overflow, double delta0)
{
	extern unsigned char **sn;
	int i = 0, j, qf, offset, cnt, tempb, tempc;
	extern int **nc;
	extern int len;
	extern float ***VideoData;
	extern union Fabs f1;
	offset = 4 << harr;
	cnt = 1 << mul;
	tempc = c; tempb = b;
	for (int k = 1; k < 4; k++) {
		memset(sn[k], 0, lg * sizeof(uchar));
		memset(nc[k], 0, lg * sizeof(int));
	}

	for (j = 0; j < lg; j++) {
		/*第1个子带*/
		f1.qcf[1][j] = VideoData[a][b][c];
		if (f1.qcf[1][j]) {
			sn[1][j] = (((f1.temp[1][j] & 0x80000000) >> 31) & 0x01);//取符号
			sn[1][j]++; //1表示正，2表示负，0表示零
		}
		f1.temp[1][j] &= 0x7fffffff;		//通过将符号位置零来求绝对值
		qf = (int)(f1.qcf[1][j] * delta0);	//量化
		if (qf > max[0]) {//求出量化后系数的最大值
			max[0] = qf;
		}
		nc[1][qf]++; //求nc

		/*第2个子带*/
		f1.qcf[2][j] = VideoData[a][b + cnt][c];
		if (f1.qcf[2][j]) {
			sn[2][j] = (((f1.temp[2][j] & 0x80000000) >> 31) & 0x01);//取符号
			sn[2][j]++; //1表示正，2表示负，0表示零
		}
		f1.temp[2][j] &= 0x7fffffff;		//通过将符号位置零来求绝对值
		qf = (int)(f1.qcf[2][j] * delta0);	//量化
		if (qf > max[1]) {//求出量化后系数的最大值
			max[1] = qf;
		}
		nc[2][qf]++; //求nc

		/*3*/
		f1.qcf[3][j] = VideoData[a][b + cnt][c - cnt];
		if (f1.qcf[3][j]) {
			sn[3][j] = (((f1.temp[3][j] & 0x80000000) >> 31) & 0x01);//取符号
			sn[3][j]++; //1表示正，2表示负，0表示零
		}
		f1.temp[3][j] &= 0x7fffffff;		//通过将符号位置零来求绝对值
		qf = (int)(f1.qcf[3][j] * delta0);	//量化
		if (qf > max[2]) {//求出量化后系数的最大值
			max[2] = qf;
		}
		nc[3][qf]++; //求nc

		/*保存sn*/
		switch (m[i])
		{
		case 1:
			a += offset;
			break;
		case -1:
			a -= offset;
			break;
		case 2:
			b += 8;
			break;
		case -2:
			b -= 8;
			break;
		case 3:
			c += 8;
			break;
		case -3:
			c -= 8;
			break;
		default:
			break;
		}i++;

		if (overflow) {
			if (b >= H) {
				tempb++;
				b = tempb;
			}
			else if (b < 0) {
				tempb--;
				b = tempb + H - 8;
			}
			else if (c >= W) {
				tempc++;
				c = tempc;
			}
			else if (c < 0) {
				tempc--;
				c = tempc + W - 8;
			}
		}
	}
	for (int i = 0; i < 3; i++)
		max[i]++;

	return;
}

void quant_L2_3(int *max, int H, int W, int a, int b, int c, int lg, int harr, int *m, double delta0)
{
	extern unsigned char **sn;
	extern int **nc;
	extern int len;
	int i = 0, j, qf, tempc, tempb, offset;
	extern float ***VideoData;
	extern union Fabs f1;
	offset = 4 << harr;
	tempc = c; tempb = b;
	for (int k = 1; k < 4; k++) {
		memset(sn[k], 0, lg * sizeof(uchar));
		memset(nc[k], 0, lg * sizeof(int));
	}

	for (j = 0; j < lg; j++) {
		/*1*/
		f1.qcf[1][j] = VideoData[a][b][c];
		if (f1.qcf[1][j]) {
			sn[1][j] = (((f1.temp[1][j] & 0x80000000) >> 31) & 0x01);//取符号
			//sn[1][j] = ((((signed char*)&f1.qcf[1][j])[3] >> 7) & 0x01);//取符号
			sn[1][j]++; //1表示正，2表示负，0表示零
		}
		f1.temp[1][j] &= 0x7fffffff;		//通过将符号位置零来求绝对值
		qf = (int)(f1.qcf[1][j] * delta0);	//量化
		if (qf > max[0]) {//求出量化后系数的最大值
			max[0] = qf;
		}
		nc[1][qf]++; //求nc

		/*2*/
		f1.qcf[2][j] = VideoData[a][b][c + 2];
		if (f1.qcf[2][j]) {
			sn[2][j] = (((f1.temp[2][j] & 0x80000000) >> 31) & 0x01);//取符号
			sn[2][j]++; //1表示正，2表示负，0表示零
		}
		f1.temp[2][j] &= 0x7fffffff;		//通过将符号位置零来求绝对值
		qf = (int)(f1.qcf[2][j] * delta0);	//量化
		if (qf > max[1]) {//求出量化后系数的最大值
			max[1] = qf;
		}
		nc[2][qf]++; //求nc

		/*3*/
		f1.qcf[3][j] = VideoData[a][b][c + 4];
		if (f1.qcf[3][j]) {
			sn[3][j] = (((f1.temp[3][j] & 0x80000000) >> 31) & 0x01);//取符号
			sn[3][j]++; //1表示正，2表示负，0表示零
		}
		f1.temp[3][j] &= 0x7fffffff;		//通过将符号位置零来求绝对值
		qf = (int)(f1.qcf[3][j] * delta0);	//量化
		if (qf > max[2]) {//求出量化后系数的最大值
			max[2] = qf;
		}
		nc[3][qf]++; //求nc

		/*保存sn*/
		switch (m[i])
		{
		case 3:
			c += 8;
			break;
		case -3:
			c -= 8;
			break;
		case 2:
			b += 8;
			break;
		case -2:
			b -= 8;
			break;
		case 1:
			a += offset;
			break;
		case -1:
			a -= offset;
			break;
		default:
			break;
		}i++;

		if (b >= H)
			b = tempb + 1;
		else if (b < 0)
			b = tempb + H - 8;
		else if (c >= W)
			c = tempc + 1;
		else if (c < 0)
			c = tempc + W - 8;
	}
	for (int i = 0; i < 3; i++)
		max[i]++;

	return;
}

void quant_L2_4(int *max, int H, int W, int a, int b, int c, int lg, int harr, int *m, double delta0)
{
	extern unsigned char **sn;
	extern int **nc;
	extern int len;
	int i = 0, j, qf, tempc, tempb, offset;
	extern float ***VideoData;
	extern union Fabs f1;
	offset = 4 << harr;
	tempc = c; tempb = b;
	for (int k = 1; k < 5; k++) {
		memset(sn[k], 0, lg * sizeof(uchar));
		memset(nc[k], 0, lg * sizeof(int));
	}

	for (j = 0; j < lg; j++) {
		/*1*/
		f1.qcf[1][j] = VideoData[a][b][c];
		if (f1.qcf[1][j]) {
			sn[1][j] = (((f1.temp[1][j] & 0x80000000) >> 31) & 0x01);//取符号
			sn[1][j]++; //1表示正，2表示负，0表示零
		}
		f1.temp[1][j] &= 0x7fffffff;		//通过将符号位置零来求绝对值
		qf = (int)(f1.qcf[1][j] * delta0);	//量化
		if (qf > max[0]) {//求出量化后系数的最大值
			max[0] = qf;
		}
		nc[1][qf]++; //求nc

		/*2*/
		f1.qcf[2][j] = VideoData[a][b][c + 2];
		if (f1.qcf[2][j]) {
			sn[2][j] = (((f1.temp[2][j] & 0x80000000) >> 31) & 0x01);//取符号
			sn[2][j]++; //1表示正，2表示负，0表示零
		}
		f1.temp[2][j] &= 0x7fffffff;		//通过将符号位置零来求绝对值
		qf = (int)(f1.qcf[2][j] * delta0);	//量化
		if (qf > max[1]) {//求出量化后系数的最大值
			max[1] = qf;
		}
		nc[2][qf]++; //求nc

		/*3*/
		f1.qcf[3][j] = VideoData[a][b][c + 4];
		if (f1.qcf[3][j]) {
			sn[3][j] = (((f1.temp[3][j] & 0x80000000) >> 31) & 0x01);//取符号
			sn[3][j]++; //1表示正，2表示负，0表示零
		}
		f1.temp[3][j] &= 0x7fffffff;		//通过将符号位置零来求绝对值
		qf = (int)(f1.qcf[3][j] * delta0);	//量化
		if (qf > max[2]) {//求出量化后系数的最大值
			max[2] = qf;
		}
		nc[3][qf]++; //求nc

		/*4*/
		f1.qcf[4][j] = VideoData[a][b][c + 6];
		if (f1.qcf[4][j]) {
			sn[4][j] = (((f1.temp[4][j] & 0x80000000) >> 31) & 0x01);//取符号
			sn[4][j]++; //1表示正，2表示负，0表示零
		}
		f1.temp[4][j] &= 0x7fffffff;		//通过将符号位置零来求绝对值
		qf = (int)(f1.qcf[4][j] * delta0);	//量化
		if (qf > max[3]) {//求出量化后系数的最大值
			max[3] = qf;
		}
		nc[4][qf]++; //求nc

		/*保存sn*/
		switch (m[i])
		{
		case 3:
			c += 8;
			break;
		case -3:
			c -= 8;
			break;
		case 2:
			b += 8;
			break;
		case -2:
			b -= 8;
			break;
		case 1:
			a += offset;
			break;
		case -1:
			a -= offset;
			break;
		default:
			break;
		}i++;

		if (b >= H)
			b = tempb + 1;
		else if (b < 0)
			b = tempb + H - 8;
		else if (c >= W)
			c = tempc + 1;
		else if (c < 0)
			c = tempc + W - 8;
	}
	for (int i = 0; i < 4; i++)
		max[i]++;

	return;
}

