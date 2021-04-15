#define _CRT_SECURE_NO_WARNINGS
#include "en_sub.h"
#include "modify_Hilbert3D.h"
#include "modify_Hilbert2D.h"
#include <stdio.h>
#include <stdlib.h>
#include "quant_rest.h"
#include "trim.h"
#include <math.h>
#include "all.h"
#include "encoding.h"
#include "quant_inv.h"
#include <trim.h>

#define ENCODE 1
#define ulong unsigned long long
#define uchar unsigned char

int quant_sub0(uchar *sgn, int max, int *nc, int lenth, int idx, double delta, float *fabsA, int *qnt0)
{
	union data {
		unsigned int a;
		uchar b[4];
	} rem;
	QT qtTHD; SFC se;
	uint tp = 1;
	int i, j, sum3 = 0, sum4 = 0, lenT = 131072, NT, lastlen, flg = 0, qctr, sumflg = 0, p = 0, qnt = 1;
	extern double *T;
	extern double *T2;
	extern double *absA;
	extern double *absA2;
	extern uchar *bin;
	extern int len;
	extern int ptr;
	double non, r61;
	double temp, min, max_xx, sum1 = 0.0, sum2 = 0.0;
	double ctr[2], x2[14], x1[15], x[29], xx[30];

	if (max >= 6) {
		qnt = 1;
		if ((max < 1200) || (delta < 0.5)) {  //因为这里的delta在外面已经进行了1/delta操作，所以是<0.5，而不>5
			non = 1.0 - (double)(nc[0] + nc[1] + nc[2] + nc[3] + nc[4]) / (double)lenth;// non=sum(nc(6:lennc))/sum(nc);
			r61 = (double)nc[5] / nc[0];
			if (max >= 30) {
				for (j = 0; j < 29; j += 2) {
					x1[(j >> 1)] = (double)nc[j] + nc[j + 1];
					if (j < 27) {
						x2[(j >> 1)] = (double)nc[j + 1] + nc[j + 2];
						x[j] = x1[(j >> 1)];
						x[j + 1] = x2[(j >> 1)];
					}
					else x[j] = x1[(j >> 1)];   //x[28]=x1[14]
				}
				min = x[0]; max_xx = 0;
				for (j = 0; j < 29; j++) {	//找出x[]中的最小值
					if (min > x[j])
						min = x[j];
				}
				if (min > 0) {
					for (j = 0; j < 28; j += 2) {
						xx[j] = x[j + 1] / x[j];
						xx[j + 1] = x[j + 2] / x[j + 1];
					}
					for (j = 15; j < 25; j++) {
						if (max_xx < xx[j]) max_xx = xx[j];  //求xx最大值
					}
					if (max_xx > 0.9) {
						for (j = 0; j < 30; j++)
							xx[j] = j + 1;
					}
				}
				else {
					for (j = 0; j < 30; j++)
						xx[j] = j + 1;
				}
			}
			else {
				for (j = 0; j < 30; j++)
					xx[j] = j + 1;
			}
			for (i = 15; i < 20; i++)//sum1为sum(xx(16:20)),sum2为sum(xx(21:25))
				sum1 += xx[i];
			for (i = 20; i < 25; i++)
				sum2 += xx[i];
			if (((non < 0.6) && (r61 < 0.71) && (sum1 < sum2)) || ((non < 0.36) && (r61 < 0.6))) {
				qctr = 1;
				if (lenth > lenT) {
					NT = (int)lenth / lenT;
					for (i = 0; i < NT; i++) {
						qtTHD = quanTHD(sgn, lenT, idx, fabsA, delta);
						if ((qtTHD.lg > 1) || ((qtTHD.lg == 1) && (fabsA[lenT + idx - 1] != 0))) {
							flg = trim_coef(sgn, qtTHD.runs, qtTHD.lg, lenT, idx, fabsA);
							sumflg += flg;
						}
						idx += lenT;
					}
					lastlen = lenth - NT * lenT;
					qtTHD = quanTHD(sgn, lastlen, idx, fabsA, delta);
					if ((qtTHD.lg > 1) || ((qtTHD.lg == 1) && (fabsA[lastlen + idx - 1] != 0))) {
						flg = trim_coef(sgn, qtTHD.runs, qtTHD.lg, lastlen, idx, fabsA);
						sumflg += flg;
					}
				}
				else {
					qtTHD = quanTHD(sgn, lenth, idx, fabsA, delta);  //cf0=sb;
					if ((qtTHD.lg > 1) || ((qtTHD.lg == 1) && (fabsA[lenth + idx - 1] != 0))) {
						flg = trim_coef(sgn, qtTHD.runs, qtTHD.lg, lenth, idx, fabsA);
						sumflg = flg;
					}
				}

				if (sumflg) { //sumflg大于0，表示max(cf0r)>0 ,bin=[11]
					p = ptr & 7;
					rem.a = 0;
					rem.b[2] = bin[ptr >> 3];
					rem.a >>= (8 - p);
					rem.b[1] = (uchar)(3 << 6);
					rem.a <<= (16 - p);
					bin[ptr >> 3] = rem.b[3];
					bin[(ptr >> 3) + 1] = rem.b[2];
					ptr += 2;
				}
				else { //bin=[0];
					qnt = 0;
					p = ptr & 7;
					rem.a = 0;
					rem.b[2] = bin[ptr >> 3];
					rem.a >>= (8 - p);
					rem.b[1] = 0;
					rem.a <<= (16 - p);
					bin[ptr >> 3] = rem.b[3];
					bin[(ptr >> 3) + 1] = rem.b[2];
					ptr += 1;
				}
			}
			else {
				for (i = 0; i < 5; i++)  sum3 += nc[i];
				for (i = 5; i < 15; i++) sum4 += nc[i];
				temp = (double)sum3 / sum4;
				if ((temp > 0.7071) || (((double)nc[0] / (double)nc[1] > 1.4) && (nc[0] > 35))) {
					lenth += idx;
					for (i = idx; i < lenth; i++)
					{
						fabsA[i] = (int)(fabsA[i] * delta + 0.48); //quanTHD量化
						if (fabsA[i] == 0)
							sgn[i] = 0;
					}
					qctr = 1;
					p = ptr & 7;
					rem.a = 0;
					rem.b[2] = bin[ptr >> 3];
					rem.a >>= (8 - p);
					rem.b[1] = (uchar)(5 << 5);
					rem.a <<= (16 - p);
					bin[ptr >> 3] = rem.b[3];
					bin[(ptr >> 3) + 1] = rem.b[2];
					ptr += 3;
				}
				else {
					ctr[0] = quanEVEN2(lenth, idx, fabsA, delta);
					if (ctr[0] > 63) ctr[0] = 63;
					if (ctr[0] < 0)  ctr[0] = 0;
					se = SFcode((uint)ctr[0] + 1, 64);
					p = ptr & 7;
					rem.a = 0;
					rem.b[1] = bin[ptr >> 3];
					rem.a >>= (16 - p);
					rem.a <<= 3;
					rem.b[0] |= (uchar)4;
					rem.a <<= se.lb;
					rem.a |= se.code;
					rem.a <<= (32 - p - 3 - se.lb);
					bin[ptr >> 3] = rem.b[3];
					bin[(ptr >> 3) + 1] = rem.b[2];
					bin[(ptr >> 3) + 2] = rem.b[1];
					bin[(ptr >> 3) + 3] = rem.b[0];
					ptr += (3 + se.lb);
					qctr = 0;
				}
			}
		}
		else {
			for (i = 0; i < 5; i++)  sum3 += nc[i];
			for (i = 5; i < 15; i++) sum4 += nc[i];
			temp = (double)sum3 / sum4;
			if (((temp > 0.7071) && max < 2600) || (((double)nc[0] / (double)nc[1] > 1.4) && nc[0] > 35)) {
				lenth += idx;
				for (i = idx; i < lenth; i++)
				{
					fabsA[i] = (int)(fabsA[i] * delta + 0.48); //quantTHD量化,相比较matlab,通过计算省去了offset
					if (!fabsA[i])
						sgn[i] = 0;
				}
				qctr = 1;
				p = ptr & 7;
				rem.a = 0;
				rem.b[2] = bin[ptr >> 3];
				rem.a >>= (8 - p);
				rem.b[1] = (uchar)(5 << 5);
				rem.a <<= (16 - p);
				bin[ptr >> 3] = rem.b[3];
				bin[(ptr >> 3) + 1] = rem.b[2];
				ptr += 3;
			}
			else {
				ctr[0] = quanEVEN2(lenth, idx, fabsA, delta);
				if (ctr[0] > 63) ctr[0] = 63;
				if (ctr[0] < 0)  ctr[0] = 0;
				se = SFcode((uint)ctr[0] + 1, 64);
				p = ptr & 7;
				rem.a = 0;
				rem.b[1] = bin[ptr >> 3];
				rem.a >>= (16 - p);
				rem.a <<= 3;
				rem.b[0] |= (uchar)4;
				rem.a <<= se.lb;
				rem.a |= se.code;
				rem.a <<= (32 - p - 3 - se.lb);
				bin[ptr >> 3] = rem.b[3];
				bin[(ptr >> 3) + 1] = rem.b[2];
				bin[(ptr >> 3) + 2] = rem.b[1];
				bin[(ptr >> 3) + 3] = rem.b[0];
				ptr += (3 + se.lb);
				qctr = 0;
			}
		}
	}
	else {
		qnt = 0;
		qctr = 1;
		p = ptr & 7;
		rem.a = 0;
		rem.b[2] = bin[ptr >> 3];
		rem.a >>= (8 - p);
		rem.b[1] = 0;
		rem.a <<= (16 - p);
		bin[ptr >> 3] = rem.b[3];
		bin[(ptr >> 3) + 1] = rem.b[2];
		ptr += 1;
	}
	*qnt0 = qnt;
	memset(nc, 0, (int)(len) * sizeof(int));
	return qctr;
}

/*
编码DC系数
*/
void en_DC(int a, int lg, int *m, double delta)
{
	int b = 0, c = 0, i = 0, qf;
	ulong max, min;
	extern int ptr;
	extern uchar *bin;
	extern float ***PTVData;
	extern union Fabs f1;
	int offset0 = 8;
	int offset1 = 32;
	max = 0; min = PTVData[a][b][c];
	for (int j = 0; j < lg; j++) {
		f1.qcf[0][j] = PTVData[a][b][c];
		f1.temp[0][j] &= 0x7fffffff; //将符号位置零求绝对值
		qf = (int)(f1.qcf[0][j] * delta + 0.5); //round的量化
		f1.qcf[0][j] = (float)qf;
		if (qf > max)
			max = qf;
		if (qf < min)
			min = qf;
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

	int nbits, nbitsMax, Nsym, idx;
	int usebit = ptr & 7; //表示当前字节用了已经多少位
	union data {
		uint a;
		uchar b[4];
	} rem;
	nbits = (int)(17.5 + log(delta) / log(2.0)); //这里delta=1/delta，所以减号变加号
	nbitsMax = (int)(log((double)max) / log(2.0)) + 1;
	//编码Maxqsb
	rem.a = max;
	rem.a <<= (32 - usebit - nbits);
	idx = (ptr >> 3);
	bin[idx] |= rem.b[3];
	bin[idx + 1] = rem.b[2];
	bin[idx + 2] = rem.b[1];
	bin[idx + 3] = rem.b[0];
	ptr += nbits;
	rem.a = 0;	//清空寄存器

	//编码Minqsb
	usebit = ptr & 7;
	rem.a = min;
	rem.a <<= (32 - usebit - nbitsMax);
	idx = (ptr >> 3);
	bin[idx] |= rem.b[3];
	bin[idx + 1] = rem.b[2];
	bin[idx + 2] = rem.b[1];
	bin[idx + 3] = rem.b[0];
	ptr += nbitsMax;
	rem.a = 0;

	//编码系数
	SFC se;
	Nsym = max - min + 1;
	for (int ii = 0; ii < lg; ii++) {
		qf = (int)(f1.qcf[0][ii] + 1 - min);
		se = SFcode(qf, Nsym);
		usebit = ptr & 7;
		rem.a = se.code;
		rem.a <<= (32 - usebit - se.lb);
		idx = (ptr >> 3);
		bin[idx] |= rem.b[3];
		bin[idx + 1] = rem.b[2];
		bin[idx + 2] = rem.b[1];
		bin[idx + 3] = rem.b[0];
		ptr += se.lb;
		rem.a = 0;
	}
	return;
}

void en_subDC_noharr(int a, int *w, int *h, int **m, double delta)
{
	extern uchar **sn;
	extern uchar *snbin;
	extern uchar *bin;
	extern int ptr;
	extern int **nc;
	extern int len;
	int max, lg, idx = 0, qctr, maxcf0 = 0, qnt = 0;
	extern union Fabs f1;
	double delta0;
	Uint8_Dat sign;
	delta0 = 1.0 / delta;
	/*PTV第二层*/
	/*1------------------------------------------------------------------*/
	lg = (w[2] * h[2]) << 3;
	max = quanEVEN_DC(a, 0, 0, lg, idx, 2, 0, m[27], delta0*10.0);
	qctr = quant_sub0(sn[0], max, nc[0], lg, idx, delta0, f1.qcf[0], &qnt);
#if ENCODE
	if (qnt) {
		sign.len = handleSn(sn[0], idx, lg, qctr, &maxcf0, f1.qcf[0]);
		sign.dat = snbin;
		en_sub3d_sub2(f1.qcf[0], &sign, lg, maxcf0);
	}
#endif
	/*2------------------------------------------------------------------------*/
	lg = (w[3] * h[2]) << 3;
	max = quanEVEN_DC(a, 0, 2, lg, idx, 2, 0, m[28], delta0*10.0);
	qctr = quant_sub0(sn[0], max, nc[0], lg, idx, delta0, f1.qcf[0], &qnt);
#if ENCODE
	if (qnt) {
		sign.len = handleSn(sn[0], idx, lg, qctr, &maxcf0, f1.qcf[0]);
		sign.dat = snbin;
		en_sub3d_sub2(f1.qcf[0], &sign, lg, maxcf0);
	}
#endif
	/*3----------------------------------------------------------------*/
	lg = (w[2] * h[3]) << 3;
	max = quanEVEN_DC(a, 2, 0, lg, idx, 2, 0, m[29], delta0*10.0);
	qctr = quant_sub0(sn[0], max, nc[0], lg, idx, delta0, f1.qcf[0], &qnt);
#if ENCODE
	if (qnt) {
		sign.len = handleSn(sn[0], idx, lg, qctr, &maxcf0, f1.qcf[0]);
		sign.dat = snbin;
		en_sub3d_sub2(f1.qcf[0], &sign, lg, maxcf0);
	}
#endif
	/*4--------------------------------------------------------*/
	lg = (w[3] * h[3]) << 3;
	max = quanEVEN_DC(a, 2, 2, lg, idx, 2, 0, m[30], delta0*10.0);
	qctr = quant_sub0(sn[0], max, nc[0], lg, idx, delta0, f1.qcf[0], &qnt);
#if ENCODE
	if (qnt) {
		sign.len = handleSn(sn[0], idx, lg, qctr, &maxcf0, f1.qcf[0]);
		sign.dat = snbin;
		en_sub3d_sub2(f1.qcf[0], &sign, lg, maxcf0);
	}
	//FILE *fp = fopen("E:\\程序代码\\量化\\2021-01-15\\test_en_sub3d_sub2_sub\\cbin.txt", "wb");
	//fwrite(bin, sizeof(uchar), (ptr>>3)+1, fp);/*4*/
	//fclose(fp);
#endif
	/*PTV第一层*/
	/*1----------------------------------------------------------------*/
	lg = (w[5] * h[4]) << 3;
	max = quanEVEN_DC(a, 0, 1, lg, idx, 1, 0, m[31], delta0*10.0);
	qctr = quant_sub0(sn[0], max, nc[0], lg, idx, delta0, f1.qcf[0], &qnt);
#if ENCODE
	if (qnt) {
		sign.len = handleSn(sn[0], idx, lg, qctr, &maxcf0, f1.qcf[0]);
		sign.dat = snbin;
		en_sub3d_sub2(f1.qcf[0], &sign, lg, maxcf0);
	}
#endif
	/*2-----------------------------------------------------------------*/
	lg = (w[5] * h[5]) << 3;
	max = quanEVEN_DC(a, 1, 1, lg, idx, 1, 0, m[33], delta0*10.0);
	qctr = quant_sub0(sn[0], max, nc[0], lg, idx, delta0, f1.qcf[0], &qnt);
#if ENCODE
	if (qnt) {
		sign.len = handleSn(sn[0], idx, lg, qctr, &maxcf0, f1.qcf[0]);
		sign.dat = snbin;
		write_en_sub2_sub_data(f1.qcf[0], &sign, lg);
		en_sub3d_sub2(f1.qcf[0], &sign, lg, maxcf0);
	}
#endif

	/*3-------------------------------------------------------------------*/
	lg = (w[4] * h[5]) << 3;
	max = quanEVEN_DC(a, 1, 0, lg, idx, 1, 0, m[32], delta0*10.0);
	qctr = quant_sub0(sn[0], max, nc[0], lg, idx, delta0, f1.qcf[0], &qnt);
#if ENCODE
	if (qnt) {
		sign.len = handleSn(sn[0], idx, lg, qctr, &maxcf0, f1.qcf[0]);
		sign.dat = snbin;
		en_sub3d_sub2(f1.qcf[0], &sign, lg, maxcf0);
	}
#endif

	return;
}

void en_sub7_noharr(int a, int *w, int *h, int **m, double delta)
{
	extern uchar **sn;
	extern uchar *snbin;
	extern uchar *bin;
	extern int ptr;
	extern int **nc;
	extern int len;
	int *max;
	int idx = 0, lg, qctr, maxcf0 = 0, qnt = 0;
	extern union Fabs f1;
	Uint8_Dat sign;
	double delta0;
	max = (int*)calloc(7, sizeof(int));
	if (!max) {
		printf("创建max数组失败！\n");
		exit(1);
	}
	delta0 = 1.0 / delta;
	/*DC后七个系数*/
	lg = (w[6] * h[6]) << 3;
	quanEVEN7(max, a, lg, idx, 0, m[37], delta0*10.0);
	/*1-------------------------------------------------------------------*/
	qctr = quant_sub0(sn[1], max[0], nc[1], lg, idx, delta0, f1.qcf[1], &qnt);
#if ENCODE
	if (qnt) {
		sign.len = handleSn(sn[1], idx, lg, qctr, &maxcf0, f1.qcf[1]);
		sign.dat = snbin;
		en_sub3d_sub2(f1.qcf[1], &sign, lg, maxcf0);
	}
#endif
	/*2------------------------------------------------------------------*/
	qctr = quant_sub0(sn[2], max[1], nc[2], lg, idx, delta0, f1.qcf[2], &qnt);
#if ENCODE
	if (qnt) {
		sign.len = handleSn(sn[2], idx, lg, qctr, &maxcf0, f1.qcf[2]);
		sign.dat = snbin;
		en_sub3d_sub2(f1.qcf[2], &sign, lg, maxcf0);
	}
#endif
	/*3-----------------------------------------------------------------------*/
	qctr = quant_sub0(sn[3], max[2], nc[3], lg, idx, delta0, f1.qcf[3], &qnt);
#if ENCODE
	if (qnt) {
		sign.len = handleSn(sn[3], idx, lg, qctr, &maxcf0, f1.qcf[3]);
		sign.dat = snbin;
		write_en_sub2_sub_data(f1.qcf[3], &sign, lg);
		en_sub3d_sub2(f1.qcf[3], &sign, lg, maxcf0);
	}
#endif
	/*4-------------------------------------------------------------------------*/
	qctr = quant_sub0(sn[4], max[3], nc[4], lg, idx, delta0, f1.qcf[4], &qnt);
#if ENCODE
	if (qnt) {
		sign.len = handleSn(sn[4], idx, lg, qctr, &maxcf0, f1.qcf[4]);
		sign.dat = snbin;
		en_sub3d_sub2(f1.qcf[4], &sign, lg, maxcf0);
	}
#endif
	/*5--------------------------------------------------------------------*/
	qctr = quant_sub0(sn[5], max[4], nc[5], lg, idx, delta0, f1.qcf[5], &qnt);
#if ENCODE
	if (qnt) {
		sign.len = handleSn(sn[5], idx, lg, qctr, &maxcf0, f1.qcf[5]);
		sign.dat = snbin;
		en_sub3d_sub2(f1.qcf[5], &sign, lg, maxcf0);
	}
#endif
	/*6------------------------------------------------------------------------*/
	qctr = quant_sub0(sn[6], max[5], nc[6], lg, idx, delta0, f1.qcf[6], &qnt);
#if ENCODE
	if (qnt) {
		sign.len = handleSn(sn[6], idx, lg, qctr, &maxcf0, f1.qcf[6]);
		sign.dat = snbin;
		en_sub3d_sub2(f1.qcf[6], &sign, lg, maxcf0);
	}
#endif
	/*7--------------------------------------------------------------------*/
	qctr = quant_sub0(sn[7], max[6], nc[7], lg, idx, delta0, f1.qcf[7], &qnt);
#if ENCODE
	if (qnt) {
		sign.len = handleSn(sn[7], idx, lg, qctr, &maxcf0, f1.qcf[7]);
		sign.dat = snbin;
		en_sub3d_sub2(f1.qcf[7], &sign, lg, maxcf0);
	}
#endif

	free(max);
	return;
}

void en_sub8_noharr(int a, int b, int *w, int *h, int **m, double delta)
{
	extern uchar **sn;
	extern uchar *snbin;
	extern uchar *bin;
	extern int ptr;
	extern int **nc;
	extern int len;
	int *max;
	int lg, idx = 0, qctr, maxcf0 = 0, qnt = 0;
	extern union Fabs f1;
	Uint8_Dat sign;
	double delta0;
	max = (int*)calloc(8, sizeof(int));
	if (!max) {
		printf("创建max数组失败！\n");
		exit(1);
	}
	delta0 = 1.0 / delta;

	lg = (w[6] * h[6]) << 3;
	quanEVEN8(max, a, b, lg, idx, 0, m[37], delta0*10.0);
	/*1-------------------------------------------------------*/
	qctr = quant_sub0(sn[0], max[0], nc[0], lg, idx, delta0, f1.qcf[0], &qnt);
#if ENCODE
	if (qnt) {
		sign.len = handleSn(sn[0], idx, lg, qctr, &maxcf0, f1.qcf[0]);
		sign.dat = snbin;
		en_sub3d_sub2(f1.qcf[0], &sign, lg, maxcf0);
	}
#endif

	/*2-----------------------------------------------------------------------*/
	qctr = quant_sub0(sn[1], max[1], nc[1], lg, idx, delta0, f1.qcf[1], &qnt);
#if ENCODE
	if (qnt) {
		sign.len = handleSn(sn[1], idx, lg, qctr, &maxcf0, f1.qcf[1]);
		sign.dat = snbin;
		en_sub3d_sub2(f1.qcf[1], &sign, lg, maxcf0);
	}
#endif

	/*3------------------------------------------------------------------------*/
	qctr = quant_sub0(sn[2], max[2], nc[2], lg, idx, delta0, f1.qcf[2], &qnt);
#if ENCODE
	if (qnt) {
		sign.len = handleSn(sn[2], idx, lg, qctr, &maxcf0, f1.qcf[2]);
		sign.dat = snbin;
		en_sub3d_sub2(f1.qcf[2], &sign, lg, maxcf0);
	}
#endif

	/*4---------------------------------------------------------------------*/
	qctr = quant_sub0(sn[3], max[3], nc[3], lg, idx, delta0, f1.qcf[3], &qnt);
#if ENCODE
	if (qnt) {
		sign.len = handleSn(sn[3], idx, lg, qctr, &maxcf0, f1.qcf[3]);
		sign.dat = snbin;
		en_sub3d_sub2(f1.qcf[3], &sign, lg, maxcf0);
	}
#endif

	/*5--------------------------------------------------------------------*/
	qctr = quant_sub0(sn[4], max[4], nc[4], lg, idx, delta0, f1.qcf[4], &qnt);
#if ENCODE
	if (qnt) {
		sign.len = handleSn(sn[4], idx, lg, qctr, &maxcf0, f1.qcf[4]);
		sign.dat = snbin;
		en_sub3d_sub2(f1.qcf[4], &sign, lg, maxcf0);
	}
#endif

	/*6-------------------------------------------------------------------------*/
	qctr = quant_sub0(sn[5], max[5], nc[5], lg, idx, delta0, f1.qcf[5], &qnt);
#if ENCODE
	if (qnt) {
		sign.len = handleSn(sn[5], idx, lg, qctr, &maxcf0, f1.qcf[5]);
		sign.dat = snbin;
		en_sub3d_sub2(f1.qcf[5], &sign, lg, maxcf0);
	}
#endif

	/*7----------------------------------------------------------------------*/
	qctr = quant_sub0(sn[6], max[6], nc[6], lg, idx, delta0, f1.qcf[6], &qnt);
#if ENCODE
	if (qnt) {
		sign.len = handleSn(sn[6], idx, lg, qctr, &maxcf0, f1.qcf[6]);
		sign.dat = snbin;
		en_sub3d_sub2(f1.qcf[6], &sign, lg, maxcf0);
	}
#endif

	/*8----------------------------------------------------------------------*/
	qctr = quant_sub0(sn[7], max[7], nc[7], lg, idx, delta0, f1.qcf[7], &qnt);
#if ENCODE
	if (qnt) {
		sign.len = handleSn(sn[7], idx, lg, qctr, &maxcf0, f1.qcf[7]);
		sign.dat = snbin;
		en_sub3d_sub2(f1.qcf[7], &sign, lg, maxcf0);
	}
#endif

	free(max);

	return;
}

void en_coef3d_dc(int a, int *w, int *h, int **m, double delta)
{
	extern uchar **sn;
	extern uchar *snbin;
	extern uchar *bin;
	extern int ptr;
	extern int **nc;
	extern int len;
	extern union Fabs f1;
	int idx = 0, max, lg, qctr, maxcf0 = 0, qnt = 0;
	double delta0;
	Uint8_Dat sign;
	delta0 = 1.0 / delta;

	/*DC部分*/
	//1-------------------------------------------------------------------
	idx = 0;
	lg = (w[2] * h[2]) << 2;	//第二次PTV
	max = quanEVEN_DC(a + 4, 0, 0, lg, idx, 2, 1, m[20], delta0*10.0);
	qctr = quant_sub0(sn[0], max, nc[0], lg, idx, delta0, f1.qcf[0], &qnt);
#if ENCODE
	if (qnt) {
		sign.len = handleSn(sn[0], idx, lg, qctr, &maxcf0, f1.qcf[0]);
		sign.dat = snbin;
		en_sub3d_sub2(f1.qcf[0], &sign, lg, maxcf0);
	}
#endif

	//2-------------------------------------------------------------------
	idx = 0;
	lg = (w[3] * h[2]) << 2;
	max = quanEVEN_DC(a + 4, 0, 2, lg, idx, 2, 1, m[21], delta0*10.0);
	qctr = quant_sub0(sn[0], max, nc[0], lg, idx, delta0, f1.qcf[0], &qnt);
#if ENCODE
	if (qnt) {
		sign.len = handleSn(sn[0], idx, lg, qctr, &maxcf0, f1.qcf[0]);
		sign.dat = snbin;
		en_sub3d_sub2(f1.qcf[0], &sign, lg, maxcf0);
	}
#endif

	//3-------------------------------------------------------------------
	idx = 0;
	lg = (w[2] * h[3]) << 2;
	max = quanEVEN_DC(a + 4, 2, 0, lg, idx, 2, 1, m[22], delta0*10.0);
	qctr = quant_sub0(sn[0], max, nc[0], lg, idx, delta0, f1.qcf[0], &qnt);
#if ENCODE
	if (qnt) {
		sign.len = handleSn(sn[0], idx, lg, qctr, &maxcf0, f1.qcf[0]);
		sign.dat = snbin;
		en_sub3d_sub2(f1.qcf[0], &sign, lg, maxcf0);
	}
#endif

	//4-------------------------------------------------------------------
	idx = 0;
	lg = (w[3] * h[3]) << 2;
	max = quanEVEN_DC(a + 4, 2, 2, lg, idx, 2, 1, m[23], delta0*10.0);
	qctr = quant_sub0(sn[0], max, nc[0], lg, idx, delta0, f1.qcf[0], &qnt);
#if ENCODE
	if (qnt) {
		sign.len = handleSn(sn[0], idx, lg, qctr, &maxcf0, f1.qcf[0]);
		sign.dat = snbin;
		en_sub3d_sub2(f1.qcf[0], &sign, lg, maxcf0);
	}
#endif

	//5-------------------------------------------------------------------
	idx = 0;
	lg = (w[5] * h[4]) << 2;	//第一次PTV
	max = quanEVEN_DC(a + 4, 0, 1, lg, idx, 1, 1, m[24], delta0*10.0);
	qctr = quant_sub0(sn[0], max, nc[0], lg, idx, delta0, f1.qcf[0], &qnt);
#if ENCODE
	if (qnt) {
		sign.len = handleSn(sn[0], idx, lg, qctr, &maxcf0, f1.qcf[0]);
		sign.dat = snbin;
		en_sub3d_sub2(f1.qcf[0], &sign, lg, maxcf0);
	}
#endif

	//6-------------------------------------------------------------------
	idx = 0;
	lg = (w[5] * h[5]) << 2;
	max = quanEVEN_DC(a + 4, 1, 1, lg, idx, 1, 1, m[26], delta0*10.0);
	qctr = quant_sub0(sn[0], max, nc[0], lg, idx, delta0, f1.qcf[0], &qnt);
#if ENCODE
	if (qnt) {
		sign.len = handleSn(sn[0], idx, lg, qctr, &maxcf0, f1.qcf[0]);
		sign.dat = snbin;
		en_sub3d_sub2(f1.qcf[0], &sign, lg, maxcf0);
	}
#endif
	//7-------------------------------------------------------------------
	idx = 0;
	lg = (w[4] * h[5]) << 2;
	max = quanEVEN_DC(a + 4, 1, 0, lg, idx, 1, 1, m[25], delta0*10.0);
	qctr = quant_sub0(sn[0], max, nc[0], lg, idx, delta0, f1.qcf[0], &qnt);
#if ENCODE
	if (qnt) {
		sign.len = handleSn(sn[0], idx, lg, qctr, &maxcf0, f1.qcf[0]);
		sign.dat = snbin;
		//write_en_sub2_sub_data(f1.qcf[0], &sign, lg);
		en_sub3d_sub2(f1.qcf[0], &sign, lg, maxcf0);
	}
#endif

	return;
}

void en_coed3d_7(int a, int *w, int *h, int **m, double delta)
{
	extern uchar **sn;
	extern uchar *snbin;
	extern uchar *bin;
	extern int ptr;
	extern int **nc;
	extern int len;
	extern union Fabs f1;
	int *max;
	int idx = 0, lg, qctr, maxcf0 = 0, qnt = 0;
	Uint8_Dat sign;
	double delta0;
	max = (int*)calloc(7, sizeof(int));
	if (!max) {
		printf("创建max数组失败！\n");
		exit(1);
	}
	delta0 = 1.0 / delta;

	/*%level 3: == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == ==*/
	lg = (w[6] * h[6]) << 2;
	quanEVEN7(max, a + 4, lg, idx, 1, m[36], delta0 * 10);
	/*1-------------------------------------------------------------------*/
	qctr = quant_sub0(sn[1], max[0], nc[1], lg, idx, delta0, f1.qcf[1], &qnt);
#if ENCODE
	if (qnt) {
		sign.len = handleSn(sn[1], idx, lg, qctr, &maxcf0, f1.qcf[1]);
		sign.dat = snbin;
		en_sub3d_sub2(f1.qcf[1], &sign, lg, maxcf0);
	}
#endif

	/*2---------------------------------------------------------------------*/
	qctr = quant_sub0(sn[2], max[1], nc[2], lg, idx, delta0, f1.qcf[2], &qnt);
#if ENCODE
	if (qnt) {
		sign.len = handleSn(sn[2], idx, lg, qctr, &maxcf0, f1.qcf[2]);
		sign.dat = snbin;
		en_sub3d_sub2(f1.qcf[2], &sign, lg, maxcf0);
	}
#endif

	/*3----------------------------------------------------------------------*/
	qctr = quant_sub0(sn[3], max[2], nc[3], lg, idx, delta0, f1.qcf[3], &qnt);
#if ENCODE
	if (qnt) {
		sign.len = handleSn(sn[3], idx, lg, qctr, &maxcf0, f1.qcf[3]);
		sign.dat = snbin;
		en_sub3d_sub2(f1.qcf[3], &sign, lg, maxcf0);
	}
#endif

	/*4--------------------------------------------------------------------*/
	qctr = quant_sub0(sn[4], max[3], nc[4], lg, idx, delta0, f1.qcf[4], &qnt);
#if ENCODE
	if (qnt) {
		sign.len = handleSn(sn[4], idx, lg, qctr, &maxcf0, f1.qcf[4]);
		sign.dat = snbin;
		en_sub3d_sub2(f1.qcf[4], &sign, lg, maxcf0);
	}
#endif

	/*5----------------------------------------------------------------------*/
	qctr = quant_sub0(sn[5], max[4], nc[5], lg, idx, delta0, f1.qcf[5], &qnt);
#if ENCODE
	if (qnt) {
		sign.len = handleSn(sn[5], idx, lg, qctr, &maxcf0, f1.qcf[5]);
		sign.dat = snbin;
		en_sub3d_sub2(f1.qcf[5], &sign, lg, maxcf0);
	}
#endif

	/*6-------------------------------------------------------------------*/
	qctr = quant_sub0(sn[6], max[5], nc[6], lg, idx, delta0, f1.qcf[6], &qnt);
#if ENCODE
	if (qnt) {
		sign.len = handleSn(sn[6], idx, lg, qctr, &maxcf0, f1.qcf[6]);
		sign.dat = snbin;
		en_sub3d_sub2(f1.qcf[6], &sign, lg, maxcf0);
	}
#endif

	/*7-------------------------------------------------------------------*/
	qctr = quant_sub0(sn[7], max[6], nc[7], lg, idx, delta0, f1.qcf[7], &qnt);
#if ENCODE
	if (qnt) {
		sign.len = handleSn(sn[7], idx, lg, qctr, &maxcf0, f1.qcf[7]);
		sign.dat = snbin;
		en_sub3d_sub2(f1.qcf[7], &sign, lg, maxcf0);
	}
#endif

	free(max);
	return;
}

void en_coed3d_8(int a, int b, int *w, int *h, int **m, double delta)
{
	extern uchar **sn;
	extern uchar *snbin;
	extern uchar *bin;
	extern int ptr;
	extern int **nc;
	extern int len;
	int *max;
	int lg, idx = 0, qctr, maxcf0 = 0, qnt = 0;
	extern union Fabs f1;
	Uint8_Dat sign;
	double delta0;
	max = (int*)calloc(8, sizeof(int));
	if (!max) {
		printf("创建max数组失败！\n");
		exit(1);
	}
	delta0 = 1.0 / delta;
	lg = (w[6] * h[6]) << 2;

	/*第一次harr变换后的高频部分*/
	quanEVEN8(max, a + 4, b, lg, idx, 1, m[36], delta0*10.0);
	/*1----------------------------------------------------------------*/
	qctr = quant_sub0(sn[0], max[0], nc[0], lg, idx, delta0, f1.qcf[0], &qnt);
#if ENCODE
	if (qnt) {
		sign.len = handleSn(sn[0], idx, lg, qctr, &maxcf0, f1.qcf[0]);
		sign.dat = snbin;
		//write_en_sub2_sub_data(f1.qcf[0], &sign, lg);
		en_sub3d_sub2(f1.qcf[0], &sign, lg, maxcf0);
	}
#endif

	/*2--------------------------------------------------------------------*/
	qctr = quant_sub0(sn[1], max[1], nc[1], lg, idx, delta0, f1.qcf[1], &qnt);
#if ENCODE
	if (qnt) {
		sign.len = handleSn(sn[1], idx, lg, qctr, &maxcf0, f1.qcf[1]);
		sign.dat = snbin;
		en_sub3d_sub2(f1.qcf[1], &sign, lg, maxcf0);
	}
#endif

	/*3--------------------------------------------------------------------*/
	qctr = quant_sub0(sn[2], max[2], nc[2], lg, idx, delta0, f1.qcf[2], &qnt);
#if ENCODE
	if (qnt) {
		sign.len = handleSn(sn[2], idx, lg, qctr, &maxcf0, f1.qcf[2]);
		sign.dat = snbin;
		en_sub3d_sub2(f1.qcf[2], &sign, lg, maxcf0);
	}
#endif

	/*4------------------------------------------------------------------*/
	qctr = quant_sub0(sn[3], max[3], nc[3], lg, idx, delta0, f1.qcf[3], &qnt);
#if ENCODE
	if (qnt) {
		sign.len = handleSn(sn[3], idx, lg, qctr, &maxcf0, f1.qcf[3]);
		sign.dat = snbin;
		en_sub3d_sub2(f1.qcf[3], &sign, lg, maxcf0);
	}
#endif

	/*5-------------------------------------------------------------------*/
	qctr = quant_sub0(sn[4], max[4], nc[4], lg, idx, delta0, f1.qcf[4], &qnt);
#if ENCODE
	if (qnt) {
		sign.len = handleSn(sn[4], idx, lg, qctr, &maxcf0, f1.qcf[4]);
		sign.dat = snbin;
		en_sub3d_sub2(f1.qcf[4], &sign, lg, maxcf0);
	}
#endif

	/*6-----------------------------------------------------------------*/
	qctr = quant_sub0(sn[5], max[5], nc[5], lg, idx, delta0, f1.qcf[5], &qnt);
#if ENCODE
	if (qnt) {
		sign.len = handleSn(sn[5], idx, lg, qctr, &maxcf0, f1.qcf[5]);
		sign.dat = snbin;
		en_sub3d_sub2(f1.qcf[5], &sign, lg, maxcf0);
	}
#endif

	/*7-----------------------------------------------------------------*/
	qctr = quant_sub0(sn[6], max[6], nc[6], lg, idx, delta0, f1.qcf[6], &qnt);
#if ENCODE
	if (qnt) {
		sign.len = handleSn(sn[6], idx, lg, qctr, &maxcf0, f1.qcf[6]);
		sign.dat = snbin;
		en_sub3d_sub2(f1.qcf[6], &sign, lg, maxcf0);
	}
#endif

	/*8------------------------------------------------------------------*/
	qctr = quant_sub0(sn[7], max[7], nc[7], lg, idx, delta0, f1.qcf[7], &qnt);
#if ENCODE
	if (qnt) {
		sign.len = handleSn(sn[7], idx, lg, qctr, &maxcf0, f1.qcf[7]);
		sign.dat = snbin;
		en_sub3d_sub2(f1.qcf[7], &sign, lg, maxcf0);
	}
#endif

	free(max);

	return;
}

void en_coef3dB_dc(int a, int *w, int *h, int **m, double delta)
{
	extern uchar **sn;
	extern uchar *snbin;
	extern uchar *bin;
	extern int ptr;
	extern int **nc;
	extern int len;
	extern union Fabs f1;
	int idx = 0, max, lg, qctr, maxcf0 = 0, qnt = 0;
	double delta0;
	Uint8_Dat sign;
	delta0 = 1.0 / delta;

	/*第二次harr变换后的高频部分*/
	//1-------------------------------------------------------------------
	lg = (w[0] * h[0]) << 1;//第三次PTV
	max = quanEVEN_DC(a + 8, 0, 0, lg, idx, 3, 2, m[10], delta0*10.0);
	qctr = quant_sub0(sn[0], max, nc[0], lg, idx, delta0, f1.qcf[0], &qnt);
#if ENCODE
	if (qnt) {
		sign.len = handleSn(sn[0], idx, lg, qctr, &maxcf0, f1.qcf[0]);
		sign.dat = snbin;
		en_sub3d_sub2(f1.qcf[0], &sign, lg, maxcf0);
	}
#endif

	//2-------------------------------------------------------------------
	lg = (w[1] * h[0]) << 1;
	max = quanEVEN_DC(a + 8, 0, 4, lg, idx, 3, 2, m[11], delta0*10.0);
	qctr = quant_sub0(sn[0], max, nc[0], lg, idx, delta0, f1.qcf[0], &qnt);
#if ENCODE
	if (qnt) {
		sign.len = handleSn(sn[0], idx, lg, qctr, &maxcf0, f1.qcf[0]);
		sign.dat = snbin;
		en_sub3d_sub2(f1.qcf[0], &sign, lg, maxcf0);
	}
#endif

	//3-------------------------------------------------------------------
	lg = (w[1] * h[1]) << 1;
	max = quanEVEN_DC(a + 8, 4, 4, lg, idx, 3, 2, m[13], delta0*10.0);
	qctr = quant_sub0(sn[0], max, nc[0], lg, idx, delta0, f1.qcf[0], &qnt);
#if ENCODE
	if (qnt) {
		sign.len = handleSn(sn[0], idx, lg, qctr, &maxcf0, f1.qcf[0]);
		sign.dat = snbin;
		en_sub3d_sub2(f1.qcf[0], &sign, lg, maxcf0);
	}
#endif
	//4-------------------------------------------------------------------
	lg = (w[0] * h[1]) << 1;
	max = quanEVEN_DC(a + 8, 4, 0, lg, idx, 3, 2, m[12], delta0*10.0);
	qctr = quant_sub0(sn[0], max, nc[0], lg, idx, delta0, f1.qcf[0], &qnt);
#if ENCODE
	if (qnt) {
		sign.len = handleSn(sn[0], idx, lg, qctr, &maxcf0, f1.qcf[0]);
		sign.dat = snbin;
		en_sub3d_sub2(f1.qcf[0], &sign, lg, maxcf0);
	}
#endif

	//5-------------------------------------------------------------------
	lg = (w[3] * h[2]) << 1;//第二次PTV
	max = quanEVEN_DC(a + 8, 0, 2, lg, idx, 2, 2, m[14], delta0*10.0);
	qctr = quant_sub0(sn[0], max, nc[0], lg, idx, delta0, f1.qcf[0], &qnt);
#if ENCODE
	if (qnt) {
		sign.len = handleSn(sn[0], idx, lg, qctr, &maxcf0, f1.qcf[0]);
		sign.dat = snbin;
		en_sub3d_sub2(f1.qcf[0], &sign, lg, maxcf0);
	}
#endif

	//6-------------------------------------------------------------------
	lg = (w[3] * h[3]) << 1;
	max = quanEVEN_DC(a + 8, 2, 2, lg, idx, 2, 2, m[16], delta0*10.0);
	qctr = quant_sub0(sn[0], max, nc[0], lg, idx, delta0, f1.qcf[0], &qnt);
#if ENCODE
	if (qnt) {
		sign.len = handleSn(sn[0], idx, lg, qctr, &maxcf0, f1.qcf[0]);
		sign.dat = snbin;
		en_sub3d_sub2(f1.qcf[0], &sign, lg, maxcf0);
	}
#endif
	//7-------------------------------------------------------------------
	lg = (w[2] * h[3]) << 1;
	max = quanEVEN_DC(a + 8, 2, 0, lg, idx, 2, 2, m[15], delta0*10.0);
	qctr = quant_sub0(sn[0], max, nc[0], lg, idx, delta0, f1.qcf[0], &qnt);
#if ENCODE
	if (qnt) {
		sign.len = handleSn(sn[0], idx, lg, qctr, &maxcf0, f1.qcf[0]);
		sign.dat = snbin;
		en_sub3d_sub2(f1.qcf[0], &sign, lg, maxcf0);
	}
#endif

	//8-------------------------------------------------------------------
	lg = (w[5] * h[4]) << 1;//第一次PTV
	max = quanEVEN_DC(a + 8, 0, 1, lg, idx, 1, 2, m[17], delta0*10.0);
	qctr = quant_sub0(sn[0], max, nc[0], lg, idx, delta0, f1.qcf[0], &qnt);
#if ENCODE
	if (qnt) {
		sign.len = handleSn(sn[0], idx, lg, qctr, &maxcf0, f1.qcf[0]);
		sign.dat = snbin;
		en_sub3d_sub2(f1.qcf[0], &sign, lg, maxcf0);
	}
#endif

	//9-------------------------------------------------------------------
	lg = (w[5] * h[5]) << 1;
	max = quanEVEN_DC(a + 8, 1, 1, lg, idx, 1, 2, m[19], delta0*10.0);
	qctr = quant_sub0(sn[0], max, nc[0], lg, idx, delta0, f1.qcf[0], &qnt);
#if ENCODE
	if (qnt) {
		sign.len = handleSn(sn[0], idx, lg, qctr, &maxcf0, f1.qcf[0]);
		sign.dat = snbin;
		en_sub3d_sub2(f1.qcf[0], &sign, lg, maxcf0);
	}
#endif

	//10-------------------------------------------------------------------
	lg = (w[4] * h[5]) << 1;
	max = quanEVEN_DC(a + 8, 1, 0, lg, idx, 1, 2, m[18], delta0*10.0);
	qctr = quant_sub0(sn[0], max, nc[0], lg, idx, delta0, f1.qcf[0], &qnt);
#if ENCODE
	if (qnt) {
		sign.len = handleSn(sn[0], idx, lg, qctr, &maxcf0, f1.qcf[0]);
		sign.dat = snbin;
		//write_en_sub2_sub_data(f1.qcf[0], &sign, lg);
		en_sub3d_sub2(f1.qcf[0], &sign, lg, maxcf0);
	}
#endif

	return;
}

void en_coef3dB_L3(int H, int W, int a, int b, int c, int *w, int *h, int **m, double delta)
{
	extern uchar **sn;
	extern uchar *snbin;
	extern uchar *bin;
	extern int ptr;
	extern int **nc;
	extern int len;
	int *max;
	int lg, idx = 0, qctr, maxcf0 = 0, qnt = 0;
	extern union Fabs f1;
	Uint8_Dat sign;
	double delta0;
	max = (int*)calloc(3, sizeof(int));
	if (!max) {
		printf("创建max数组失败！\n");
		exit(1);
	}
	delta0 = 1.0 / delta;

	lg = (w[6] * h[6]) << 1;
	quant3(max, H, W, a + 8, b, c, lg, 2, 0, m[35], 0, delta0 * 10);
	/*1-------------------------------------------------------------------------------*/
	qctr = quant_sub0(sn[1], max[0], nc[1], lg, idx, delta0, f1.qcf[1], &qnt); max[0] = 0;
#if ENCODE
	if (qnt) {
		sign.len = handleSn(sn[1], idx, lg, qctr, &maxcf0, f1.qcf[1]);
		sign.dat = snbin;
		en_sub3d_sub2(f1.qcf[1], &sign, lg, maxcf0);
	}
#endif

	/*2-------------------------------------------------------------------------------*/
	qctr = quant_sub0(sn[2], max[1], nc[2], lg, idx, delta0, f1.qcf[2], &qnt); max[1] = 0;
#if ENCODE
	if (qnt) {
		sign.len = handleSn(sn[2], idx, lg, qctr, &maxcf0, f1.qcf[2]);
		sign.dat = snbin;
		en_sub3d_sub2(f1.qcf[2], &sign, lg, maxcf0);
	}
#endif

	/*3----------------------------------------------------------------------------*/
	qctr = quant_sub0(sn[3], max[2], nc[3], lg, idx, delta0, f1.qcf[3], &qnt); max[2] = 0;
#if ENCODE
	if (qnt) {
		sign.len = handleSn(sn[3], idx, lg, qctr, &maxcf0, f1.qcf[3]);
		sign.dat = snbin;
		en_sub3d_sub2(f1.qcf[3], &sign, lg, maxcf0);
	}
#endif

	return;
}

void en_coef3dB_L2_3(int H, int W, int a, int b, int c, int *w, int *h, int **m, double delta)
{
	extern uchar **sn;
	extern uchar *snbin;
	extern uchar *bin;
	extern int ptr;
	extern int **nc;
	extern int len;
	int *max;
	int lg, idx = 0, qctr, maxcf0 = 0, qnt = 0;
	extern union Fabs f1;
	Uint8_Dat sign;
	double delta0;
	max = (int*)calloc(3, sizeof(int));
	if (!max) {
		printf("创建max数组失败！\n");
		exit(1);
	}
	delta0 = 1.0 / delta;

	lg = (2 * w[6] * 2 * h[6]) << 1;
	quant_L2_3(max, H, W, a + 8, b, c, lg, 2, m[38], delta0 * 10);
	/*1-------------------------------------------------------------------------------*/
	qctr = quant_sub0(sn[1], max[0], nc[1], lg, idx, delta0, f1.qcf[1], &qnt); max[0] = 0;
#if ENCODE
	if (qnt) {
		sign.len = handleSn(sn[1], idx, lg, qctr, &maxcf0, f1.qcf[1]);
		sign.dat = snbin;
		en_sub3d_sub2(f1.qcf[1], &sign, lg, maxcf0);
	}
#endif

	/*2-------------------------------------------------------------------------------*/
	qctr = quant_sub0(sn[2], max[1], nc[2], lg, idx, delta0, f1.qcf[2], &qnt); max[1] = 0;
#if ENCODE
	if (qnt) {
		sign.len = handleSn(sn[2], idx, lg, qctr, &maxcf0, f1.qcf[2]);
		sign.dat = snbin;
		en_sub3d_sub2(f1.qcf[2], &sign, lg, maxcf0);
	}
#endif

	/*3----------------------------------------------------------------------------*/
	qctr = quant_sub0(sn[3], max[2], nc[3], lg, idx, delta0, f1.qcf[3], &qnt); max[2] = 0;
#if ENCODE
	if (qnt) {
		sign.len = handleSn(sn[3], idx, lg, qctr, &maxcf0, f1.qcf[3]);
		sign.dat = snbin;
		write_en_sub2_sub_data(f1.qcf[3], &sign, lg);
		en_sub3d_sub2(f1.qcf[3], &sign, lg, maxcf0);
	}
#endif

	return;
}

void en_coef3dB_L2_4(int H, int W, int a, int b, int c, int *w, int *h, int **m, double delta)
{
	extern uchar **sn;
	extern uchar *snbin;
	extern uchar *bin;
	extern int ptr;
	extern int **nc;
	extern int len;
	int *max;
	int lg, idx = 0, qctr, maxcf0 = 0, qnt = 0;
	extern union Fabs f1;
	Uint8_Dat sign;
	double delta0;
	max = (int*)calloc(4, sizeof(int));
	if (!max) {
		printf("创建max数组失败！\n");
		exit(1);
	}
	delta0 = 1.0 / delta;

	lg = (2 * w[6] * 2 * h[6]) << 1;
	quant_L2_4(max, H, W, a + 8, b, c, lg, 2, m[38], delta0 * 10);

	/*1-------------------------------------------------------------------------------*/
	qctr = quant_sub0(sn[1], max[0], nc[1], lg, idx, delta0, f1.qcf[1], &qnt); max[0] = 0;
#if ENCODE
	if (qnt) {
		sign.len = handleSn(sn[1], idx, lg, qctr, &maxcf0, f1.qcf[1]);
		sign.dat = snbin;
		en_sub3d_sub2(f1.qcf[1], &sign, lg, maxcf0);
	}
#endif
	/*2-------------------------------------------------------------------------------*/
	qctr = quant_sub0(sn[2], max[1], nc[2], lg, idx, delta0, f1.qcf[2], &qnt); max[1] = 0;
#if ENCODE
	if (qnt) {
		sign.len = handleSn(sn[2], idx, lg, qctr, &maxcf0, f1.qcf[2]);
		sign.dat = snbin;
		en_sub3d_sub2(f1.qcf[2], &sign, lg, maxcf0);
	}
#endif
	/*3----------------------------------------------------------------------------*/
	qctr = quant_sub0(sn[3], max[2], nc[3], lg, idx, delta0, f1.qcf[3], &qnt); max[2] = 0;
#if ENCODE
	if (qnt) {
		sign.len = handleSn(sn[3], idx, lg, qctr, &maxcf0, f1.qcf[3]);
		sign.dat = snbin;
		en_sub3d_sub2(f1.qcf[3], &sign, lg, maxcf0);
	}
#endif

	/*4----------------------------------------------------------------------------*/
	qctr = quant_sub0(sn[4], max[3], nc[4], lg, idx, delta0, f1.qcf[4], &qnt); max[3] = 0;
#if ENCODE
	if (qnt) {
		sign.len = handleSn(sn[4], idx, lg, qctr, &maxcf0, f1.qcf[4]);
		sign.dat = snbin;
		en_sub3d_sub2(f1.qcf[4], &sign, lg, maxcf0);
	}
#endif

	free(max);
	return;
}

void en_coef3d_dc5B_dc(int a, int *w, int *h, int **m, double delta, int AC)
{
	extern uchar **sn;
	extern uchar *snbin;
	extern uchar *bin;
	extern int ptr;
	extern int **nc;
	extern int len;
	extern union Fabs f1;
	int idx = 0, max, lg, qctr, maxcf0 = 0, qnt = 0;
	double delta0;
	Uint8_Dat sign;
	delta0 = 1.0 / delta;

	if (AC) {
		/*第三次harr变换后的高频部分*/
		lg = w[0] * h[0];//第三次PTV
		max = quanEVEN_DC(a, 0, 0, lg, idx, 3, 3, m[0], delta0*10.0);
		qctr = quant_sub0(sn[0], max, nc[0], lg, idx, delta0, f1.qcf[0], &qnt);
#if ENCODE
		if (qnt) {
			sign.len = handleSn(sn[0], idx, lg, qctr, &maxcf0, f1.qcf[0]);
			sign.dat = snbin;
			en_sub3d_sub2(f1.qcf[0], &sign, lg, maxcf0);
		}
#endif
	}
	else {
		lg = w[0] * h[0];//第三次PTV，对应matlab里Ldc=coef(:,:,2)
		en_DC(a, lg, m[0], delta0);
	}
	//2-------------------------------------------------------------------
	lg = w[1] * h[0];
	max = quanEVEN_DC(a, 0, 4, lg, idx, 3, 3, m[1], delta0*10.0);
	qctr = quant_sub0(sn[0], max, nc[0], lg, idx, delta0, f1.qcf[0], &qnt);
#if ENCODE
	if (qnt) {
		sign.len = handleSn(sn[0], idx, lg, qctr, &maxcf0, f1.qcf[0]);
		sign.dat = snbin;
		en_sub3d_sub2(f1.qcf[0], &sign, lg, maxcf0);
	}
#endif
	//3-------------------------------------------------------------------
	lg = w[1] * h[1];
	max = quanEVEN_DC(a, 4, 4, lg, idx, 3, 3, m[3], delta0*10.0);
	qctr = quant_sub0(sn[0], max, nc[0], lg, idx, delta0, f1.qcf[0], &qnt);
#if ENCODE
	if (qnt) {
		sign.len = handleSn(sn[0], idx, lg, qctr, &maxcf0, f1.qcf[0]);
		sign.dat = snbin;
		en_sub3d_sub2(f1.qcf[0], &sign, lg, maxcf0);
	}
#endif
	//4---------------------------------------------------------------------
	lg = w[0] * h[1];
	max = quanEVEN_DC(a, 4, 0, lg, idx, 3, 3, m[2], delta0*10.0);
	qctr = quant_sub0(sn[0], max, nc[0], lg, idx, delta0, f1.qcf[0], &qnt);
#if ENCODE
	if (qnt) {
		sign.len = handleSn(sn[0], idx, lg, qctr, &maxcf0, f1.qcf[0]);
		sign.dat = snbin;
		en_sub3d_sub2(f1.qcf[0], &sign, lg, maxcf0);
	}
#endif
	//5--------------------------------------------------------------------
	lg = w[3] * h[2];//第二次PTV
	max = quanEVEN_DC(a, 0, 2, lg, idx, 2, 3, m[4], delta0*10.0);
	qctr = quant_sub0(sn[0], max, nc[0], lg, idx, delta0, f1.qcf[0], &qnt);
#if ENCODE
	if (qnt) {
		sign.len = handleSn(sn[0], idx, lg, qctr, &maxcf0, f1.qcf[0]);
		sign.dat = snbin;
		en_sub3d_sub2(f1.qcf[0], &sign, lg, maxcf0);
	}
#endif
	//6-------------------------------------------------------------------
	lg = w[3] * h[3];
	max = quanEVEN_DC(a, 2, 2, lg, idx, 2, 3, m[6], delta0*10.0);
	qctr = quant_sub0(sn[0], max, nc[0], lg, idx, delta0, f1.qcf[0], &qnt);
#if ENCODE
	if (qnt) {
		sign.len = handleSn(sn[0], idx, lg, qctr, &maxcf0, f1.qcf[0]);
		sign.dat = snbin;
		en_sub3d_sub2(f1.qcf[0], &sign, lg, maxcf0);
	}
#endif
	//7-------------------------------------------------------------------
	lg = w[2] * h[3];
	max = quanEVEN_DC(a, 2, 0, lg, idx, 2, 3, m[5], delta0*10.0);
	qctr = quant_sub0(sn[0], max, nc[0], lg, idx, delta0, f1.qcf[0], &qnt);
#if ENCODE
	if (qnt) {
		sign.len = handleSn(sn[0], idx, lg, qctr, &maxcf0, f1.qcf[0]);
		sign.dat = snbin;
		en_sub3d_sub2(f1.qcf[0], &sign, lg, maxcf0);
	}
#endif
	//8-------------------------------------------------------------------
	lg = w[5] * h[4];//第一次PTV
	max = quanEVEN_DC(a, 0, 1, lg, idx, 1, 3, m[7], delta0*10.0);
	qctr = quant_sub0(sn[0], max, nc[0], lg, idx, delta0, f1.qcf[0], &qnt);
#if ENCODE
	if (qnt) {
		sign.len = handleSn(sn[0], idx, lg, qctr, &maxcf0, f1.qcf[0]);
		sign.dat = snbin;
		en_sub3d_sub2(f1.qcf[0], &sign, lg, maxcf0);
	}
#endif
	//9-------------------------------------------------------------------
	lg = w[5] * h[5];
	max = quanEVEN_DC(a, 1, 1, lg, idx, 1, 3, m[9], delta0*10.0);
	qctr = quant_sub0(sn[0], max, nc[0], lg, idx, delta0, f1.qcf[0], &qnt);
#if ENCODE
	if (qnt) {
		sign.len = handleSn(sn[0], idx, lg, qctr, &maxcf0, f1.qcf[0]);
		sign.dat = snbin;
		en_sub3d_sub2(f1.qcf[0], &sign, lg, maxcf0);
	}
#endif
	//10-------------------------------------------------------------------
	lg = w[4] * h[5];
	max = quanEVEN_DC(a, 1, 0, lg, idx, 1, 3, m[8], delta0*10.0);
	qctr = quant_sub0(sn[0], max, nc[0], lg, idx, delta0, f1.qcf[0], &qnt);
#if ENCODE
	if (qnt) {
		sign.len = handleSn(sn[0], idx, lg, qctr, &maxcf0, f1.qcf[0]);
		sign.dat = snbin;
		write_en_sub2_sub_data(f1.qcf[0], &sign, lg);
		en_sub3d_sub2(f1.qcf[0], &sign, lg, maxcf0);
	}
#endif

	return;
}

void en_coef3d_dc5B_L3(int H, int W, int a, int b, int c, int *w, int *h, int **m, double delta)
{
	extern uchar **sn;
	extern uchar *snbin;
	extern uchar *bin;
	extern int ptr;
	extern int **nc;
	extern int len;
	int *max;
	int lg, idx = 0, qctr, maxcf0 = 0, qnt = 0;
	extern union Fabs f1;
	Uint8_Dat sign;
	double delta0;
	max = (int*)calloc(3, sizeof(int));
	if (!max) {
		printf("创建max数组失败！\n");
		exit(1);
	}
	delta0 = 1.0 / delta;

	lg = (w[6] * h[6]);
	quant3(max, H, W, a, b, c, lg, 3, 0, m[34], 0, delta0 * 10);
	/*1-------------------------------------------------------------------------------*/
	qctr = quant_sub0(sn[1], max[0], nc[1], lg, idx, delta0, f1.qcf[1], &qnt); max[0] = 0;
#if ENCODE
	if (qnt) {
		sign.len = handleSn(sn[1], idx, lg, qctr, &maxcf0, f1.qcf[1]);
		sign.dat = snbin;
		en_sub3d_sub2(f1.qcf[1], &sign, lg, maxcf0);
	}
#endif

	/*2-------------------------------------------------------------------------------*/
	qctr = quant_sub0(sn[2], max[1], nc[2], lg, idx, delta0, f1.qcf[2], &qnt); max[1] = 0;
#if ENCODE
	if (qnt) {
		sign.len = handleSn(sn[2], idx, lg, qctr, &maxcf0, f1.qcf[2]);
		sign.dat = snbin;
		en_sub3d_sub2(f1.qcf[2], &sign, lg, maxcf0);
	}
#endif

	/*3----------------------------------------------------------------------------*/
	qctr = quant_sub0(sn[3], max[2], nc[3], lg, idx, delta0, f1.qcf[3], &qnt); max[2] = 0;
#if ENCODE
	if (qnt) {
		sign.len = handleSn(sn[3], idx, lg, qctr, &maxcf0, f1.qcf[3]);
		sign.dat = snbin;
		en_sub3d_sub2(f1.qcf[3], &sign, lg, maxcf0);
	}
#endif

	return;
}

void en_coef3d_dc5B_L2(int H, int W, int a, int b, int c, int *w, int *h, int **m, double delta)
{
	extern uchar **sn;
	extern uchar *snbin;
	extern uchar *bin;
	extern int ptr;
	extern int **nc;
	extern int len;
	int *max;
	int lg, idx = 0, qctr, maxcf0 = 0, qnt = 0;
	extern union Fabs f1;
	Uint8_Dat sign;
	double delta0;
	max = (int*)calloc(3, sizeof(int));
	if (!max) {
		printf("创建max数组失败！\n");
		exit(1);
	}
	delta0 = 1.0 / delta;

	lg = (2 * w[6] * 2 * h[6]);
	quant3(max, H, W, a, b, c, lg, 3, 1, m[39], 1, delta0 * 10);
	/*1-------------------------------------------------------------------------------*/
	qctr = quant_sub0(sn[1], max[0], nc[1], lg, idx, delta0, f1.qcf[1], &qnt); max[0] = 0;
#if ENCODE
	if (qnt) {
		sign.len = handleSn(sn[1], idx, lg, qctr, &maxcf0, f1.qcf[1]);
		sign.dat = snbin;
		en_sub3d_sub2(f1.qcf[1], &sign, lg, maxcf0);
	}
#endif

	/*2-------------------------------------------------------------------------------*/
	qctr = quant_sub0(sn[2], max[1], nc[2], lg, idx, delta0, f1.qcf[2], &qnt); max[1] = 0;
#if ENCODE
	if (qnt) {
		sign.len = handleSn(sn[2], idx, lg, qctr, &maxcf0, f1.qcf[2]);
		sign.dat = snbin;
		en_sub3d_sub2(f1.qcf[2], &sign, lg, maxcf0);
	}
#endif

	/*3----------------------------------------------------------------------------*/
	qctr = quant_sub0(sn[3], max[2], nc[3], lg, idx, delta0, f1.qcf[3], &qnt); max[2] = 0;
#if ENCODE
	if (qnt) {
		sign.len = handleSn(sn[3], idx, lg, qctr, &maxcf0, f1.qcf[3]);
		sign.dat = snbin;
		en_sub3d_sub2(f1.qcf[3], &sign, lg, maxcf0);
	}
#endif

	return;

}

void en_coef3d_dc5B_L1(int H, int W, int a, int b, int c, int *w, int *h, int **m, double delta)
{
	extern uchar **sn;
	extern uchar *snbin;
	extern uchar *bin;
	extern int ptr;
	extern int **nc;
	extern int len;
	int *max;
	int lg, idx = 0, qctr, maxcf0 = 0, qnt = 0;
	extern union Fabs f1;
	Uint8_Dat sign;
	double delta0;
	max = (int*)calloc(3, sizeof(int));
	if (!max) {
		printf("创建max数组失败！\n");
		exit(1);
	}
	delta0 = 1.0 / delta;

	lg = (4 * w[6] * 4 * h[6]);
	quant3(max, H, W, a, b, c, lg, 3, 2, m[40], 1, delta0 * 10);
	/*1-------------------------------------------------------------------------------*/
	qctr = quant_sub0(sn[1], max[0], nc[1], lg, idx, delta0, f1.qcf[1], &qnt); max[0] = 0;
#if ENCODE
	if (qnt) {
		sign.len = handleSn(sn[1], idx, lg, qctr, &maxcf0, f1.qcf[1]);
		sign.dat = snbin;
		en_sub3d_sub2(f1.qcf[1], &sign, lg, maxcf0);
	}
#endif

	/*2-------------------------------------------------------------------------------*/
	qctr = quant_sub0(sn[2], max[1], nc[2], lg, idx, delta0, f1.qcf[2], &qnt); max[1] = 0;
#if ENCODE
	if (qnt) {
		sign.len = handleSn(sn[2], idx, lg, qctr, &maxcf0, f1.qcf[2]);
		sign.dat = snbin;
		en_sub3d_sub2(f1.qcf[2], &sign, lg, maxcf0);
	}
#endif

	/*3----------------------------------------------------------------------------*/
	qctr = quant_sub0(sn[3], max[2], nc[3], lg, idx, delta0, f1.qcf[3], &qnt); max[2] = 0;
#if ENCODE
	if (qnt) {
		sign.len = handleSn(sn[3], idx, lg, qctr, &maxcf0, f1.qcf[3]);
		sign.dat = snbin;
		en_sub3d_sub2(f1.qcf[3], &sign, lg, maxcf0);
	}
#endif

	return;

}




/*输入为W,H; W为宽，H为高，例如1024*1024*96的数据，W=1024，H=1024*
*输出为w,h,m; w和h分别为PTV每个小块的宽和高，m为希尔伯特扫描路径*
*w,h为大小为7的一维数组，m是38行的二维数组，对应38个希尔伯特路径*/
int GetM(int W, int H, int *w, int *h, int **m)
{
	int h1h, w1h, h2h, w2h, h3h, w3h;
	int h3l, w3l, h2l, w2l, h1l, w1l;
	W >>= 3; H >>= 3;
	/*求PTV后每个小块的大小*/
	if (!(W % 2)) {
		w3h = W >> 1; w3l = W >> 1;
	}
	else {
		w3h = (W + 1) >> 1; w3l = w3h - 1;
	}
	if (!(w3h % 2)) {
		w2h = w3h >> 1; w2l = w3h >> 1;
	}
	else {
		w2h = (w3h + 1) >> 1; w2l = w2h - 1;
	}
	if (!(w2h % 2)) {
		w1h = w2h >> 1; w1l = w1h;
	}
	else {
		w1h = (w2h + 1) >> 1; w1l = w1h - 1;
	}
	if (!(H % 2)) {
		h3h = H >> 1; h3l = H >> 1;
	}
	else {
		h3h = (H + 1) >> 1; h3l = h3h - 1;
	}
	if (!(h3h % 2)) {
		h2h = h3h >> 1; h2l = h2h;
	}
	else {
		h2h = (h3h + 1) >> 1; h2l = h2h - 1;
	}
	if (!(h2h % 2)) {
		h1h = h2h >> 1; h1l = h1h;
	}
	else {
		h1h = (h2h + 1) >> 1; h1l = h1h - 1;
	}
	w[0] = w1h; w[1] = w1l;
	w[2] = w2h; w[3] = w2l;
	w[4] = w3h; w[5] = w3l;
	w[6] = W;
	h[0] = h1h; h[1] = h1l;
	h[2] = h2h; h[3] = h2l;
	h[4] = h3h; h[5] = h3l;
	h[6] = H;

	/*求路径M*/
	m[0] = modify_Hilbert2D(w1h, h1h);
	m[1] = modify_Hilbert2D(w1l, h1h);
	m[2] = modify_Hilbert2D(w1h, h1l);
	m[3] = modify_Hilbert2D(w1l, h1l);

	m[4] = modify_Hilbert2D(w2l, h2h);
	m[5] = modify_Hilbert2D(w2h, h2l);
	m[6] = modify_Hilbert2D(w2l, h2l);

	m[7] = modify_Hilbert2D(w3l, h3h);
	m[8] = modify_Hilbert2D(w3h, h3l);
	m[9] = modify_Hilbert2D(w3l, h3l);
	/**/
	m[10] = modify_Hilbert3D(w1h, h1h, 1);
	m[11] = modify_Hilbert3D(w1l, h1h, 1);
	m[12] = modify_Hilbert3D(w1h, h1l, 1);
	m[13] = modify_Hilbert3D(w1l, h1l, 1);

	m[14] = modify_Hilbert3D(w2l, h2h, 1);
	m[15] = modify_Hilbert3D(w2h, h2l, 1);
	m[16] = modify_Hilbert3D(w2l, h2l, 1);

	m[17] = modify_Hilbert3D(w3l, h3h, 1);
	m[18] = modify_Hilbert3D(w3h, h3l, 1);
	m[19] = modify_Hilbert3D(w3l, h3l, 1);
	/**/
	m[20] = modify_Hilbert3D(w2h, h2h, 2);//32*32*4
	m[21] = modify_Hilbert3D(w2l, h2h, 2);
	m[22] = modify_Hilbert3D(w2h, h2l, 2);
	m[23] = modify_Hilbert3D(w2l, h2l, 2);

	m[24] = modify_Hilbert3D(w3l, h3h, 2);//64*64*4
	m[25] = modify_Hilbert3D(w3h, h3l, 2);
	m[26] = modify_Hilbert3D(w3l, h3l, 2);
	/**/
	m[27] = modify_Hilbert3D(w2h, h2h, 3);//32*32*8
	m[28] = modify_Hilbert3D(w2l, h2h, 3);
	m[29] = modify_Hilbert3D(w2h, h2l, 3);
	m[30] = modify_Hilbert3D(w2l, h2l, 3);

	m[31] = modify_Hilbert3D(w3l, h3h, 3);//64*64*8
	m[32] = modify_Hilbert3D(w3h, h3l, 3);
	m[33] = modify_Hilbert3D(w3l, h3l, 3);

	m[34] = modify_Hilbert2D(W, H);//128*128
	m[35] = modify_Hilbert3D(W, H, 1);//128*128*2
	m[36] = modify_Hilbert3D(W, H, 2);//128*128*4
	m[37] = modify_Hilbert3D(W, H, 3);//128*128*8

	m[38] = modify_Hilbert3D(2 * W, 2 * H, 1);
	m[39] = modify_Hilbert2D(2 * W, 2 * H);
	m[40] = modify_Hilbert2D(4 * W, 4 * H);

	return 0;
}

int handleSn(uchar *sgn, int idx, int lg, int qctr, int *maxcf, float *fabsA)
{
	extern unsigned char *snbin;
	unsigned char temp = 0;
	int cnt = 7, num = 0, max = 0;
	int ia = 0;
	int len_snbin;
	if (qctr) { //quantTHD量化的
		for (int i = idx; i < idx + lg; i++) {
			if (fabsA[i] > max)
				max = fabsA[i];
			if (sgn[i]) {
				sgn[i]--;
				temp |= (sgn[i] << cnt);
				cnt--;
				ia++;
			}
			if (cnt == -1) {
				cnt = 7;
				snbin[num] = temp;
				num++;
				temp = 0;
			}
		}
	}
	else { //quantEVEN量化，这里量化后是没有0值的。因为是向上取整
		for (int i = idx; i < idx + lg; i++) {
			fabsA[i]--;
			if (fabsA[i] > max)
				max = fabsA[i];
			if (!sgn[i]) {
				temp <<= 1;
				cnt--;
			}
			else {
				sgn[i]--;
				temp |= (sgn[i] << cnt);
				cnt--;
			}
			if (cnt == -1) {
				cnt = 7;
				snbin[num] = temp;
				num++;
				temp = 0;
			}
		}
	}
	if (num) {
		if (cnt == 7) {
			len_snbin = num * 8;
		}
		else {
			snbin[num] = temp; num++;
			len_snbin = num * 8 - cnt - 1;
		}
	}
	else {
		len_snbin = 7 - cnt;
		snbin[num] = temp;
	}
	*maxcf = max;
	return len_snbin; //此处的长度是snbin的总位数，不是数组的长度。
}