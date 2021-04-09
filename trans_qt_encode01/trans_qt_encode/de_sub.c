#define _CRT_SECURE_NO_WARNINGS
#include "de_sub.h"
#include "modify_Hilbert3D.h"
#include "modify_Hilbert2D.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "all.h"
#include "encoding.h"
#include "quant_inv.h"

float* de_DCinv(char* bin, int a, int lg, int* m, double delta)
{
	union data {
		unsigned long long a;
		uchar b[8];
	} rem1, rem2;
	rem1.a = 0; rem2.a = 0;
	extern int ptr;
	extern float ***rePTVData;
	float *deqcf = (float*)calloc(lg, sizeof(float));
	int nbits, maxqsb, minqsb, nbitsMax, Nsym, x;
	int b = 0, c = 0;
	int idx = ptr >> 3;
	x = ptr & 7;
	rem1.b[5] = bin[idx];
	rem1.b[4] = bin[idx + 1];
	rem1.b[3] = bin[idx + 2];
	rem1.b[2] = bin[idx + 3];
	rem1.b[1] = bin[idx + 4];
	rem1.b[0] = bin[idx + 5];
	rem1.a <<= x;
	rem1.b[7] = 0;
	rem1.b[6] = 0;

	/*解码maxqsb*/
	nbits = (int)(17.5 - log(delta) / log(2.0));
	rem1.a <<= nbits;
	rem2.b[1] = rem1.b[7];
	rem2.b[0] = rem1.b[6];
	maxqsb = (int)rem2.a;
	rem1.b[7] = 0;
	rem1.b[6] = 0;
	ptr += nbits;

	/*解码minqsb*/
	nbitsMax = (int)(log(maxqsb) / log(2)) + 1;
	rem1.a <<= nbitsMax;
	rem2.b[1] = rem1.b[7];
	rem2.b[0] = rem1.b[6];
	minqsb = (int)rem2.a;
	rem1.b[7] = 0;
	rem1.b[6] = 0;
	ptr += nbitsMax;

	Nsym = maxqsb - minqsb + 1;
	/*解码bin*/
	int j = 0;
	int offset0 = 8;
	int offset1 = 32;
	for (int i = 0; i < lg; i++){
		DES de = deSFcode(bin, Nsym);
		rePTVData[a][b][c] = (float)((de.sym - 1 + minqsb) * delta);
		deqcf[i] = rePTVData[a][b][c];
		switch (m[j])
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
		}j++;
	}

	return deqcf;
}

void de_subDC_noharr(int a, int *w, int *h, int **m, double delta, int lenbits)
{
	extern uchar *bin;
	extern int ptr;
	extern int **nc;
	extern int len;
	int lg, qctr;
	double delta0;

	delta0 = 1.0 / delta;
	DE_S_SUB desub;
	float *qcf0;
	uchar qnt = 0, trim = 0;
	uint ctr1 = 0;

	/*PTV第二层*/
	/*1------------------------------------------------------------------*/
	lg = (w[2] * h[2]) << 3;
	qctr = find_qctr(bin, &qnt, &trim, &ctr1);
	if (qnt){
		desub = de_sub3d_sub2(bin, lg, qctr, lenbits);	//解码
		free(desub.sn.dat); free(desub.cfw.dat);
		qcf0 = rst_sub(desub.cf.dat, delta, lg, qnt, qctr, trim, ctr1);
	}
	else{
		qcf0 = (float*)calloc(lg, sizeof(float));
		printf("子带全为零\n");
	}
	scan_inv(a, 0, 0, lg, m[27], 2, 0, 1, qcf0);	//反向扫描

	/*2------------------------------------------------------------------------*/
	lg = (w[3] * h[2]) << 3;
	qctr = find_qctr(bin, &qnt, &trim, &ctr1);
	if (qnt){
		desub = de_sub3d_sub2(bin, lg, qctr, lenbits);	//解码
		free(desub.sn.dat); free(desub.cfw.dat);
		qcf0 = rst_sub(desub.cf.dat, delta, lg, qnt, qctr, trim, ctr1);
	}
	else{
		qcf0 = (float*)calloc(lg, sizeof(float));
		printf("子带全为零\n");
	}
	scan_inv(a, 0, 2, lg, m[28], 2, 0, 1, qcf0);	//反向扫描

	/*3----------------------------------------------------------------*/
	lg = (w[2] * h[3]) << 3;
	qctr = find_qctr(bin, &qnt, &trim, &ctr1);
	if (qnt){
		desub = de_sub3d_sub2(bin, lg, qctr, lenbits);	//解码
		free(desub.sn.dat); free(desub.cfw.dat);
		qcf0 = rst_sub(desub.cf.dat, delta, lg, qnt, qctr, trim, ctr1);
	}
	else{
		qcf0 = (float*)calloc(lg, sizeof(float));
		printf("子带全为零\n");
	}
	scan_inv(a, 2, 0, lg, m[29], 2, 0, 1, qcf0);	//反向扫描

	/*4--------------------------------------------------------*/
	lg = (w[3] * h[3]) << 3;
	qctr = find_qctr(bin, &qnt, &trim, &ctr1);
	if (qnt){
		desub = de_sub3d_sub2(bin, lg, qctr, lenbits);	//解码
		free(desub.sn.dat); free(desub.cfw.dat);
		qcf0 = rst_sub(desub.cf.dat, delta, lg, qnt, qctr, trim, ctr1);
	}
	else{
		qcf0 = (float*)calloc(lg, sizeof(float));
		printf("子带全为零\n");
	}
	scan_inv(a, 2, 2, lg, m[30], 2, 0, 1, qcf0);	//反向扫描

	/*PTV第一层*/
	/*1----------------------------------------------------------------*/
	lg = (w[5] * h[4]) << 3;
	qctr = find_qctr(bin, &qnt, &trim, &ctr1);
	if (qnt){
		desub = de_sub3d_sub2(bin, lg, qctr, lenbits);	//解码
		free(desub.sn.dat); free(desub.cfw.dat);
		qcf0 = rst_sub(desub.cf.dat, delta, lg, qnt, qctr, trim, ctr1);
	}
	else{
		qcf0 = (float*)calloc(lg, sizeof(float));
		printf("子带全为零\n");
	}

	scan_inv(a, 0, 1, lg, m[31], 1, 0, 1, qcf0);	//反向扫描
	/*3-----------------------------------------------------------------*/
	lg = (w[5] * h[5]) << 3;
	qctr = find_qctr(bin, &qnt, &trim, &ctr1);
	if (qnt){
		desub = de_sub3d_sub2(bin, lg, qctr, lenbits);	//解码
		free(desub.sn.dat); free(desub.cfw.dat);
	qcf0 = rst_sub(desub.cf.dat, delta, lg, qnt, qctr, trim, ctr1);
	}
	else{
		qcf0 = (float*)calloc(lg, sizeof(float));
		printf("子带全为零\n");
	}

	scan_inv(a, 1, 1, lg, m[33], 1, 0, 1, qcf0);	//反向扫描
	/*2-------------------------------------------------------------------*/
	lg = (w[4] * h[5]) << 3;
	qctr = find_qctr(bin, &qnt, &trim, &ctr1);
	if (qnt){
		desub = de_sub3d_sub2(bin, lg, qctr, lenbits);	//解码
		free(desub.sn.dat); free(desub.cfw.dat);
	qcf0 = rst_sub(desub.cf.dat, delta, lg, qnt, qctr, trim, ctr1);
	}
	else{
		qcf0 = (float*)calloc(lg, sizeof(float));
		printf("子带全为零\n");
	}

	scan_inv(a, 1, 0, lg, m[32], 1, 0, 1, qcf0);	//反向扫描


	return;
}

void de_sub7_noharr(int a, int *w, int *h, int **m, double delta, int lenbits)
{
	extern uchar *bin;
	extern int ptr;
	extern int len;
	int lg, qctr;
	double delta0;
	delta0 = 1.0 / delta;
	DE_S_SUB desub;
	uchar qnt = 0, trim = 0;
	uint ctr1 = 0;
	float *qcf0;
	/*DC后七个系数*/
	lg = (w[6] * h[6]) << 3;
	/*1-------------------------------------------------------------------*/
	qctr = find_qctr(bin, &qnt, &trim, &ctr1);
	if (qnt){
		desub = de_sub3d_sub2(bin, lg, qctr, lenbits);	//解码
		free(desub.sn.dat); free(desub.cfw.dat);
	qcf0 = rst_sub(desub.cf.dat, delta, lg, qnt, qctr, trim, ctr1);
	}
	else{
		qcf0 = (float*)calloc(lg, sizeof(float));
		printf("子带全为零\n");
	}

	scan_inv(a, 0, 1, lg, m[37], 3, 0, 0, qcf0);	//反向扫描

	/*2------------------------------------------------------------------*/
	qctr = find_qctr(bin, &qnt, &trim, &ctr1);
	if (qnt){
		desub = de_sub3d_sub2(bin, lg, qctr, lenbits);	//解码
		free(desub.sn.dat); free(desub.cfw.dat);
	qcf0 = rst_sub(desub.cf.dat, delta, lg, qnt, qctr, trim, ctr1);
	}
	else{
		qcf0 = (float*)calloc(lg, sizeof(float));
		printf("子带全为零\n");
	}

	scan_inv(a, 0, 2, lg, m[37], 3, 0, 0, qcf0);	//反向扫描

	/*3-----------------------------------------------------------------------*/
	qctr = find_qctr(bin, &qnt, &trim, &ctr1);
	if (qnt){
		desub = de_sub3d_sub2(bin, lg, qctr, lenbits);	//解码
		free(desub.sn.dat); free(desub.cfw.dat);
	qcf0 = rst_sub(desub.cf.dat, delta, lg, qnt, qctr, trim, ctr1);
	}
	else{
		qcf0 = (float*)calloc(lg, sizeof(float));
		printf("子带全为零\n");
	}

	scan_inv(a, 0, 3, lg, m[37], 3, 0, 0, qcf0);	//反向扫描

	/*4-------------------------------------------------------------------------*/
	qctr = find_qctr(bin, &qnt, &trim, &ctr1);
	if (qnt){
		desub = de_sub3d_sub2(bin, lg, qctr, lenbits);	//解码
		free(desub.sn.dat); free(desub.cfw.dat);
	qcf0 = rst_sub(desub.cf.dat, delta, lg, qnt, qctr, trim, ctr1);
	}
	else{
		qcf0 = (float*)calloc(lg, sizeof(float));
		printf("子带全为零\n");
	}

	scan_inv(a, 0, 4, lg, m[37], 3, 0, 0, qcf0);	//反向扫描

	/*5--------------------------------------------------------------------*/
	qctr = find_qctr(bin, &qnt, &trim, &ctr1);
	if (qnt){
		desub = de_sub3d_sub2(bin, lg, qctr, lenbits);	//解码
		free(desub.sn.dat); free(desub.cfw.dat);
	qcf0 = rst_sub(desub.cf.dat, delta, lg, qnt, qctr, trim, ctr1);
	}
	else{
		qcf0 = (float*)calloc(lg, sizeof(float));
		printf("子带全为零\n");
	}

	scan_inv(a, 0, 5, lg, m[37], 3, 0, 0, qcf0);	//反向扫描

	/*6------------------------------------------------------------------------*/
	qctr = find_qctr(bin, &qnt, &trim, &ctr1);
	if (qnt){
		desub = de_sub3d_sub2(bin, lg, qctr, lenbits);	//解码
		free(desub.sn.dat); free(desub.cfw.dat);
	qcf0 = rst_sub(desub.cf.dat, delta, lg, qnt, qctr, trim, ctr1);
	}
	else{
		qcf0 = (float*)calloc(lg, sizeof(float));
		printf("子带全为零\n");
	}

	scan_inv(a, 0, 6, lg, m[37], 3, 0, 0, qcf0);	//反向扫描

	/*7--------------------------------------------------------------------*/
	qctr = find_qctr(bin, &qnt, &trim, &ctr1);
	if (qnt){
		desub = de_sub3d_sub2(bin, lg, qctr, lenbits);	//解码
		free(desub.sn.dat); free(desub.cfw.dat);
	qcf0 = rst_sub(desub.cf.dat, delta, lg, qnt, qctr, trim, ctr1);
	}
	else{
		qcf0 = (float*)calloc(lg, sizeof(float));
		printf("子带全为零\n");
	}

	scan_inv(a, 0, 7, lg, m[37], 3, 0, 0, qcf0);	//反向扫描

	return;
}

void de_sub8_noharr(int a, int b, int *w, int *h, int **m, double delta, int lenbits)
{
	extern uchar *bin;
	extern int ptr;
	extern int len;
	int lg, qctr;
	double delta0;
	float *qcf0;
	delta0 = 1.0 / delta;
	DE_S_SUB desub;
	uchar qnt = 0, trim = 0;
	uint ctr1 = 0;

	lg = (w[6] * h[6]) << 3;
	/*1-------------------------------------------------------*/
	qctr = find_qctr(bin, &qnt, &trim, &ctr1);
	if (qnt){
		desub = de_sub3d_sub2(bin, lg, qctr, lenbits);	//解码
		free(desub.sn.dat); free(desub.cfw.dat);
	qcf0 = rst_sub(desub.cf.dat, delta, lg, qnt, qctr, trim, ctr1);
	}
	else{
		qcf0 = (float*)calloc(lg, sizeof(float));
		printf("子带全为零\n");
	}

	scan_inv(a, b, 0, lg, m[37], 3, 0, 0, qcf0);	//反向扫描

	/*2-----------------------------------------------------------------------*/
	qctr = find_qctr(bin, &qnt, &trim, &ctr1);
	if (qnt){
		desub = de_sub3d_sub2(bin, lg, qctr, lenbits);	//解码
		free(desub.sn.dat); free(desub.cfw.dat);
	qcf0 = rst_sub(desub.cf.dat, delta, lg, qnt, qctr, trim, ctr1);
	}
	else{
		qcf0 = (float*)calloc(lg, sizeof(float));
		printf("子带全为零\n");
	}

	scan_inv(a, b, 1, lg, m[37], 3, 0, 0, qcf0);	//反向扫描

	/*3------------------------------------------------------------------------*/
	qctr = find_qctr(bin, &qnt, &trim, &ctr1);
	if (qnt){
		desub = de_sub3d_sub2(bin, lg, qctr, lenbits);	//解码
		free(desub.sn.dat); free(desub.cfw.dat);
	qcf0 = rst_sub(desub.cf.dat, delta, lg, qnt, qctr, trim, ctr1);
	}
	else{
		qcf0 = (float*)calloc(lg, sizeof(float));
		printf("子带全为零\n");
	}

	scan_inv(a, b, 2, lg, m[37], 3, 0, 0, qcf0);	//反向扫描

	/*4---------------------------------------------------------------------*/
	qctr = find_qctr(bin, &qnt, &trim, &ctr1);
	if (qnt){
		desub = de_sub3d_sub2(bin, lg, qctr, lenbits);	//解码
		free(desub.sn.dat); free(desub.cfw.dat);
	qcf0 = rst_sub(desub.cf.dat, delta, lg, qnt, qctr, trim, ctr1);
	}
	else{
		qcf0 = (float*)calloc(lg, sizeof(float));
		printf("子带全为零\n");
	}

	scan_inv(a, b, 3, lg, m[37], 3, 0, 0, qcf0);	//反向扫描

	/*5--------------------------------------------------------------------*/
	qctr = find_qctr(bin, &qnt, &trim, &ctr1);
	if (qnt){
		desub = de_sub3d_sub2(bin, lg, qctr, lenbits);	//解码
		free(desub.sn.dat); free(desub.cfw.dat);
	qcf0 = rst_sub(desub.cf.dat, delta, lg, qnt, qctr, trim, ctr1);
	}
	else{
		qcf0 = (float*)calloc(lg, sizeof(float));
		printf("子带全为零\n");
	}

	scan_inv(a, b, 4, lg, m[37], 3, 0, 0, qcf0);	//反向扫描

	/*6-------------------------------------------------------------------------*/
	qctr = find_qctr(bin, &qnt, &trim, &ctr1);
	if (qnt){
		desub = de_sub3d_sub2(bin, lg, qctr, lenbits);	//解码
		free(desub.sn.dat); free(desub.cfw.dat);
	qcf0 = rst_sub(desub.cf.dat, delta, lg, qnt, qctr, trim, ctr1);
	}
	else{
		qcf0 = (float*)calloc(lg, sizeof(float));
		printf("子带全为零\n");
	}

	scan_inv(a, b, 5, lg, m[37], 3, 0, 0, qcf0);	//反向扫描

	/*7----------------------------------------------------------------------*/
	qctr = find_qctr(bin, &qnt, &trim, &ctr1);
	if (qnt){
		desub = de_sub3d_sub2(bin, lg, qctr, lenbits);	//解码
		free(desub.sn.dat); free(desub.cfw.dat);
	qcf0 = rst_sub(desub.cf.dat, delta, lg, qnt, qctr, trim, ctr1);
	}
	else{
		qcf0 = (float*)calloc(lg, sizeof(float));
		printf("子带全为零\n");
	}

	scan_inv(a, b, 6, lg, m[37], 3, 0, 0, qcf0);	//反向扫描

	/*8----------------------------------------------------------------------*/
	qctr = find_qctr(bin, &qnt, &trim, &ctr1);
	if (qnt){
		desub = de_sub3d_sub2(bin, lg, qctr, lenbits);	//解码
		free(desub.sn.dat); free(desub.cfw.dat);
	qcf0 = rst_sub(desub.cf.dat, delta, lg, qnt, qctr, trim, ctr1);
	}
	else{
		qcf0 = (float*)calloc(lg, sizeof(float));
		printf("子带全为零\n");
	}

	scan_inv(a, b, 7, lg, m[37], 3, 0, 0, qcf0);	//反向扫描

	return;
}

void de_coef3d_dc(int a, int *w, int *h, int **m, double delta, int lenbits)
{
	extern uchar *bin;
	extern int ptr;
	extern int **nc;
	extern int len;
	int lg, qctr;
	double delta0;

	delta0 = 1.0 / delta;
	DE_S_SUB desub;
	float *qcf0;
	uchar qnt = 0, trim = 0;
	uint ctr1 = 0;

	/*DC部分*/
	//1-------------------------------------------------------------------
	lg = (w[2] * h[2]) << 2;	//第二次PTV
	qctr = find_qctr(bin, &qnt, &trim, &ctr1);
	if (qnt){
		desub = de_sub3d_sub2(bin, lg, qctr, lenbits);	//解码
		free(desub.sn.dat); free(desub.cfw.dat);
	qcf0 = rst_sub(desub.cf.dat, delta, lg, qnt, qctr, trim, ctr1);	//反量化
	}
	else{
		qcf0 = (float*)calloc(lg, sizeof(float));
		printf("子带全为零\n");
	}

	scan_inv(a + 4, 0, 0, lg, m[20], 2, 1, 1, qcf0);	//反向扫描

	//2-------------------------------------------------------------------
	lg = (w[3] * h[2]) << 2;
	qctr = find_qctr(bin, &qnt, &trim, &ctr1);
	if (qnt){
		desub = de_sub3d_sub2(bin, lg, qctr, lenbits);	//解码
		free(desub.sn.dat); free(desub.cfw.dat);
	qcf0 = rst_sub(desub.cf.dat, delta, lg, qnt, qctr, trim, ctr1);	//反量化
	}
	else{
		qcf0 = (float*)calloc(lg, sizeof(float));
		printf("子带全为零\n");
	}

	scan_inv(a + 4, 0, 2, lg, m[21], 2, 1, 1, qcf0);	//反向扫描

	//3-------------------------------------------------------------------
	lg = (w[2] * h[3]) << 2;
	qctr = find_qctr(bin, &qnt, &trim, &ctr1);
	if (qnt){
		desub = de_sub3d_sub2(bin, lg, qctr, lenbits);	//解码
		free(desub.sn.dat); free(desub.cfw.dat);
	qcf0 = rst_sub(desub.cf.dat, delta, lg, qnt, qctr, trim, ctr1);	//反量化
	}
	else{
		qcf0 = (float*)calloc(lg, sizeof(float));
		printf("子带全为零\n");
	}

	scan_inv(a + 4, 2, 0, lg, m[22], 2, 1, 1, qcf0);	//反向扫描

	//4-------------------------------------------------------------------
	lg = (w[3] * h[3]) << 2;
	qctr = find_qctr(bin, &qnt, &trim, &ctr1);
	if (qnt){
		desub = de_sub3d_sub2(bin, lg, qctr, lenbits);	//解码
		free(desub.sn.dat); free(desub.cfw.dat);
	qcf0 = rst_sub(desub.cf.dat, delta, lg, qnt, qctr, trim, ctr1);	//反量化
	}
	else{
		qcf0 = (float*)calloc(lg, sizeof(float));
		printf("子带全为零\n");
	}

	scan_inv(a + 4, 2, 2, lg, m[23], 2, 1, 1, qcf0);	//反向扫描

	//5-------------------------------------------------------------------
	lg = (w[5] * h[4]) << 2;	//第一次PTV
	qctr = find_qctr(bin, &qnt, &trim, &ctr1);
	if (qnt){
		desub = de_sub3d_sub2(bin, lg, qctr, lenbits);	//解码
		free(desub.sn.dat); free(desub.cfw.dat);
	qcf0 = rst_sub(desub.cf.dat, delta, lg, qnt, qctr, trim, ctr1);	//反量化
	}
	else{
		qcf0 = (float*)calloc(lg, sizeof(float));
		printf("子带全为零\n");
	}

	scan_inv(a + 4, 0, 1, lg, m[24], 1, 1, 1, qcf0);	//反向扫描
	//6-------------------------------------------------------------------
	lg = (w[5] * h[5]) << 2;
	qctr = find_qctr(bin, &qnt, &trim, &ctr1);
	if (qnt){
		desub = de_sub3d_sub2(bin, lg, qctr, lenbits);	//解码
		free(desub.sn.dat); free(desub.cfw.dat);
	qcf0 = rst_sub(desub.cf.dat, delta, lg, qnt, qctr, trim, ctr1);	//反量化
	}
	else{
		qcf0 = (float*)calloc(lg, sizeof(float));
		printf("子带全为零\n");
	}

	scan_inv(a + 4, 1, 1, lg, m[26], 1, 1, 1, qcf0);	//反向扫描
	//7-------------------------------------------------------------------
	lg = (w[4] * h[5]) << 2;
	qctr = find_qctr(bin, &qnt, &trim, &ctr1);
	if (qnt){
		desub = de_sub3d_sub2(bin, lg, qctr, lenbits);	//解码
		free(desub.sn.dat); free(desub.cfw.dat);
	qcf0 = rst_sub(desub.cf.dat, delta, lg, qnt, qctr, trim, ctr1);	//反量化
	}
	else{
		qcf0 = (float*)calloc(lg, sizeof(float));
		printf("子带全为零\n");
	}

	scan_inv(a + 4, 1, 0, lg, m[25], 1, 1, 1, qcf0);	//反向扫描

	return;
}

void de_coed3d_7(int a, int *w, int *h, int **m, double delta, int lenbits)
{
	extern uchar *bin;
	extern int ptr;
	extern int len;
	int lg, qctr;
	double delta0;
	delta0 = 1.0 / delta;
	DE_S_SUB desub;
	float *qcf0;
	uchar qnt = 0, trim = 0;
	uint ctr1 = 0;

	lg = (w[6] * h[6]) << 2;
	/*1-------------------------------------------------------------------*/
	qctr = find_qctr(bin, &qnt, &trim, &ctr1);
	if (qnt){
		desub = de_sub3d_sub2(bin, lg, qctr, lenbits);	//解码
		free(desub.sn.dat); free(desub.cfw.dat);
	qcf0 = rst_sub(desub.cf.dat, delta, lg, qnt, qctr, trim, ctr1);
	}
	else{
		qcf0 = (float*)calloc(lg, sizeof(float));
		printf("子带全为零\n");
	}

	scan_inv(a + 4, 0, 1, lg, m[36], 3, 1, 0, qcf0);	//反向扫描

	/*2---------------------------------------------------------------------*/
	qctr = find_qctr(bin, &qnt, &trim, &ctr1);
	if (qnt){
		desub = de_sub3d_sub2(bin, lg, qctr, lenbits);	//解码
		free(desub.sn.dat); free(desub.cfw.dat);
	qcf0 = rst_sub(desub.cf.dat, delta, lg, qnt, qctr, trim, ctr1);
	}
	else{
		qcf0 = (float*)calloc(lg, sizeof(float));
		printf("子带全为零\n");
	}

	scan_inv(a + 4, 0, 2, lg, m[36], 3, 1, 0, qcf0);	//反向扫描

	/*3----------------------------------------------------------------------*/
	qctr = find_qctr(bin, &qnt, &trim, &ctr1);
	if (qnt){
		desub = de_sub3d_sub2(bin, lg, qctr, lenbits);	//解码
		free(desub.sn.dat); free(desub.cfw.dat);
	qcf0 = rst_sub(desub.cf.dat, delta, lg, qnt, qctr, trim, ctr1);
	}
	else{
		qcf0 = (float*)calloc(lg, sizeof(float));
		printf("子带全为零\n");
	}

	scan_inv(a + 4, 0, 3, lg, m[36], 3, 1, 0, qcf0);	//反向扫描

	/*4--------------------------------------------------------------------*/
	qctr = find_qctr(bin, &qnt, &trim, &ctr1);
	if (qnt){
		desub = de_sub3d_sub2(bin, lg, qctr, lenbits);	//解码
		free(desub.sn.dat); free(desub.cfw.dat);
	qcf0 = rst_sub(desub.cf.dat, delta, lg, qnt, qctr, trim, ctr1);
	}
	else{
		qcf0 = (float*)calloc(lg, sizeof(float));
		printf("子带全为零\n");
	}

	scan_inv(a + 4, 0, 4, lg, m[36], 3, 1, 0, qcf0);	//反向扫描

	/*5----------------------------------------------------------------------*/
	qctr = find_qctr(bin, &qnt, &trim, &ctr1);
	if (qnt){
		desub = de_sub3d_sub2(bin, lg, qctr, lenbits);	//解码
		free(desub.sn.dat); free(desub.cfw.dat);
	qcf0 = rst_sub(desub.cf.dat, delta, lg, qnt, qctr, trim, ctr1);
	}
	else{
		qcf0 = (float*)calloc(lg, sizeof(float));
		printf("子带全为零\n");
	}

	scan_inv(a + 4, 0, 5, lg, m[36], 3, 1, 0, qcf0);	//反向扫描

	/*6-------------------------------------------------------------------*/
	qctr = find_qctr(bin, &qnt, &trim, &ctr1);
	if (qnt){
		desub = de_sub3d_sub2(bin, lg, qctr, lenbits);	//解码
		free(desub.sn.dat); free(desub.cfw.dat);
	qcf0 = rst_sub(desub.cf.dat, delta, lg, qnt, qctr, trim, ctr1);
	}
	else{
		qcf0 = (float*)calloc(lg, sizeof(float));
		printf("子带全为零\n");
	}

	scan_inv(a + 4, 0, 6, lg, m[36], 3, 1, 0, qcf0);	//反向扫描

	/*7-------------------------------------------------------------------*/
	qctr = find_qctr(bin, &qnt, &trim, &ctr1);
	if (qnt){
		desub = de_sub3d_sub2(bin, lg, qctr, lenbits);	//解码
		free(desub.sn.dat); free(desub.cfw.dat);
	qcf0 = rst_sub(desub.cf.dat, delta, lg, qnt, qctr, trim, ctr1);
	}
	else{
		qcf0 = (float*)calloc(lg, sizeof(float));
		printf("子带全为零\n");
	}

	scan_inv(a + 4, 0, 7, lg, m[36], 3, 1, 0, qcf0);	//反向扫描

	return;
}

void de_coed3d_8(int a, int b, int *w, int *h, int **m, double delta, int lenbits)
{
	extern uchar *bin;
	extern int ptr;
	extern int len;
	int lg, qctr;
	double delta0;
	delta0 = 1.0 / delta;
	DE_S_SUB desub;
	uchar qnt = 0, trim = 0;
	uint ctr1 = 0;
	float *qcf0;

	lg = (w[6] * h[6]) << 2;
	/*1----------------------------------------------------------------*/
	qctr = find_qctr(bin, &qnt, &trim, &ctr1);
	if (qnt){
		desub = de_sub3d_sub2(bin, lg, qctr, lenbits);	//解码
		free(desub.sn.dat); free(desub.cfw.dat);
	qcf0 = rst_sub(desub.cf.dat, delta, lg, qnt, qctr, trim, ctr1);
	}
	else{
		qcf0 = (float*)calloc(lg, sizeof(float));
		printf("子带全为零\n");
	}

	scan_inv(a + 4, b, 0, lg, m[36], 3, 1, 0, qcf0);	//反向扫描

	/*2--------------------------------------------------------------------*/
	qctr = find_qctr(bin, &qnt, &trim, &ctr1);
	if (qnt){
		desub = de_sub3d_sub2(bin, lg, qctr, lenbits);	//解码
		free(desub.sn.dat); free(desub.cfw.dat);
	qcf0 = rst_sub(desub.cf.dat, delta, lg, qnt, qctr, trim, ctr1);
	}
	else{
		qcf0 = (float*)calloc(lg, sizeof(float));
		printf("子带全为零\n");
	}

	scan_inv(a + 4, b, 1, lg, m[36], 3, 1, 0, qcf0);	//反向扫描

	/*3--------------------------------------------------------------------*/
	qctr = find_qctr(bin, &qnt, &trim, &ctr1);
	if (qnt){
		desub = de_sub3d_sub2(bin, lg, qctr, lenbits);	//解码
		free(desub.sn.dat); free(desub.cfw.dat);
	qcf0 = rst_sub(desub.cf.dat, delta, lg, qnt, qctr, trim, ctr1);
	}
	else{
		qcf0 = (float*)calloc(lg, sizeof(float));
		printf("子带全为零\n");
	}

	scan_inv(a + 4, b, 2, lg, m[36], 3, 1, 0, qcf0);	//反向扫描

	/*4------------------------------------------------------------------*/
	qctr = find_qctr(bin, &qnt, &trim, &ctr1);
	if (qnt){
		desub = de_sub3d_sub2(bin, lg, qctr, lenbits);	//解码
		free(desub.sn.dat); free(desub.cfw.dat);
	qcf0 = rst_sub(desub.cf.dat, delta, lg, qnt, qctr, trim, ctr1);
	}
	else{
		qcf0 = (float*)calloc(lg, sizeof(float));
		printf("子带全为零\n");
	}

	scan_inv(a + 4, b, 3, lg, m[36], 3, 1, 0, qcf0);	//反向扫描

	/*5-------------------------------------------------------------------*/
	qctr = find_qctr(bin, &qnt, &trim, &ctr1);
	if (qnt){
		desub = de_sub3d_sub2(bin, lg, qctr, lenbits);	//解码
		free(desub.sn.dat); free(desub.cfw.dat);
	qcf0 = rst_sub(desub.cf.dat, delta, lg, qnt, qctr, trim, ctr1);
	}
	else{
		qcf0 = (float*)calloc(lg, sizeof(float));
		printf("子带全为零\n");
	}

	scan_inv(a + 4, b, 4, lg, m[36], 3, 1, 0, qcf0);	//反向扫描

	/*6-----------------------------------------------------------------*/
	qctr = find_qctr(bin, &qnt, &trim, &ctr1);
	if (qnt){
		desub = de_sub3d_sub2(bin, lg, qctr, lenbits);	//解码
		free(desub.sn.dat); free(desub.cfw.dat);
	qcf0 = rst_sub(desub.cf.dat, delta, lg, qnt, qctr, trim, ctr1);
	}
	else{
		qcf0 = (float*)calloc(lg, sizeof(float));
		printf("子带全为零\n");
	}

	scan_inv(a + 4, b, 5, lg, m[36], 3, 1, 0, qcf0);	//反向扫描

	/*7-----------------------------------------------------------------*/
	qctr = find_qctr(bin, &qnt, &trim, &ctr1);
	if (qnt){
		desub = de_sub3d_sub2(bin, lg, qctr, lenbits);	//解码
		free(desub.sn.dat); free(desub.cfw.dat);
	qcf0 = rst_sub(desub.cf.dat, delta, lg, qnt, qctr, trim, ctr1);
	}
	else{
		qcf0 = (float*)calloc(lg, sizeof(float));
		printf("子带全为零\n");
	}

	scan_inv(a + 4, b, 6, lg, m[36], 3, 1, 0, qcf0);	//反向扫描

	/*8------------------------------------------------------------------*/
	qctr = find_qctr(bin, &qnt, &trim, &ctr1);
	if (qnt){
		desub = de_sub3d_sub2(bin, lg, qctr, lenbits);	//解码
		free(desub.sn.dat); free(desub.cfw.dat);
	qcf0 = rst_sub(desub.cf.dat, delta, lg, qnt, qctr, trim, ctr1);
	}
	else{
		qcf0 = (float*)calloc(lg, sizeof(float));
		printf("子带全为零\n");
	}

	scan_inv(a + 4, b, 7, lg, m[36], 3, 1, 0, qcf0);	//反向扫描

	return;
}

void de_coef3dB_dc(int a, int *w, int *h, int **m, double delta, int lenbits)
{
	extern uchar *bin;
	extern int ptr;
	extern float ***rePTVData;
	extern int len;
	int lg, qctr;
	double delta0;
	DE_S_SUB desub;
	uchar qnt = 0, trim = 0;
	uint ctr1 = 0;
	float *qcf0;
	delta0 = 1.0 / delta;

	//1-------------------------------------------------------------------
	lg = (w[0] * h[0]) << 1;//第三次PTV
	qctr = find_qctr(bin, &qnt, &trim, &ctr1);
	if (qnt){
		desub = de_sub3d_sub2(bin, lg, qctr, lenbits);	//解码
		free(desub.sn.dat); free(desub.cfw.dat);
	qcf0 = rst_sub(desub.cf.dat, delta, lg, qnt, qctr, trim, ctr1);	//反量化
	}
	else{
		qcf0 = (float*)calloc(lg, sizeof(float));
		printf("子带全为零\n");
	}

	scan_inv(a + 8, 0, 0, lg, m[10], 3, 2, 1, qcf0);	//反向扫描

	//2-------------------------------------------------------------------
	lg = (w[1] * h[0]) << 1;
	qctr = find_qctr(bin, &qnt, &trim, &ctr1);
	if (qnt){
		desub = de_sub3d_sub2(bin, lg, qctr, lenbits);	//解码
		free(desub.sn.dat); free(desub.cfw.dat);
	qcf0 = rst_sub(desub.cf.dat, delta, lg, qnt, qctr, trim, ctr1);	//反量化
	}
	else{
		qcf0 = (float*)calloc(lg, sizeof(float));
		printf("子带全为零\n");
	}

	scan_inv(a + 8, 0, 4, lg, m[11], 3, 2, 1, qcf0);	//反向扫描
	//3-------------------------------------------------------------------
	lg = (w[1] * h[1]) << 1;
	qctr = find_qctr(bin, &qnt, &trim, &ctr1);
	if (qnt){
		desub = de_sub3d_sub2(bin, lg, qctr, lenbits);	//解码
		free(desub.sn.dat); free(desub.cfw.dat);
	qcf0 = rst_sub(desub.cf.dat, delta, lg, qnt, qctr, trim, ctr1);	//反量化
	}
	else{
		qcf0 = (float*)calloc(lg, sizeof(float));
		printf("子带全为零\n");
	}

	scan_inv(a + 8, 4, 4, lg, m[13], 3, 2, 1, qcf0);	//反向扫描
	//4-------------------------------------------------------------------
	lg = (w[0] * h[1]) << 1;
	qctr = find_qctr(bin, &qnt, &trim, &ctr1);
	if (qnt){
		desub = de_sub3d_sub2(bin, lg, qctr, lenbits);	//解码
		free(desub.sn.dat); free(desub.cfw.dat);
	qcf0 = rst_sub(desub.cf.dat, delta, lg, qnt, qctr, trim, ctr1);	//反量化
	}
	else{
		qcf0 = (float*)calloc(lg, sizeof(float));
		printf("子带全为零\n");
	}

	scan_inv(a + 8, 4, 0, lg, m[12], 3, 2, 1, qcf0);	//反向扫描

	//5-------------------------------------------------------------------
	lg = (w[3] * h[2]) << 1;//第二次PTV
	qctr = find_qctr(bin, &qnt, &trim, &ctr1);
	if (qnt){
		desub = de_sub3d_sub2(bin, lg, qctr, lenbits);	//解码
		free(desub.sn.dat); free(desub.cfw.dat);
	qcf0 = rst_sub(desub.cf.dat, delta, lg, qnt, qctr, trim, ctr1);	//反量化
	}
	else{
		qcf0 = (float*)calloc(lg, sizeof(float));
		printf("子带全为零\n");
	}

	scan_inv(a + 8, 0, 2, lg, m[14], 2, 2, 1, qcf0);	//反向扫描
	//6-------------------------------------------------------------------
	lg = (w[3] * h[3]) << 1;
	qctr = find_qctr(bin, &qnt, &trim, &ctr1);
	if (qnt){
		desub = de_sub3d_sub2(bin, lg, qctr, lenbits);	//解码
		free(desub.sn.dat); free(desub.cfw.dat);
	qcf0 = rst_sub(desub.cf.dat, delta, lg, qnt, qctr, trim, ctr1);	//反量化
	}
	else{
		qcf0 = (float*)calloc(lg, sizeof(float));
		printf("子带全为零\n");
	}

	scan_inv(a + 8, 2, 2, lg, m[16], 2, 2, 1, qcf0);	//反向扫描
	//7-------------------------------------------------------------------
	lg = (w[2] * h[3]) << 1;
	qctr = find_qctr(bin, &qnt, &trim, &ctr1);
	if (qnt){
		desub = de_sub3d_sub2(bin, lg, qctr, lenbits);	//解码
		free(desub.sn.dat); free(desub.cfw.dat);
	qcf0 = rst_sub(desub.cf.dat, delta, lg, qnt, qctr, trim, ctr1);	//反量化
	}
	else{
		qcf0 = (float*)calloc(lg, sizeof(float));
		printf("子带全为零\n");
	}

	scan_inv(a + 8, 2, 0, lg, m[15], 2, 2, 1, qcf0);	//反向扫描

	//8-------------------------------------------------------------------
	lg = (w[5] * h[4]) << 1;//第一次PTV
	qctr = find_qctr(bin, &qnt, &trim, &ctr1);
	if (qnt){
		desub = de_sub3d_sub2(bin, lg, qctr, lenbits);	//解码
		free(desub.sn.dat); free(desub.cfw.dat);
	qcf0 = rst_sub(desub.cf.dat, delta, lg, qnt, qctr, trim, ctr1);	//反量化
	}
	else{
		qcf0 = (float*)calloc(lg, sizeof(float));
		printf("子带全为零\n");
	}

	scan_inv(a + 8, 0, 1, lg, m[17], 1, 2, 1, qcf0);	//反向扫描
	//9-------------------------------------------------------------------
	lg = (w[5] * h[5]) << 1;
	qctr = find_qctr(bin, &qnt, &trim, &ctr1);
	if (qnt){
		desub = de_sub3d_sub2(bin, lg, qctr, lenbits);	//解码
		free(desub.sn.dat); free(desub.cfw.dat);
	qcf0 = rst_sub(desub.cf.dat, delta, lg, qnt, qctr, trim, ctr1);	//反量化
	}
	else{
		qcf0 = (float*)calloc(lg, sizeof(float));
		printf("子带全为零\n");
	}

	scan_inv(a + 8, 1, 1, lg, m[19], 1, 2, 1, qcf0);	//反向扫描
	//10-------------------------------------------------------------------
	lg = (w[4] * h[5]) << 1;
	qctr = find_qctr(bin, &qnt, &trim, &ctr1);
	if (qnt){
		desub = de_sub3d_sub2(bin, lg, qctr, lenbits);	//解码
		free(desub.sn.dat); free(desub.cfw.dat);
	qcf0 = rst_sub(desub.cf.dat, delta, lg, qnt, qctr, trim, ctr1);	//反量化
	}
	else{
		qcf0 = (float*)calloc(lg, sizeof(float));
		printf("子带全为零\n");
	}

	scan_inv(a + 8, 1, 0, lg, m[18], 1, 2, 1, qcf0);	//反向扫描


	return;
}

void de_coef3dB_L3(int H, int W, int a, int b, int c, int *w, int *h, int **m, double delta, int lenbits)
{
	extern uchar *bin;
	extern int ptr;
	extern int len;
	int lg, qctr;
	double delta0;
	delta0 = 1.0 / delta;
	DE_S_SUB desub;
	uchar qnt = 0, trim = 0;
	uint ctr1 = 0;
	float *qcf0;

	lg = (w[6] * h[6]) << 1;
	/*1-------------------------------------------------------------------------------*/
	qctr = find_qctr(bin, &qnt, &trim, &ctr1);
	if (qnt){
		desub = de_sub3d_sub2(bin, lg, qctr, lenbits);	//解码
		free(desub.sn.dat); free(desub.cfw.dat);
	qcf0 = rst_sub(desub.cf.dat, delta, lg, qnt, qctr, trim, ctr1);
	}
	else{
		qcf0 = (float*)calloc(lg, sizeof(float));
		printf("子带全为零\n");
	}

	scan_inv(a + 8, b, c, lg, m[35], 3, 2, 0, qcf0);	//反向扫描
	/*2-------------------------------------------------------------------------------*/
	qctr = find_qctr(bin, &qnt, &trim, &ctr1);
	if (qnt){
		desub = de_sub3d_sub2(bin, lg, qctr, lenbits);	//解码
		free(desub.sn.dat); free(desub.cfw.dat);
	qcf0 = rst_sub(desub.cf.dat, delta, lg, qnt, qctr, trim, ctr1);
	}
	else{
		qcf0 = (float*)calloc(lg, sizeof(float));
		printf("子带全为零\n");
	}

	scan_inv(a + 8, b+1, c, lg, m[35], 3, 2, 0, qcf0);	//反向扫描
	/*3-------------------------------------------------------------------------------*/
	qctr = find_qctr(bin, &qnt, &trim, &ctr1);
	if (qnt){
		desub = de_sub3d_sub2(bin, lg, qctr, lenbits);	//解码
		free(desub.sn.dat); free(desub.cfw.dat);
	qcf0 = rst_sub(desub.cf.dat, delta, lg, qnt, qctr, trim, ctr1);
	}
	else{
		qcf0 = (float*)calloc(lg, sizeof(float));
		printf("子带全为零\n");
	}

	scan_inv(a + 8, b+1, c-1, lg, m[35], 3, 2, 0, qcf0);	//反向扫描


	return;
}

void de_coef3dB_L2_3(int H, int W, int a, int b, int c, int *w, int *h, int **m, double delta, int lenbits)
{
	extern uchar *bin;
	extern int ptr;
	extern int len;
	int lg, qctr;
	double delta0;
	delta0 = 1.0 / delta;
	DE_S_SUB desub;
	uchar qnt = 0, trim = 0;
	uint ctr1 = 0;
	float *qcf0;

	lg = (2 * w[6] * 2 * h[6]) << 1;
	/*1-------------------------------------------------------------------------------*/
	qctr = find_qctr(bin, &qnt, &trim, &ctr1);
	if (qnt){
		desub = de_sub3d_sub2(bin, lg, qctr, lenbits);	//解码
		free(desub.sn.dat); free(desub.cfw.dat);
	qcf0 = rst_sub(desub.cf.dat, delta, lg, qnt, qctr, trim, ctr1);
	}
	else{
		qcf0 = (float*)calloc(lg, sizeof(float));
		printf("子带全为零\n");
	}

	scan_inv_L2(H, W, a + 8, b, c, lg, m[38], 2, qcf0);
	/*2-------------------------------------------------------------------------------*/
	qctr = find_qctr(bin, &qnt, &trim, &ctr1);
	if (qnt){
		desub = de_sub3d_sub2(bin, lg, qctr, lenbits);	//解码
		free(desub.sn.dat); free(desub.cfw.dat);
	qcf0 = rst_sub(desub.cf.dat, delta, lg, qnt, qctr, trim, ctr1);
	}
	else{
		qcf0 = (float*)calloc(lg, sizeof(float));
		printf("子带全为零\n");
	}

	scan_inv_L2(H, W, a + 8, b, c+2, lg, m[38], 2, qcf0);	//反向扫描
	/*3-------------------------------------------------------------------------------*/
	qctr = find_qctr(bin, &qnt, &trim, &ctr1);
	if (qnt){
		desub = de_sub3d_sub2(bin, lg, qctr, lenbits);	//解码
		free(desub.sn.dat); free(desub.cfw.dat);
	qcf0 = rst_sub(desub.cf.dat, delta, lg, qnt, qctr, trim, ctr1);
	}
	else{
		qcf0 = (float*)calloc(lg, sizeof(float));
		printf("子带全为零\n");
	}

	scan_inv_L2(H, W, a + 8, b, c+4, lg, m[38], 2, qcf0);	//反向扫描

	return;
}

void de_coef3dB_L2_4(int H, int W, int a, int b, int c, int *w, int *h, int **m, double delta, int lenbits)
{
	extern uchar *bin;
	extern int ptr;
	extern int len;
	int lg, qctr;
	double delta0;
	delta0 = 1.0 / delta;
	DE_S_SUB desub;
	uchar qnt = 0, trim = 0;
	uint ctr1 = 0;
	float *qcf0;

	lg = (2 * w[6] * 2 * h[6]) << 1;
	/*1-------------------------------------------------------------------------------*/
	qctr = find_qctr(bin, &qnt, &trim, &ctr1);
	if (qnt){
		desub = de_sub3d_sub2(bin, lg, qctr, lenbits);	//解码
		free(desub.sn.dat); free(desub.cfw.dat);
	qcf0 = rst_sub(desub.cf.dat, delta, lg, qnt, qctr, trim, ctr1);
	}
	else{
		qcf0 = (float*)calloc(lg, sizeof(float));
		printf("子带全为零\n");
	}

	scan_inv_L2(H, W, a + 8, b, c, lg, m[38], 2, qcf0);
	/*2-------------------------------------------------------------------------------*/
	qctr = find_qctr(bin, &qnt, &trim, &ctr1);
	if (qnt){
		desub = de_sub3d_sub2(bin, lg, qctr, lenbits);	//解码
		free(desub.sn.dat); free(desub.cfw.dat);
	qcf0 = rst_sub(desub.cf.dat, delta, lg, qnt, qctr, trim, ctr1);
	}
	else{
		qcf0 = (float*)calloc(lg, sizeof(float));
		printf("子带全为零\n");
	}

	scan_inv_L2(H, W, a + 8, b, c + 2, lg, m[38], 2, qcf0);	//反向扫描
	/*3-------------------------------------------------------------------------------*/
	qctr = find_qctr(bin, &qnt, &trim, &ctr1);
	if (qnt){
		desub = de_sub3d_sub2(bin, lg, qctr, lenbits);	//解码
		free(desub.sn.dat); free(desub.cfw.dat);
	qcf0 = rst_sub(desub.cf.dat, delta, lg, qnt, qctr, trim, ctr1);
	}
	else{
		qcf0 = (float*)calloc(lg, sizeof(float));
		printf("子带全为零\n");
	}

	scan_inv_L2(H, W, a + 8, b, c + 4, lg, m[38], 2, qcf0);	//反向扫描
	/*4-------------------------------------------------------------------------------*/
	qctr = find_qctr(bin, &qnt, &trim, &ctr1);
	if (qnt){
		desub = de_sub3d_sub2(bin, lg, qctr, lenbits);	//解码
		free(desub.sn.dat); free(desub.cfw.dat);
	qcf0 = rst_sub(desub.cf.dat, delta, lg, qnt, qctr, trim, ctr1);
	}
	else{
		qcf0 = (float*)calloc(lg, sizeof(float));
		printf("子带全为零\n");
	}

	scan_inv_L2(H, W, a + 8, b, c + 6, lg, m[38], 2, qcf0);	//反向扫描

	return;
}

void de_coef3d_dc5B_dc(int a, int *w, int *h, int **m, double delta, int AC, int lenbits)
{
	extern uchar *bin;
	extern int ptr;
	extern float ***rePTVData;
	extern int **nc;
	extern int len;
	int lg, qctr;
	double delta0;

	DE_S_SUB desub;
	uchar qnt = 0, trim = 0;
	uint ctr1 = 0;
	float *qcf0;
	delta0 = 1.0 / delta;

	if (!AC){
		//DC--------------------------------------------------------------------
		lg = w[0] * h[0];//第三次PTV，对应matlab里Ldc=coef(:,:,2)
		qcf0 = de_DCinv(bin, a, lg, m[0], delta);
	}
	else{
		lg = w[0] * h[0];
		qctr = find_qctr(bin, &qnt, &trim, &ctr1);
		if (qnt){
			desub = de_sub3d_sub2(bin, lg, qctr, lenbits);	//解码
			free(desub.sn.dat); free(desub.cfw.dat);
		qcf0 = rst_sub(desub.cf.dat, delta, lg, qnt, qctr, trim, ctr1);	//反量化
		}
		else{
			qcf0 = (float*)calloc(lg, sizeof(float));
			printf("子带全为零\n");
		}

		scan_inv(a, 0, 0, lg, m[0], 3, 3, 1, qcf0);	//反向扫描
	}

	//2---------------------------------------------------------------------
	lg = w[1] * h[0];
	qctr = find_qctr(bin, &qnt, &trim, &ctr1);
	if (qnt){
		desub = de_sub3d_sub2(bin, lg, qctr, lenbits);
		free(desub.sn.dat); free(desub.cfw.dat);
	qcf0 = rst_sub(desub.cf.dat, delta, lg, qnt, qctr, trim, ctr1);	//反量化
	}
	else{
		qcf0 = (float*)calloc(lg, sizeof(float));
		printf("子带全为零\n");
	}

	scan_inv(a, 0, 4, lg, m[1], 3, 3, 1, qcf0);	//反向扫描

	//3--------------------------------------------------------------------------
	lg = w[1] * h[1];
	qctr = find_qctr(bin, &qnt, &trim, &ctr1);
	if (qnt){
		desub = de_sub3d_sub2(bin, lg, qctr, lenbits);	//解码
		free(desub.sn.dat); free(desub.cfw.dat);
	qcf0 = rst_sub(desub.cf.dat, delta, lg, qnt, qctr, trim, ctr1);	//反量化
	}
	else{
		qcf0 = (float*)calloc(lg, sizeof(float));
		printf("子带全为零\n");
	}

	scan_inv(a, 4, 4, lg, m[3], 3, 3, 1, qcf0);	//反向扫描

	//4-------------------------------------------------------------------
	lg = w[0] * h[1];
	qctr = find_qctr(bin, &qnt, &trim, &ctr1);
	if (qnt){
		desub = de_sub3d_sub2(bin, lg, qctr, lenbits);	//解码
		free(desub.sn.dat); free(desub.cfw.dat);
	qcf0 = rst_sub(desub.cf.dat, delta, lg, qnt, qctr, trim, ctr1);	//反量化
	}
	else{
		qcf0 = (float*)calloc(lg, sizeof(float));
		printf("子带全为零\n");
	}

	scan_inv(a, 4, 0, lg, m[2], 3, 3, 1, qcf0);	//反向扫描

	//5--------------------------------------------------------------------
	lg = w[3] * h[2];//第二次PTV，对应matlab里L5=coef(:,:,5)
	qctr = find_qctr(bin, &qnt, &trim, &ctr1);
	if (qnt){
		desub = de_sub3d_sub2(bin, lg, qctr, lenbits);	//解码
		free(desub.sn.dat); free(desub.cfw.dat);
	qcf0 = rst_sub(desub.cf.dat, delta, lg, qnt, qctr, trim, ctr1);	//反量化
	}
	else{
		qcf0 = (float*)calloc(lg, sizeof(float));
		printf("子带全为零\n");
	}

	scan_inv(a, 0, 2, lg, m[4], 2, 3, 1, qcf0);	//反向扫描
	//7--------------------------------------------------------------------
	lg = w[3] * h[3];
	qctr = find_qctr(bin, &qnt, &trim, &ctr1);
	if (qnt){
		desub = de_sub3d_sub2(bin, lg, qctr, lenbits);	//解码
		free(desub.sn.dat); free(desub.cfw.dat);
	qcf0 = rst_sub(desub.cf.dat, delta, lg, qnt, qctr, trim, ctr1);	//反量化
	}
	else{
		qcf0 = (float*)calloc(lg, sizeof(float));
		printf("子带全为零\n");
	}

	scan_inv(a, 2, 2, lg, m[6], 2, 3, 1, qcf0);	//反向扫描
	//6-------------------------------------------------------------------
	lg = w[2] * h[3];
	qctr = find_qctr(bin, &qnt, &trim, &ctr1);
	if (qnt){
		desub = de_sub3d_sub2(bin, lg, qctr, lenbits);	//解码
		free(desub.sn.dat); free(desub.cfw.dat);
	qcf0 = rst_sub(desub.cf.dat, delta, lg, qnt, qctr, trim, ctr1);	//反量化
	}
	else{
		qcf0 = (float*)calloc(lg, sizeof(float));
		printf("子带全为零\n");
	}

	scan_inv(a, 2, 0, lg, m[5], 2, 3, 1, qcf0);	//反向扫描

	//8-------------------------------------------------------------------
	lg = w[5] * h[4];//第一次PTV
	qctr = find_qctr(bin, &qnt, &trim, &ctr1);
	if (qnt){
		desub = de_sub3d_sub2(bin, lg, qctr, lenbits);	//解码
		free(desub.sn.dat); free(desub.cfw.dat);
	qcf0 = rst_sub(desub.cf.dat, delta, lg, qnt, qctr, trim, ctr1);	//反量化
	}
	else{
		qcf0 = (float*)calloc(lg, sizeof(float));
		printf("子带全为零\n");
	}

	scan_inv(a, 0, 1, lg, m[7], 1, 3, 1, qcf0);	//反向扫描
	//10---------------------------------------------------------------------
	lg = w[5] * h[5];
	qctr = find_qctr(bin, &qnt, &trim, &ctr1);
	if (qnt){
		desub = de_sub3d_sub2(bin, lg, qctr, lenbits);	//解码
		free(desub.sn.dat); free(desub.cfw.dat);
	qcf0 = rst_sub(desub.cf.dat, delta, lg, qnt, qctr, trim, ctr1);	//反量化
	}
	else{
		qcf0 = (float*)calloc(lg, sizeof(float));
		printf("子带全为零\n");
	}

	scan_inv(a, 1, 1, lg, m[9], 1, 3, 1, qcf0);	//反向扫描
	//9---------------------------------------------------------------------
	lg = w[4] * h[5];
	qctr = find_qctr(bin, &qnt, &trim, &ctr1);
	if (qnt){
		desub = de_sub3d_sub2(bin, lg, qctr, lenbits);	//解码
		free(desub.sn.dat); free(desub.cfw.dat);
	qcf0 = rst_sub(desub.cf.dat, delta, lg, qnt, qctr, trim, ctr1);	//反量化
	}
	else{
		qcf0 = (float*)calloc(lg, sizeof(float));
		printf("子带全为零\n");
	}

	scan_inv(a, 1, 0, lg, m[8], 1, 3, 1, qcf0);	//反向扫描



	return;
}

void de_coef3d_dc5B_L3(int H, int W, int a, int b, int c, int *w, int *h, int **m, double delta, int lenbits)
{
	extern uchar *bin;
	extern int ptr;
	extern int **nc;
	extern int len;
	int lg, qctr;
	double delta0;

	DE_S_SUB desub;
	uchar qnt = 0, trim = 0;
	uint ctr1 = 0;
	float *qcf0;
	delta0 = 1.0 / delta;

	lg = (w[6] * h[6]);
	/*1------------------------------------------------------------------------*/
	qctr = find_qctr(bin, &qnt, &trim, &ctr1);
	if (qnt){
		desub = de_sub3d_sub2(bin, lg, qctr, lenbits);
		free(desub.sn.dat); free(desub.cfw.dat);
	qcf0 = rst_sub(desub.cf.dat, delta, lg, qnt, qctr, trim, ctr1);	//反量化
	}
	else{
		qcf0 = (float*)calloc(lg, sizeof(float));
		printf("子带全为零\n");
	}

	scan_inv(a, b, c, lg, m[34], 3, 3, 0, qcf0);	//反向扫描
	/*2------------------------------------------------------------------------*/
	qctr = find_qctr(bin, &qnt, &trim, &ctr1);
	if (qnt){
		desub = de_sub3d_sub2(bin, lg, qctr, lenbits);
		free(desub.sn.dat); free(desub.cfw.dat);
	qcf0 = rst_sub(desub.cf.dat, delta, lg, qnt, qctr, trim, ctr1);	//反量化
	}
	else{
		qcf0 = (float*)calloc(lg, sizeof(float));
		printf("子带全为零\n");
	}

	scan_inv(a, b+1, c, lg, m[34], 3, 3, 0, qcf0);	//反向扫描
	/*3------------------------------------------------------------------------*/
	qctr = find_qctr(bin, &qnt, &trim, &ctr1);
	if (qnt){
		desub = de_sub3d_sub2(bin, lg, qctr, lenbits);
		free(desub.sn.dat); free(desub.cfw.dat);
	qcf0 = rst_sub(desub.cf.dat, delta, lg, qnt, qctr, trim, ctr1);	//反量化
	}
	else{
		qcf0 = (float*)calloc(lg, sizeof(float));
		printf("子带全为零\n");
	}

	scan_inv(a, b+1, c-1, lg, m[34], 3, 3, 0, qcf0);	//反向扫描

	return;
}

void de_coef3d_dc5B_L2(int H, int W, int a, int b, int c, int *w, int *h, int **m, double delta, int lenbits)
{
	extern uchar *bin;
	extern int ptr;
	extern int **nc;
	extern int len;
	int lg, qctr;
	double delta0;

	DE_S_SUB desub;
	uchar qnt = 0, trim = 0;
	uint ctr1 = 0;
	float *qcf0;
	delta0 = 1.0 / delta;

	lg = (2 * w[6] * 2 * h[6]);
	/*1------------------------------------------------------------------------*/
	qctr = find_qctr(bin, &qnt, &trim, &ctr1);
	if (qnt){
		desub = de_sub3d_sub2(bin, lg, qctr, lenbits);
		free(desub.sn.dat); free(desub.cfw.dat);
	qcf0 = rst_sub(desub.cf.dat, delta, lg, qnt, qctr, trim, ctr1);	//反量化
	}
	else{
		qcf0 = (float*)calloc(lg, sizeof(float));
		printf("子带全为零\n");
	}

	scan_inv_L2(H, W, a, b, c, lg, m[39], 3, qcf0);
	/*2------------------------------------------------------------------------*/
	qctr = find_qctr(bin, &qnt, &trim, &ctr1);
	if (qnt){
		desub = de_sub3d_sub2(bin, lg, qctr, lenbits);
		free(desub.sn.dat); free(desub.cfw.dat);
	qcf0 = rst_sub(desub.cf.dat, delta, lg, qnt, qctr, trim, ctr1);	//反量化
	}
	else{
		qcf0 = (float*)calloc(lg, sizeof(float));
		printf("子带全为零\n");
	}

	scan_inv_L2(H, W, a, b+2, c, lg, m[39], 3, qcf0);	//反向扫描
	/*3------------------------------------------------------------------------*/
	qctr = find_qctr(bin, &qnt, &trim, &ctr1);
	if (qnt){
		desub = de_sub3d_sub2(bin, lg, qctr, lenbits);
		free(desub.sn.dat); free(desub.cfw.dat);
	qcf0 = rst_sub(desub.cf.dat, delta, lg, qnt, qctr, trim, ctr1);	//反量化
	}
	else{
		qcf0 = (float*)calloc(lg, sizeof(float));
		printf("子带全为零\n");
	}

	scan_inv_L2(H, W, a, b+2, c-2, lg, m[39], 3, qcf0);	//反向扫描

	return;
}

void de_coef3d_dc5B_L1(int H, int W, int a, int b, int c, int *w, int *h, int **m, double delta, int lenbits)
{
	extern uchar *bin;
	extern int ptr;
	extern int **nc;
	extern int len;
	int lg, qctr;
	double delta0;

	DE_S_SUB desub;
	uchar qnt = 0, trim = 0;
	uint ctr1 = 0;
	float *qcf0;
	delta0 = 1.0 / delta;

	lg = (4 * w[6] * 4 * h[6]);
	/*1------------------------------------------------------------------------*/
	qctr = find_qctr(bin, &qnt, &trim, &ctr1);
	if (qnt){
		desub = de_sub3d_sub2(bin, lg, qctr, lenbits);
		free(desub.sn.dat); free(desub.cfw.dat);
	qcf0 = rst_sub(desub.cf.dat, delta, lg, qnt, qctr, trim, ctr1);	//反量化
	}
	else{
		qcf0 = (float*)calloc(lg, sizeof(float));
		printf("子带全为零\n");
	}

	scan_inv_L2(H, W, a, b, c, lg, m[40], 3, qcf0);
	/*2------------------------------------------------------------------------*/
	qctr = find_qctr(bin, &qnt, &trim, &ctr1);
	if (qnt){
		desub = de_sub3d_sub2(bin, lg, qctr, lenbits);
		free(desub.sn.dat); free(desub.cfw.dat);
	qcf0 = rst_sub(desub.cf.dat, delta, lg, qnt, qctr, trim, ctr1);	//反量化
	}
	else{
		qcf0 = (float*)calloc(lg, sizeof(float));
		printf("子带全为零\n");
	}

	scan_inv_L2(H, W, a, b + 4, c, lg, m[40], 3, qcf0);	//反向扫描
	/*3------------------------------------------------------------------------*/
	qctr = find_qctr(bin, &qnt, &trim, &ctr1);
	if (qnt){
		desub = de_sub3d_sub2(bin, lg, qctr, lenbits);
		free(desub.sn.dat); free(desub.cfw.dat);
	qcf0 = rst_sub(desub.cf.dat, delta, lg, qnt, qctr, trim, ctr1);	//反量化
	}
	else{
		qcf0 = (float*)calloc(lg, sizeof(float));
		printf("子带全为零\n");
	}

	scan_inv_L2(H, W, a, b + 4, c - 4, lg, m[40], 3, qcf0);	//反向扫描

	return;
}
