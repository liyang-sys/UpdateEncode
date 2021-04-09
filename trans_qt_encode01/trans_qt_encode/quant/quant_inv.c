#include "quant_inv.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "all.h"

float* rst_sub(int* qcf, double delta, int lg, uchar qnt, int qctr, uchar trim, uint ctr1)
{
	extern int len;
	extern float ***rePTVData;
	int max=0, temp;
	int *nc = (int*)calloc(len, sizeof(int));
	char *sgn = (char*)calloc(lg, sizeof(char));
	float *qf = NULL;
	double T;

	if (qnt){
		if (trim){	//trim
			T = 0.55 * delta;
			for (int i = 0; i < lg; i++){
				if (qcf[i] & 0x80000000){
					sgn[i]--;
					qcf[i] *= -1;
				}
				else
					sgn[i]++;
				temp = qcf[i] + 1;
				if (temp>max)
					max = temp;
				nc[temp]++;
			}
			qf = rstTHDctr(qcf, sgn, T, delta, nc, max, lg);
		}
		else{
			if (qctr){	//quant center
				for (int i = 0; i < lg; i++){
					if (qcf[i] & 0x80000000){
						sgn[i]--;
						qcf[i] *= -1;
					}
					else
						sgn[i]++;
					temp = qcf[i] + 1;
					if (temp>max)
						max = temp;
					nc[temp]++;
				}
				qf = rstTHDctr(qcf, sgn, 0.52*delta, delta, nc, max, lg);
			}
			else{	//quant even
				for (int i = 0; i < lg; i++){
					if (qcf[i] & 0x80000000){
						sgn[i]--;
						qcf[i] *= -1;
					}
					else
						sgn[i]++;
					temp = qcf[i];
					if (temp>max)
						max = temp;
					nc[temp]++;
				}
				qf = rstEVENctr(qcf, sgn, delta, nc, max, (double)ctr1-1.0, lg);
			}
		}

	}
	else{
		qf = (float*)calloc(lg, sizeof(float));
	}
	free(nc);
	free(qcf); free(sgn);
	return qf;
}

float* rstTHDctr(int* qcf, char *sgn, double T, double delta, int* nc, int maxcf, int lg)
{
	double a3 = 1.0 / 3.0, b3 = 2.0 / 3.0, na, nb;
	int maxval = maxcf - 1;
	int valT = 1, val, idx=1;
	int minnc = nc[1];
	float* qf = (float*)malloc(lg*sizeof(float));
	for (int i = 2; i < maxcf + 1; i++){
		if (nc[i] < minnc){
			minnc = nc[i];
			idx = i;
		}
	}
	
	for (int i = 1; i < idx - 1; i++){
		na = (double)nc[i + 1]; nb = (double)nc[i + 2];
		if (nb>na){
			valT = i;
			break;
		}
		else
			valT = i + 1;
	}
	
	for (int i = 0; i < lg; i++){
		val = qcf[i];
		if (!val){
			qf[i] = 0;
		}
		else if (val < valT){
			na = (double)nc[val + 1]; nb = (double)nc[val + 2];
			double ctr_geo = (a3*na + b3*nb) / (na + nb)*delta;
			qf[i] = (T + ((double)val - 1)*delta + ctr_geo);
		}
		else{
			if (maxval>1){
				qf[i] = (T + (val - 0.5)*delta);
			}
			else{
				qf[i] = (T + 0.37*delta);
			}
		}
		qf[i] *= sgn[i];
	}

	return qf;
}

float* rstEVENctr(int* qcf, char *sgn, double delta, int* nc, int maxcf, double ctr1, int lg)
{
	double na, nb, a3 = 1.0 / 3.0, b3 = 2.0 / 3.0;
	int maxval = maxcf, minnc = nc[1], val, idx = 1;
	float* qf = (float*)malloc(lg*sizeof(float));

	for (int i = 2; i < maxcf + 1; i++){
		if (nc[i] < minnc){
			minnc = nc[i];
			idx = i;
		}
	}
	ctr1 = (float)((float)ctr1*0.125 / 64.0 + 0.375)*delta;
	for (int i = 0; i < lg; i++){
		val = qcf[i];
		qf[i] = ((float)val - 1.0)*delta;
		if (val == 1){
			qf[i] += ctr1;
		}
		else if (val < idx - 1){
			na = (double)nc[val]; nb = (double)nc[val + 1];
			double ctr_geo = (a3*na + b3*nb) / (na + nb)*delta;
			qf[i] += ctr_geo;
		}
		else{
			qf[i] += (0.5*delta);
		}
		qf[i] *= sgn[i];
	}


	return qf;
}

void scan_inv(int a, int b, int c, int lg, int *m, int ptv, int harr, int dc, float *qcf)
{
	extern float ***rePTVData;
	extern float ***reVideoData;
	float ***data;
	int offset0, offset1, i = 0;
	if (!qcf){
		printf("qcf为空，或qcf全为零\n");
		return;
	}
	if (dc){
		data = rePTVData;
		offset0 = 1 << ptv;
	}
	else{
		data = reVideoData;
		offset0 = 8;
	}
	offset1 = 4 << harr;
	for (int j = 0; j < lg; j++){
		data[a][b][c] = qcf[j];
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
	free(qcf);
	
	return;
}

void scan_inv_L2(int H, int W, int a, int b, int c, int lg, int *m, int harr, float *qcf)
{
	extern float ***reVideoData;
	float ***data;
	int offset0, offset1, i = 0, tempb, tempc;
	if (!qcf){
		printf("qcf为空，或qcf全为零\n");
		return;
	}
	data = reVideoData;
	offset0 = 8;
	offset1 = 4 << harr;
	tempb = b; tempc = c;
	for (int j = 0; j < lg; j++){
		data[a][b][c] = qcf[j];
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

		if (b >= H){
			tempb++;
			b = tempb;
		}
		else if (b < 0){
			tempb--;
			b = tempb + H - 8;
		}
		else if (c >= W){
			tempc++;
			c = tempc;
		}
		else if (c < 0){
			tempc--;
			c = tempc + W - 8;
		}

	}
	free(qcf);

	return;
}

int find_qctr(uchar *bin, uchar *qnt, uchar *trim, uint *ctr1)
{
	union data {
		unsigned long long a;
		uchar b[8];
	} rem1;
	extern int ptr;
	int qctr, x;
	x = ptr & 7;
	rem1.a = 0;
	uchar qnt0, trim0;
	int idx = ptr >> 3;

	rem1.b[5] = bin[idx];
	rem1.b[4] = bin[idx + 1];
	rem1.b[3] = bin[idx + 2];
	rem1.b[2] = bin[idx + 3];
	rem1.b[1] = bin[idx + 4];
	rem1.b[0] = bin[idx + 5];
	rem1.a <<= x;
	rem1.b[6] = 0;
	rem1.a <<= 1; ptr++;
	qnt0 = (uchar)rem1.b[6];
	rem1.b[6] = 0;
	
	*qnt = qnt0;
	if (qnt0){
		rem1.a <<= 1; ptr++;
		trim0 = (uchar)rem1.b[6];
		*trim = trim0;
		rem1.b[6] = 0;
		if (trim0){
			qctr = 1;
		}
		else{
			rem1.a <<= 1; ptr++;
			qctr = (int)rem1.b[6];
			rem1.b[6] = 0;
			if (!qctr){
				DES de = deSFcode(bin, 64);
				*ctr1 = de.sym;
			}
		}
	}
	else{
		qctr = 1;
		*trim = 0;
		*ctr1 = 0;
	}

	return qctr;
}