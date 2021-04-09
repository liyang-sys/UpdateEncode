#ifndef _QUANT_INV_H
#define _QUANT_INV_H

#define uchar unsigned char
#define uint unsigned int

float* rst_sub(int* qcf, double delta, int lg, uchar qnt, int qctr, uchar trim, uint ctr1);
float* rstTHDctr(int *qcf, char *sgn, double T, double delta, int* nc, int maxcf, int lg);
float* rstEVENctr(int* qcf, char *sgn, double delta, int* nc, int maxcf, double ctr1, int lg);
void scan_inv(int a, int b, int c, int lg, int *m, int ptv, int harr, int dc, float *qcf);
void scan_inv_L2(int H, int W, int a, int b, int c, int lg, int *m, int harr, float *qcf);
int find_qctr(uchar *bin, uchar *qnt, uchar *trim, uint *ctr1);

#endif