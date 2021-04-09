#ifndef _QUANT_REST_H
#define _QUANT_REST_H

struct QuantTHD{
	int *runs;
	int lg;
};
typedef struct QuantTHD QT;

#define uint unsigned int
#define uchar unsigned char

int quanEVEN7(int *max, int a, int lg, int ptr, int harr, int *m, double delta0);
int quanEVEN8(int *max, int a, int b, int lg, int ptr, int harr, int *m, double delta0);
int quanEVEN2(int len, int ptr, float *absA, double delta);
QT quanTHD(uchar *sgn, int len, int ptr, float *absA, double delta);
int quanTHD2(float *absA, double delta);
int quanEVEN_DC(int a, int b, int c, int lg, int ptr, int ptv, int harr, int *m, double delta);

void quant3(int *max, int H, int W, int a, int b, int c, int lg, int harr, int mul, int *m, int overflow, double delta0);
void quant_L2_3(int *max, int H, int W, int a, int b, int c, int lg, int harr, int *m, double delta0);
void quant_L2_4(int *max, int H, int W, int a, int b, int c, int lg, int harr, int *m, double delta0);


#endif