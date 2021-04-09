#ifndef _EN_SUB_H
#define _EN_SUB_H

union Fabs{
	float **qcf;
	unsigned int **temp;
};

//*整理后的总的函数*//
void en_subDC_noharr(int a, int *w, int *h, int **m, double delta);
void en_sub7_noharr(int a, int *w, int *h, int **m, double delta);
void en_sub8_noharr(int a, int b, int *w, int *h, int **m, double delta);

void en_coef3d_dc(int a, int *w, int *h, int **m, double delta);
void en_coed3d_7(int a, int *w, int *h, int **m, double delta);
void en_coed3d_8(int a, int b, int *w, int *h, int **m, double delta);

void en_coef3dB_dc(int a, int *w, int *h, int **m, double delta);
void en_coef3dB_L3(int H, int W, int a, int b, int c, int *w, int *h, int **m, double delta);
void en_coef3dB_L2_3(int H, int W, int a, int b, int c, int *w, int *h, int **m, double delta);
void en_coef3dB_L2_4(int H, int W, int a, int b, int c, int *w, int *h, int **m, double delta);

void en_coef3d_dc5B_dc(int a, int *w, int *h, int **m, double delta, int AC);
void en_coef3d_dc5B_L3(int H, int W, int a, int b, int c, int *w, int *h, int **m, double delta);
void en_coef3d_dc5B_L2(int H, int W, int a, int b, int c, int *w, int *h, int **m, double delta);
void en_coef3d_dc5B_L1(int H, int W, int a, int b, int c, int *w, int *h, int **m, double delta);

//*GetM*//获取路径 
int GetM(int W, int H, int *w, int *h, int **m);
int handleSn(unsigned char *sgn, int idx, int lg, int qctr, int *maxcf, float *fabsA);


#endif