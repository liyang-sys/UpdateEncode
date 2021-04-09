#ifndef _DE_SUB_H
#define _DE_SUB_H

float* de_DCinv(char* bin, int a, int lg, int* m, double delta);
void de_subDC_noharr(int a, int *w, int *h, int **m, double delta, int lenbits);
void de_sub7_noharr(int a, int *w, int *h, int **m, double delta, int lenbits);
void de_sub8_noharr(int a, int b, int *w, int *h, int **m, double delta, int lenbits);

void de_coef3d_dc(int a, int *w, int *h, int **m, double delta, int lenbits);
void de_coed3d_7(int a, int *w, int *h, int **m, double delta, int lenbits);
void de_coed3d_8(int a, int b, int *w, int *h, int **m, double delta, int lenbits);
void de_coef3dB_dc(int a, int *w, int *h, int **m, double delta, int lenbits);
void de_coef3dB_L3(int H, int W, int a, int b, int c, int *w, int *h, int **m, double delta, int lenbits);
void de_coef3dB_L2_3(int H, int W, int a, int b, int c, int *w, int *h, int **m, double delta, int lenbits);
void de_coef3dB_L2_4(int H, int W, int a, int b, int c, int *w, int *h, int **m, double delta, int lenbits);
void de_coef3d_dc5B_dc(int a, int *w, int *h, int **m, double delta, int AC, int lenbits);
void de_coef3d_dc5B_L3(int H, int W, int a, int b, int c, int *w, int *h, int **m, double delta, int lenbits);
void de_coef3d_dc5B_L2(int H, int W, int a, int b, int c, int *w, int *h, int **m, double delta, int lenbits);
void de_coef3d_dc5B_L1(int H, int W, int a, int b, int c, int *w, int *h, int **m, double delta, int lenbits);


#endif