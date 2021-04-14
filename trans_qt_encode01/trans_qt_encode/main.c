#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "encode.h"
#include "en_sub.h"
#include "decode.h"
#include "transform.h"
#include <math.h>
#include "quant_rest.h"
#include "encoding.h"

unsigned char **sn;
unsigned char *snbin;
unsigned char *bin;
int ptr = 0;
int **nc;
int len;
float ***VideoData;
float ***reVideoData;
float ***PTVData;
float ***rePTVData;
double *absA2;
double *absA;
double *T;
double *T2;
union Fabs f1;

void memoryAllocation(int H, int W, int m, int n)
{
	if ((VideoData = (float***)malloc(sizeof(float**)* 96)) == NULL){
		printf("创建VideoData数组第一维失败！\n");
		exit(1);
	}
	for (int i = 0; i < 96; i++){
		if ((VideoData[i] = (float**)malloc(sizeof(float*)* H)) == NULL){
			printf("创建VideoData数组第二维失败！\n");
			exit(1);
		}
		for (int j = 0; j < H; j++){
			if ((VideoData[i][j] = (float*)malloc(sizeof(float)* W)) == NULL){
				printf("创建VideoData数组第三维失败！\n");
				exit(1);
			}
		}
	}
	if ((reVideoData = (float***)malloc(sizeof(float**)* 96)) == NULL){
		printf("创建VideoData数组第一维失败！\n");
		exit(1);
	}
	for (int i = 0; i < 96; i++){
		if ((reVideoData[i] = (float**)malloc(sizeof(float*)* H)) == NULL){
			printf("创建VideoData数组第二维失败！\n");
			exit(1);
		}
		for (int j = 0; j < H; j++){
			if ((reVideoData[i][j] = (float*)malloc(sizeof(float)* W)) == NULL){
				printf("创建VideoData数组第三维失败！\n");
				exit(1);
			}
		}
	}
	if ((PTVData = (float***)malloc(sizeof(float**)* 96)) == NULL){
		printf("创建PTVData数组第一维失败！\n");
		exit(1);
	}
	for (int i = 0; i < 96; i++){
		if ((PTVData[i] = (float**)malloc(sizeof(float*)* m)) == NULL){
			printf("创建PTVData数组第二维失败！\n");
			exit(1);
		}
		for (int j = 0; j < m; j++){
			if ((PTVData[i][j] = (float*)malloc(sizeof(float)* n)) == NULL){
				printf("创建PTVData数组第三维失败！\n");
				exit(1);
			}
		}
	}
	if ((rePTVData = (float***)malloc(sizeof(float**)* 96)) == NULL){
		printf("创建PTVData数组第一维失败！\n");
		exit(1);
	}
	for (int i = 0; i < 96; i++){
		if ((rePTVData[i] = (float**)malloc(sizeof(float*)* m)) == NULL){
			printf("创建PTVData数组第二维失败！\n");
			exit(1);
		}
		for (int j = 0; j < m; j++){
			if ((rePTVData[i][j] = (float*)malloc(sizeof(float)* n)) == NULL){
				printf("创建PTVData数组第三维失败！\n");
				exit(1);
			}
		}
	}
	if ((absA2 = (double*)malloc(sizeof(double)*len)) == NULL){
		printf("创建absA2数组失败！\n");
		exit(1);
	}
	if ((absA = (double*)malloc(sizeof(double)*len)) == NULL){
		printf("创建absA2数组失败！\n");
		exit(1);
	}
	if ((T = (double*)malloc(sizeof(double)* 131072)) == NULL){
		printf("创建T数组失败！\n");
		exit(1);
	}
	if ((T2 = (double*)malloc(sizeof(double)* 131073)) == NULL){
		printf("创建T2数组失败！\n");
		exit(1);
	}
	if ((f1.qcf = (float**)malloc(sizeof(float*)* 8)) == NULL){
		printf("创建qcf数组第一维失败！\n");
		exit(1);
	}
	for (int i = 0; i < 8; i++){
		if ((f1.qcf[i] = (float*)malloc(sizeof(float)*len)) == NULL){
			printf("创建qcf数组第二维失败！\n");
			exit(1);
		}
	}
	if ((sn = (unsigned char**)calloc(8, sizeof(char*))) == NULL){
		printf("创建sn数组失败！\n");
		exit(1);
	}
	for (int i = 0; i < 8; i++){
		if ((sn[i] = (unsigned char*)calloc(len, sizeof(char))) == NULL){
			printf("创建sn数组第二维失败！\n");
			exit(1);
		}
	}
	if ((snbin = (unsigned char*)calloc((len >> 3) + 4, sizeof(char))) == NULL){
		printf("创建snbin数组失败！\n");
		exit(1);
	}
	if ((bin = (unsigned char*)calloc(H*W*120, sizeof(unsigned char))) == NULL){
		printf("创建bin数组失败！\n");
		exit(1);
	}

	if ((nc = (int**)calloc(8, sizeof(int*))) == NULL){
		printf("创建nc数组第一维失败！\n");
		exit(1);
	}
	for (int i = 0; i < 8; i++){
		if ((nc[i] = (int*)calloc(len * 10, sizeof(int))) == NULL){
			printf("创建nc数组第二维失败！\n");
			exit(1);
		}
	}
}
//int main()
//{
//	test_en_sub3d_sub2();
//}
int main()
{
	int W, H, n, m;
	int *m0[42], w[7], h[7];
	double T0, delta;
	FILE *fp3;
	errno_t err;
	printf("Please input the width\n");   //W列
	scanf_s("%d", &W);
	while (W % 8 != 0)
	{
		printf("Your input is illegal,Please try again\n");
		scanf_s("%d", &W);
	}

	printf("Please input the height\n");   //H行
	scanf_s("%d", &H);
	while (H % 8 != 0)
	{
		printf("Your input is illegal,Please try again\n");
		scanf_s("%d", &H);
	}

	printf("请输入delta：\n");
	scanf_s("%lf", &delta, sizeof(double));
	m = H / 8;
	n = W / 8;
	len = n*m * 8*4;
	printf("准备工作开始:\n");
	//内存分配
	memoryAllocation(H, W, m, n);
	//读取相关数据
	printf("开始数据读取！\n");

	if ((err = fopen_s(&fp3, "Tcoef.dat", "rb")) != 0){
		printf("打开T数据失败！\n");
		exit(1);
	}
	fseek(fp3, 0, SEEK_SET);
	fread(T, sizeof(double), 131071, fp3);
	fclose(fp3);
	printf("数据读取完成！\n");

	T0 = 0.55*delta;
	T2[1] = T0*T0; T2[0] = 0; //这里的T和T2的长度都增加了一个，是因为和matlab匹配，从1开始索引，T[0]相当于为空
	for (int i = 131070; i >= 0; i--){
		T[i + 1] = T[i] * delta;
		T2[i + 2] = T[i+1] * T[i+1];
	}
	printf("准备工作完成！\n\n");
	//.....正变换开始：.....
	clock_t st2, end2;
	unsigned char*** imageData = forward_transform(W, H); //变换

	//测试部分
	st2 = clock();
	GetM(W, H, w, h, m0);
	printf("\nHilbert路径获取完毕！\n\n");
	printf(".....量化编码开始：.....\n");
	
	//=====================================编码=============================================
	en_video_full(H, W, w, h, m0, delta); //编码

	//if ((err = fopen_s(&fp3, "E:\\程序代码\\量化\\2021-01-15\\test_en_sub3d_sub2_sub\\cbin.dat", "wb")) == 0){
	//	fwrite(bin, sizeof(uchar), (ptr>>3)+1, fp3);/*4*/
	//	fclose(fp3);
	//}

	writeBinToFile(bin, ptr);
	end2 = clock();
	double duration2 = (double)(end2 - st2) / CLOCKS_PER_SEC;
	printf("\n\n时间 %f seconds\n", duration2);
	printf("ptr=%d\n\n\n", ptr);
	printf("bitrate=%f\n", (double)ptr / (H*W * 96));
	printf("编码结束-----------------------------\n解码开始......................\n");
	//=====================================解码=============================================
	int lenbits = ptr;
	st2 = clock();

	de_video_full(H, W, w, h, m0, delta, lenbits); //解码

	end2 = clock();
	duration2 = (double)(end2 - st2) / CLOCKS_PER_SEC;
	printf("\n\n时间 %f seconds\n", duration2);

	//rearrange
#if 0
	FILE *fp2;
	if ((err = fopen_s(&fp2, "E:\\程序代码\\量化\\2021-01-15\\test_en_sub3d_sub2_sub\\PTVoutput_rearrange.dat", "wb")) == 0){
		for (int i = 0; i < 96; i++)
		{
			for (int j = 0; j < m; j++)
			{
				fwrite(&PTVData3[i][j][0], sizeof(float), n, fp2);
			}
		}
		fclose(fp2);
	}
#endif
	/*--------------------------------------------------------------------------------*/
	//.....逆变换开始：.....
	inverse_transformation(W, H); //逆变换

	for (int i = 0; i < 8; i++){
		for (int j = 0; j < 8; j++)
			printf("%f ", reVideoData[0][i][j]);
		printf("\n");
	}
	printf("\n\n");
#if 0
	//保存Video数据和PTV数据
	if ((err = fopen_s(&fp3, "E:\\程序代码\\量化\\2021-01-15\\test_en_sub3d_sub2_sub\\reVideo.dat", "wb")) != 0){
		printf("打开数据失败！\n");
		system("pause");
		return 0;
	}
	for (int i = 0; i < 96; i++)
	{
		for (int j = 0; j<1080; j++)
		{
			fwrite(&reVideoData[i][j][0], sizeof(float), 1920, fp3);
		}
	}
	fclose(fp3);

	//rearrange==================================================================
	float ***PTVData2, ***PTVData3;
	int i, j, k;
	PTVData2 = (float ***)malloc(sizeof(float **)* 96);
	for (i = 0; i < 96; i++)
	{
		PTVData2[i] = (float**)malloc(sizeof(float *)* 135);
		for (j = 0; j < 135; j++)
		{
			PTVData2[i][j] = (float*)malloc(sizeof(float)* 240);
		}
	}
	PTVData3 = (float ***)malloc(sizeof(float **)* 96);
	for (i = 0; i < 96; i++)
	{
		PTVData3[i] = (float**)malloc(sizeof(float *)* 135);
		for (j = 0; j < 135; j++)
		{
			PTVData3[i][j] = (float*)malloc(sizeof(float)* 240);
		}
	}
	for (i = 0; i<96; i++)
	{
		for (j = 0; j<m; j++)
		{
			for (k = 0; k<n; k++)
			{
				PTVData2[i][j][k] = rePTVData[i][j][k];
			}
		}
	}
	int m2 = m;
	int n2 = n;
	if (m % 2 != 0)
	{
		m2 += 1;
	}
	if (n % 2 != 0)
	{
		n2 += 1;
	}
	m2 /= 2;
	n2 /= 2;
	for (i = 0; i<96; i++)
	{
		for (j = 0; j<m; j++)
		{
			for (k = 0; k<n; k++)
			{
				PTVData2[i][j / 2 + (j % 2)*m2][k / 2 + (k % 2)*n2] = rePTVData[i][j][k];
			}
		}
	}
	if (m2 % 2 != 0)
	{
		m2 += 1;
	}
	if (n2 % 2 != 0)
	{
		n2 += 1;
	}
	m2 /= 2;
	n2 /= 2;
	for (i = 0; i<96; i++)
	{
		for (j = 0; j<m; j += 2)
		{
			for (k = 0; k<n; k += 2)
			{
				PTVData2[i][j / 4 + ((j / 2) % 2)*m2][k / 4 + ((k / 2) % 2)*n2] = rePTVData[i][j][k];
			}
		}
	}
	if (m2 % 2 != 0)
	{
		m2 += 1;
	}
	if (n2 % 2 != 0)
	{
		n2 += 1;
	}
	m2 /= 2;
	n2 /= 2;
	for (i = 0; i<96; i += 8)
	{
		for (j = 0; j<m; j += 4)
		{
			for (k = 0; k<n; k += 4)
			{
				PTVData2[i][j / 8 + (j / 4) % 2 * m2][k / 8 + (k / 4) % 2 * n2] = rePTVData[i][j][k];
			}
		}
	}
	for (i = 0; i<96; i++)
	{
		for (j = 0; j<m; j++)
		{
			for (k = 0; k<n; k++)
			{
				PTVData3[i / 4 + (i % 4) * 24][j][k] = PTVData2[i][j][k];
			}
		}
	}
	for (i = 0; i<96; i += 4)
	{
		for (j = 0; j<m; j++)
		{
			for (k = 0; k<n; k++)
			{
				PTVData3[i / 8 + (i / 4) % 2 * 12][j][k] = PTVData2[i][j][k];
			}
		}
	}
	for (i = 0; i<96; i += 8)
	{
		for (j = 0; j<m; j++)
		{
			for (k = 0; k<n; k++)
			{
				PTVData3[i / 16 + (i / 8) % 2 * 6][j][k] = PTVData2[i][j][k];
			}
		}
	}
	for (i = 0; i<96; i += 16)
	{
		for (j = 0; j<m; j++)
		{
			for (k = 0; k<n; k++)
			{
				PTVData3[i / 32 + (i / 16) % 2 * 3][j][k] = PTVData2[i][j][k];
			}
		}
	}

	if ((err = fopen_s(&fp3, "E:\\程序代码\\量化\\2021-01-15\\test_en_sub3d_sub2_sub\\rePTV.dat", "wb")) != 0){
		printf("打开数据失败！\n");
		system("pause");
		return 0;
	}
	for (i = 0; i<96; i++)
	{
		for (j = 0; j<m; j++)
		{
			fwrite(&PTVData3[i][j][0], sizeof(float), n, fp3);
		}
	}
	fclose(fp3);
#endif
#if 1
	/*-----------------------------------计算误差---------------------------------------------*/
	float *snr;
	float Npix = (float)H*W;
	float e, sume = 0.0;
	float sumsnr = 0, minsnr = 100.0;
	snr = (float*)calloc(96, sizeof(float));
	for (int i = 0; i < 96; i++){
		sume = 0.0;
		for (int j = 0; j < 1080; j++){
			for (int k = 0; k < 1920; k++){
				e = (float)imageData[i][j][k] - reVideoData[i][j][k];
				e *= e;
				sume += e;
			}
		}
		double er = sqrt((double)(sume / Npix));
		snr[i] = 20.0*log10((double)(255.0 / er));
		//printf("%f\n", snr);
		sumsnr += snr[i];
		if (snr[i] < minsnr)
			minsnr = snr[i];
	}
	double meansnr = sumsnr / 96.0;

	//if ((err = fopen_s(&fp3, "E:\\程序代码\\量化\\2021-01-15\\test_en_sub3d_sub2_sub\\snr_c5.dat", "wb")) != 0){
	//	printf("打开数据失败！\n");
	//	system("pause");
	//	return 0;
	//}
	//fwrite(snr, sizeof(float), 96, fp3);
	//fclose(fp3);

	printf("mean(snr)=%f , min(snr)=%f\n", meansnr, minsnr);

#endif

	printf("程序结束！\n\n");
	system("pause");
	return 0;
}

