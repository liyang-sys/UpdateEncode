#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>
#include <math.h>
#include "all.h"
#include "encoding.h"
#include "parameter_setting.h"

int ptr;
uchar *bin;
void test_en_sub3d_sub2()
{
	ptr = 0;
	bin = calloc(4096000, sizeof(uchar));
	union data
	{
		unsigned short int a;
		uchar b[4];
	} rem;
	int snLen = 877;
	unsigned char *bigZ = (unsigned char *)calloc(snLen, sizeof(unsigned char));
	FILE *fp = fopen("SN.txt", "rb");
	fread(bigZ, sizeof(unsigned char), snLen / 8 + 1, fp);
	fclose(fp);
	//int index = 0;
	//int indexZ = 0;
	Uint8_Dat sn;
	sn.len = snLen;
	sn.dat = bigZ;
	//int sumZ = 0;
	//for (int i = 0; i < snLen; i++)
	//{
	//	if (bigZ[i] == 1)
	//	{
	//		rem.a = 1;
	//		rem.a = rem.a << (15 - index);
	//		sn.dat[indexZ >> 3] |= rem.b[1];
	//		sn.dat[(indexZ >> 3) + 1] |= rem.b[0];
	//		sumZ++;
	//	}
	//	index++;
	//	indexZ++;
	//	index &= 7;
	//	printf("%d ", bigZ[i]);
	//}
	//printf("\n------------------------------------------------------------------------------------------\n");
	///从matlab读取runs数组
	//Int32_Dat cf0;
	//cf0.len = 259200;
	//cf0.dat = (int*)calloc(cf0.len, sizeof(int));

	//fp = fopen("wunsDAYUlenf1C2N.txt","rb");

	//    printf("c2n = %p",c2n.dat);

	//for (int i = 0; i < cf0.len; i++)
	//{
	//	printf("%d ", cf0.dat[i]);
	//}
	fp = fopen("CF0Float.txt", "rb");
	int lg = 8040;
	float *cf0 = (float *)calloc(lg, sizeof(float));
	fread(cf0, sizeof(float), lg, fp);
	fclose(fp);

	int maxcf0 = 42;
	en_sub3d_sub2(cf0, &sn, lg, maxcf0);
	printf("\n编码结束后ptr = %d\n", ptr);
	fp = fopen("encodeZThd4.txt", "wb");
	fwrite(bin, sizeof(unsigned char), (ptr / 8) + 1, fp);
	fclose(fp);
	int lenbit = ptr;
	ptr = 0;
	DE_S_SUB res = de_sub3d_sub2(bin, lg, 1, lenbit);

	//fp = fopen("CF0SN.txt", "rb");
	////int lg = 259200;
	//int *cf1 = (int *)calloc(lg, sizeof(int));
	//fread(cf1, sizeof(int), lg, fp);
	//fclose(fp);
	//for (int i = 0; i < res.cf.len; i++)
	//{
	//	if (res.cf.dat[i] != cf1[i])
	//	{
	//		printf("解码不正确 idx = %d\n", i);
	//	}
	//}
	//DEC dec = de_z0_r0(bin, 8192);
	fp = fopen("cQcf.txt", "wb");
	fwrite(res.cf.dat, sizeof(int), res.cf.len, fp);
	fclose(fp);

	//DE_S_S_SUB ret = de_sub3d_sub2_sub(bin, 1,8040 , lenbit);
	int x = 0;
}