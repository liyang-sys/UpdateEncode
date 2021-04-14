#include "encode.h"
#include "en_sub.h"
#include <stdio.h>

#define uchar unsigned char

void en_video_full(int H, int W, int *w, int *h, int **m0, double delta)
{
	extern uchar **sn;
	extern uchar *snbin;
	extern uchar *bin;
	extern int ptr;
	extern int **nc;
	extern int len;
	extern float ***PTVData;
	extern float ***VideoData;
	extern union Fabs f1;

	ptr = 0;
	/*中间的32*/
	for (int j = 35; j > 32; j--){ //L2a,L2b,L2c
		en_subDC_noharr(j, w, h, m0, delta);//不需要做harr的DC
		en_sub7_noharr(j, w, h, m0, delta);//DC后的7个子带
		//if (j == 34)
		//{
		//	int tempTest = 0;
		//}
		for (int i = 1; i < 8; i++)
		{
			//if (i == 7)
			//{
			//	int tempTest = 0;
			//}
			en_sub8_noharr(j, i, w, h, m0, delta);//8个子带
		}
			
	}
	/*L3 == == == == == == == == == == == == == == == == == ==*/
	en_coef3d_dc(32, w, h, m0, delta);
	en_coed3d_7(32, w, h, m0, delta);
	for (int i = 1; i < 8; i++)
		en_coed3d_8(32, i, w, h, m0, delta);
	/*L4 == == == == == == == == == == == == == == == == == ==*/
	en_coef3dB_dc(32, w, h, m0, delta);
	en_coef3dB_L3(H, W, 32, 0, 1, w, h, m0, delta);
	en_coef3dB_L2_3(H, W, 32, 0, 2, w, h, m0, delta);
	for (int i = 2; i < 8; i += 2)
		en_coef3dB_L2_4(H, W, 32, i, 0, w, h, m0, delta);
	/*L5 == == == == == == == == == == == == == == == == == ==*/
	en_coef3d_dc5B_dc(32 + 16, w, h, m0, delta, 1);
	en_coef3d_dc5B_L3(H, W, 32 + 16, 0, 1, w, h, m0, delta);
	en_coef3d_dc5B_L2(H, W, 32 + 16, 0, 2, w, h, m0, delta);
	en_coef3d_dc5B_L1(H, W, 32 + 16, 0, 4, w, h, m0, delta);
	/*DC == == == == == == == == == == == == == == == == == ==*/
	en_coef3d_dc5B_dc(32, w, h, m0, delta, 0);
	en_coef3d_dc5B_L3(H, W, 32, 0, 1, w, h, m0, delta);
	en_coef3d_dc5B_L2(H, W, 32, 0, 2, w, h, m0, delta);
	en_coef3d_dc5B_L1(H, W, 32, 0, 4, w, h, m0, delta);

	/*前面的32*/
	for (int j = 3; j > 0; j--){
		en_subDC_noharr(j, w, h, m0, delta);//不需要做harr的DC
		en_sub7_noharr(j, w, h, m0, delta);//DC后的7个子带
		for (int i = 1; i < 8; i++)
			en_sub8_noharr(j, i, w, h, m0, delta);//8个子带
	}
	/*L3 == == == == == == == == == == == == == == == == == ==*/
	en_coef3d_dc(0, w, h, m0, delta);
	en_coed3d_7(0, w, h, m0, delta);
	for (int i = 1; i < 8; i++)
	{
		//if (i == 5)
		//{
		//	int tempTest = 0;
		//}
		en_coed3d_8(0, i, w, h, m0, delta);
	}
		
	/*L4 == == == == == == == == == == == == == == == == == ==*/
	en_coef3dB_dc(0, w, h, m0, delta);
	en_coef3dB_L3(H, W, 0, 0, 1, w, h, m0, delta);
	en_coef3dB_L2_3(H, W, 0, 0, 2, w, h, m0, delta);
	for (int i = 2; i < 8; i += 2)
		en_coef3dB_L2_4(H, W, 0, i, 0, w, h, m0, delta);
	/*L5 == == == == == == == == == == == == == == == == == ==*/
	en_coef3d_dc5B_dc(0 + 16, w, h, m0, delta, 1);
	en_coef3d_dc5B_L3(H, W, 0 + 16, 0, 1, w, h, m0, delta);
	en_coef3d_dc5B_L2(H, W, 0 + 16, 0, 2, w, h, m0, delta);
	en_coef3d_dc5B_L1(H, W, 0 + 16, 0, 4, w, h, m0, delta);
	/*DC == == == == == == == == == == == == == == == == == ==*/
	en_coef3d_dc5B_dc(0, w, h, m0, delta, 0);
	en_coef3d_dc5B_L3(H, W, 0, 0, 1, w, h, m0, delta);
	en_coef3d_dc5B_L2(H, W, 0, 0, 2, w, h, m0, delta);
	en_coef3d_dc5B_L1(H, W, 0, 0, 4, w, h, m0, delta);

	/*后面的32*/
	for (int j = 67; j > 64; j--){
		en_subDC_noharr(j, w, h, m0, delta);//不需要做harr的DC
		en_sub7_noharr(j, w, h, m0, delta);//DC后的7个子带
		for (int i = 1; i < 8; i++)
			en_sub8_noharr(j, i, w, h, m0, delta);//8个子带
	}
	/*L3 == == == == == == == == == == == == == == == == == ==*/
	en_coef3d_dc(64, w, h, m0, delta);
	en_coed3d_7(64, w, h, m0, delta);
	for (int i = 1; i < 8; i++)
		en_coed3d_8(64, i, w, h, m0, delta);
	/*L4 == == == == == == == == == == == == == == == == == ==*/
	en_coef3dB_dc(64, w, h, m0, delta);
	en_coef3dB_L3(H, W, 64, 0, 1, w, h, m0, delta);
	en_coef3dB_L2_3(H, W, 64, 0, 2, w, h, m0, delta);
	for (int i = 2; i < 8; i += 2)
		en_coef3dB_L2_4(H, W, 64, i, 0, w, h, m0, delta);
	/*L5 == == == == == == == == == == == == == == == == == ==*/
	en_coef3d_dc5B_dc(64 + 16, w, h, m0, delta, 1);
	en_coef3d_dc5B_L3(H, W, 64 + 16, 0, 1, w, h, m0, delta);
	en_coef3d_dc5B_L2(H, W, 64 + 16, 0, 2, w, h, m0, delta);
	en_coef3d_dc5B_L1(H, W, 64 + 16, 0, 4, w, h, m0, delta);
	/*DC == == == == == == == == == == == == == == == == == ==*/
	en_coef3d_dc5B_dc(64, w, h, m0, delta, 0);
	en_coef3d_dc5B_L3(H, W, 64, 0, 1, w, h, m0, delta);
	en_coef3d_dc5B_L2(H, W, 64, 0, 2, w, h, m0, delta);
	en_coef3d_dc5B_L1(H, W, 64, 0, 4, w, h, m0, delta);

	return;
}