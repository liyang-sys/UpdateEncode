#include "decode.h"
#include "de_sub.h"

#define uchar unsigned char

void de_video_full(int H, int W, int *w, int *h, int **m0, double delta, int lenbits)
{
	extern uchar *bin;
	extern int ptr;
	extern int **nc;
	extern int len;
	extern float ***rePTVData;
	extern float ***reVideoData;

	ptr = 0;
	/*中间的32*/
	///*不需要harr部分*/
	for (int j = 35; j > 32; j--){
		de_subDC_noharr(j, w, h, m0, delta, lenbits);//不需要做harr的DC
		de_sub7_noharr(j, w, h, m0, delta, lenbits);//DC后的7个子带
		for (int i = 1; i < 8; i++)
			de_sub8_noharr(j, i, w, h, m0, delta, lenbits);//8个子带
	}
	/*L3 == == == == == == == == == == == == == == == == == ==*/
	de_coef3d_dc(32, w, h, m0, delta, lenbits);
	de_coed3d_7(32, w, h, m0, delta, lenbits);
	for (int i = 1; i < 8; i++)
		de_coed3d_8(32, i, w, h, m0, delta, lenbits);
	/*L4 == == == == == == == == == == == == == == == == == ==*/
	de_coef3dB_dc(32, w, h, m0, delta, lenbits);
	de_coef3dB_L3(H, W, 32, 0, 1, w, h, m0, delta, lenbits);
	de_coef3dB_L2_3(H, W, 32, 0, 2, w, h, m0, delta, lenbits);
	for (int i = 2; i < 8; i += 2)
		de_coef3dB_L2_4(H, W, 32, i, 0, w, h, m0, delta, lenbits);
	/*L5 == == == == == == == == == == == == == == == == == ==*/
	de_coef3d_dc5B_dc(32 + 16, w, h, m0, delta, 1, lenbits);
	de_coef3d_dc5B_L3(H, W, 32 + 16, 0, 1, w, h, m0, delta, lenbits);
	de_coef3d_dc5B_L2(H, W, 32 + 16, 0, 2, w, h, m0, delta, lenbits);
	de_coef3d_dc5B_L1(H, W, 32 + 16, 0, 4, w, h, m0, delta, lenbits);
	/*DC == == == == == == == == == == == == == == == == == ==*/
	de_coef3d_dc5B_dc(32, w, h, m0, delta, 0, lenbits);
	de_coef3d_dc5B_L3(H, W, 32, 0, 1, w, h, m0, delta, lenbits);
	de_coef3d_dc5B_L2(H, W, 32, 0, 2, w, h, m0, delta, lenbits);
	de_coef3d_dc5B_L1(H, W, 32, 0, 4, w, h, m0, delta, lenbits);

	/*前面的32*/
	///*不需要harr部分*/
	for (int j = 3; j > 0; j--){
		de_subDC_noharr(j, w, h, m0, delta, lenbits);//不需要做harr的DC
		de_sub7_noharr(j, w, h, m0, delta, lenbits);//DC后的7个子带
		for (int i = 1; i < 8; i++)
			de_sub8_noharr(j, i, w, h, m0, delta, lenbits);//8个子带
	}
	/*L3 == == == == == == == == == == == == == == == == == ==*/
	de_coef3d_dc(0, w, h, m0, delta, lenbits);
	de_coed3d_7(0, w, h, m0, delta, lenbits);
	for (int i = 1; i < 8; i++)
		de_coed3d_8(0, i, w, h, m0, delta, lenbits);
	/*L4 == == == == == == == == == == == == == == == == == ==*/
	de_coef3dB_dc(0, w, h, m0, delta, lenbits);
	de_coef3dB_L3(H, W, 0, 0, 1, w, h, m0, delta, lenbits);
	de_coef3dB_L2_3(H, W, 0, 0, 2, w, h, m0, delta, lenbits);
	for (int i = 2; i < 8; i += 2)
		de_coef3dB_L2_4(H, W, 0, i, 0, w, h, m0, delta, lenbits);
	/*L5 == == == == == == == == == == == == == == == == == ==*/
	de_coef3d_dc5B_dc(0 + 16, w, h, m0, delta, 1, lenbits);
	de_coef3d_dc5B_L3(H, W, 0 + 16, 0, 1, w, h, m0, delta, lenbits);
	de_coef3d_dc5B_L2(H, W, 0 + 16, 0, 2, w, h, m0, delta, lenbits);
	de_coef3d_dc5B_L1(H, W, 0 + 16, 0, 4, w, h, m0, delta, lenbits);
	/*DC == == == == == == == == == == == == == == == == == ==*/
	de_coef3d_dc5B_dc(0, w, h, m0, delta, 0, lenbits);
	de_coef3d_dc5B_L3(H, W, 0, 0, 1, w, h, m0, delta, lenbits);
	de_coef3d_dc5B_L2(H, W, 0, 0, 2, w, h, m0, delta, lenbits);
	de_coef3d_dc5B_L1(H, W, 0, 0, 4, w, h, m0, delta, lenbits);

	/*最后的32*/
	///*不需要harr部分*/
	for (int j = 67; j > 64; j--){
		de_subDC_noharr(j, w, h, m0, delta, lenbits);//不需要做harr的DC
		de_sub7_noharr(j, w, h, m0, delta, lenbits);//DC后的7个子带
		for (int i = 1; i < 8; i++)
			de_sub8_noharr(j, i, w, h, m0, delta, lenbits);//8个子带
	}
	/*L3 == == == == == == == == == == == == == == == == == ==*/
	de_coef3d_dc(64, w, h, m0, delta, lenbits);
	de_coed3d_7(64, w, h, m0, delta, lenbits);
	for (int i = 1; i < 8; i++)
		de_coed3d_8(64, i, w, h, m0, delta, lenbits);
	/*L4 == == == == == == == == == == == == == == == == == ==*/
	de_coef3dB_dc(64, w, h, m0, delta, lenbits);
	de_coef3dB_L3(H, W, 64, 0, 1, w, h, m0, delta, lenbits);
	de_coef3dB_L2_3(H, W, 64, 0, 2, w, h, m0, delta, lenbits);
	for (int i = 2; i < 8; i += 2)
		de_coef3dB_L2_4(H, W, 64, i, 0, w, h, m0, delta, lenbits);
	/*L5 == == == == == == == == == == == == == == == == == ==*/
	de_coef3d_dc5B_dc(64 + 16, w, h, m0, delta, 1, lenbits);
	de_coef3d_dc5B_L3(H, W, 64 + 16, 0, 1, w, h, m0, delta, lenbits);
	de_coef3d_dc5B_L2(H, W, 64 + 16, 0, 2, w, h, m0, delta, lenbits);
	de_coef3d_dc5B_L1(H, W, 64 + 16, 0, 4, w, h, m0, delta, lenbits);
	/*DC == == == == == == == == == == == == == == == == == ==*/
	de_coef3d_dc5B_dc(64, w, h, m0, delta, 0, lenbits);
	de_coef3d_dc5B_L3(H, W, 64, 0, 1, w, h, m0, delta, lenbits);
	de_coef3d_dc5B_L2(H, W, 64, 0, 2, w, h, m0, delta, lenbits);
	de_coef3d_dc5B_L1(H, W, 64, 0, 4, w, h, m0, delta, lenbits);


	return;
}