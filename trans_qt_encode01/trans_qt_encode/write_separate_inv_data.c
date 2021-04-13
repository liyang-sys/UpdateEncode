#define _CRT_SECURE_NO_WARNINGS
#include "all.h"
#include "encoding.h"
#include "parameter_setting.h"
void write_separate_inv_data(Uint32_Dat rw1, Uint32_Dat rk1 )
{
	//printf("cf0µÿ÷∑ = %p\n", cf0);
	FILE*fp  = fopen("cRw1.txt", "wb");
	fwrite(rw1.dat, sizeof(unsigned  int), rw1.len, fp);
	fclose(fp);

	fp = fopen("cRk1.txt", "wb");
	fwrite(rk1.dat, sizeof(unsigned  int), rk1.len, fp);
	fclose(fp);
}