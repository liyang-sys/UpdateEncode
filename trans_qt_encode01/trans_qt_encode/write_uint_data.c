#define _CRT_SECURE_NO_WARNINGS
#include "all.h"
#include "encoding.h"
#include "parameter_setting.h"
void write_uint_data(unsigned int *data, int len, char*name)
{
	FILE* fp = fopen(name, "wb");
	fwrite(data, sizeof(unsigned int), len, fp);
	fclose(fp);

}