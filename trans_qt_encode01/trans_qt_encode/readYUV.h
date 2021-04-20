#ifndef _READYUV_H
#define _READYUV_H

unsigned char ***readYUV8(char* filesourse, int h, int w, int f1, int Nframe, int fmt);
unsigned short ***readYUV10(char* filesourse, int h, int w, int f1, int Nframe, int fmt);

#endif
//#ifndef _READYUV_H
//#define _READYUV_H
//
//unsigned char ***readYUV(char* filesourse, int h, int w, int f1, int Nframe, int fmt);
//
//#endif