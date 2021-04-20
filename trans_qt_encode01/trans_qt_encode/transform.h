#ifndef _TRANSFORM_H
#define _TRANSFORM_H

float*** forward_transform(int N, int M, int depth); //N是宽 M是高
void inverse_transformation(int N, int M);

#endif
//#ifndef _TRANSFORM_H
//#define _TRANSFORM_H
//
//unsigned char*** forward_transform(int N, int M); //N是宽 M是高
//void inverse_transformation(int N, int M);
//
//#endif