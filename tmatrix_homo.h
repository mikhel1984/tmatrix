#ifndef T_MATRIX_HOMO_H
#define T_MATRIX_HOMO_H

#include "tmatrix.h"

#define TM_HOMO_SIZE 16

#define h_Tx(dst,x,err)  h_Txyz(dst,x,0,0,err)
#define h_Ty(dst,y,err)  h_Txyz(dst,0,y,0,err)
#define h_Tz(dst,z,err)  h_Txyz(dst,0,0,z,err)

#define h_new(err)       tm_new(4,4,err)
#define h_static(d,err)  tm_static(4,4,d,err)

int h_Txyz(tMat* dst, tmVal x, tmVal y, tmVal z, int *err);

int h_Rz(tMat* dst, tmVal a, int *err);

int h_Ry(tMat* dst, tmVal b, int *err);

int h_Rx(tMat* dst, tmVal v, int *err);

int h_DH(tMat* dst, tmVal a, tmVal alpha, tmVal d, tmVal theta, int *err);

int h_mul(tMat *dst, tMat *m, int* err);

int h_inv(tMat *dst, int *err);

int h_T(tMat *dst, int *err);

int h_eye(tMat *dst, int *err);

#endif 
