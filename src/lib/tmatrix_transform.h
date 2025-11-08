
#ifndef T_MATRIX_TRANSFORM_H
#define T_MATRIX_TRANSFORM_H

#include "tmatrix.h"

int tf_chol(tMat* dst, tMat* m, int* err);

void tf_lup(tMat* L, tMat* U, tMat* P, tMat* m, int* err);

#endif 
