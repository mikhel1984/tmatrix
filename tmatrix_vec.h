#ifndef T_MATRIX_VECTOR_H
#define T_MATRIX_VECTOR_H

#include "tmatrix.h"

#define vec_new(n,err)        tm_new(n,1,err)
#define vec_static(n,dat,err) tm_static(n,1,dat,err)

tmSize vec_len(tMat *m);

tmVal vec_get(tMat *m, tmSize k, int *err);

void vec_set(tMat *m, tmSize k, tmVal v, int *err);

tmVal vec_dot(tMat *a, tMat *b, int *err);

int vec_cross(tMat *res, tMat *a, tMat *b, int *err);



#endif /* T_MATRIX_VECTOR_H */