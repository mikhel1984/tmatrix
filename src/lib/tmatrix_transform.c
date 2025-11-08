/**
 * @file tmatrix_transform.c
 * @author Stanislav Mikhel
 * @date 2025
 * @brief Matrix decompositions.
 */
#include <stdlib.h>
#include <math.h>
#include "tmatrix.h"
#include "tmatrix_priv.h"

#define MEM_TMP_VEC 8

/* Cholesky decomposition */
void tf_chol(tMat* dst, tMat* m, int* err)
{
  int e = 0, i, j, k, n;
  tmVal s = 0, *row_i, *row_j;

  TM_ASSERT_ARGS(m && dst && m != dst, e, end_chol);

  if (m->rows == m->cols) {
    n = m->rows;
    if (tm_relevant(dst,n,n,&e)) {
      /* clear */
      tm_zeros(dst);
      /* find lower left part */
      for (i = 0; i < n; i++) {
        row_i = dst->data + i*n;
        for (j = 0; j <= i; j++) {
          row_j = dst->data + j*n;
          s = *tm_at(m, i, j);
          for (k = 0; k < j; k++) 
            s -= row_i[k] * row_j[k];
          if (j < i) {
            row_i[j] = s / row_j[j];
          } else {
            if (s < 0) {
              e = TM_ERR_NOT_POS_DEF;
              goto end_chol;
            }
            row_i[j] = sqrt(s);
          }
        }
      }
    }
  } else 
    e = TM_ERR_NOT_DEF;
  
end_chol:
  if(err) *err = e;
}

/* LUP decomposition */
void tf_lup(tMat* L, tMat* U, tMat* P, tMat* m, int* err)
{
  int e = 0, n, i, j;
  tmVal tmp = 0;
  int arr[MEM_TMP_VEC] = {0}, *idx = NULL, *acc = NULL, swap;

  TM_ASSERT_ARGS(m && L && U && P && IS_UNIQUE4(m,L,U,P), e, end_lup);

  if (m->rows == m->cols) {
    n = m->rows;
    if (tm_relevant(L,n,n,&e) && tm_relevant(U,n,n,&e) && tm_relevant(P,n,n,&e)) {
      /* allocate memory */
      idx = (2*n <= MEM_TMP_VEC) ? arr : (int*) malloc(2*n * sizeof(int));
      if (!idx) {
        e = TM_ERR_NO_MEMORY;
        goto end_lup;
      }
      acc = idx + n;  /* store order in the end of the list */
      /* copy source matrix */
      for (i = 0; i < n; i++) {
        acc[i] = i;
        for (j = 0; j < n; j++)
          *tm_at(U, i, j) = *tm_at(m, i, j);
      }
      /* find decomposition */
      ludcmp(U,idx,&tmp,&e); 
      if (e) goto end_lup;
      for (i = 0; i < n; i++) {
        /* fill L and U */
        for (j = 0; j < n; j++) {
          if (j > i) {
            *tm_at(L, i, j) = 0;
          } else if (j == i) {
            *tm_at(L, i, j) = 1;
          } else {
            *tm_at(L, i, j) = *tm_at(U, i, j);
            *tm_at(U, i, j) = 0;
          }
        }
        /* row permutations */
        if (idx[i] != i) {
          swap = acc[i];
          acc[i] = acc[idx[i]]; 
          acc[idx[i]] = swap;
        }      
      }
      tm_zeros(P);
      for (i = 0; i < n; i++) {
        *tm_at(P, i, acc[i]) = 1;
      }
    }
  } else
    e = TM_ERR_NOT_DEF;

end_lup:
  if (err) *err = e;
  if (n > MEM_TMP_VEC) free(idx);  
}
