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
  if (2*n > MEM_TMP_VEC) free(idx);  
}

/* LU decomposition */
void tf_lu(tMat* L, tMat* U, tMat* m, int* err)
{
  int e = 0, n, i, j, k;
  tmVal s, uii, *L_i;

  TM_ASSERT_ARGS(m && L && U && IS_UNIQUE3(m,L,U), e, end_lu);

  if (m->rows == m->cols) {
    n = m->rows;
    if (tm_relevant(L,n,n,&e) && tm_relevant(U,n,n,&e)) {
      /* prepare */
      tm_zeros(U);
      tm_eye(L); 
      /* fill */
      for (i = 0; i < n; i++) {
        L_i = L->data + i*n;
        /* U */
        for (j = i; j < n; j++) {
          s = *tm_at(m,i,j);
          for (k = 0; k < i; k++) 
            s -= L_i[k] * (*tm_at(U,k,j));
          *tm_at(U,i,j) = s;
        }
        uii = *tm_at(U,i,i);
        if (uii == 0.0) {
          e = TM_ERR_NO_SOLUTN;
          goto end_lu;
        }
        /* L */
        for (j = i+1; j < n; j++) {
          s = *tm_at(m,j,i);
          L_i = L->data + j*n;
          for (k = 0; k < i; k++) 
             s -= L_i[k] * (*tm_at(U,k,i));
          *tm_at(L,j,i) = s/uii;
        }
      }
    }
  } else
    e = TM_ERR_NOT_DEF;

end_lu:
  if (err) *err = e;
}

void householder(tMat* H, tMat* vec, int* err)
{
  int n, i, j; 
  tmVal s = 0, tmp, rowi;

  /* expected column-vector */
  if (vec->cols == 1) {
    n = vec->rows;
    if (H->cols == n && H->rows == n) {
      for (i = 0; i < n; i++) {
        tmp = *tm_at(vec,i,0);
        s += tmp*tmp;
      }
      if (s > 0) {
        s = 2.0 / s;
        for (i = 0; i < n; i++) {
          rowi = *tm_at(vec,i,0) * s;
          for (j = i; j < n; j++) {
            tmp = (i == j) ? 1 : 0;
            tmp -= rowi * (*tm_at(vec,j,0));
            *tm_at(H,i,j) = tmp;
            *tm_at(H,j,i) = tmp;
          }
        }
      } else 
        *err = TM_ERR_NO_SOLUTN;
    } else 
      *err = TM_ERR_NOT_COMPAT;
  } else
    *err = TM_ERR_NOT_VEC;
}

int qrdcmp(tMat* Qt, tMat* R)
{
  int i, j, k, sing=0;
  int nr = R->rows, nc = R->cols, ntot;
  tmVal scale, sigma, sum, tau, tmp, *ref;
  ntot = ((nr-1) < nc) ? (nr-1) : nc;

  for (k = 0; k < ntot; k++) {
    scale = 0;
    for (i=k; i < nc; i++) {
      tmp = fabs(*tm_at(R,i,k));
      scale = (tmp > scale) ? tmp : scale;
    }
    if (scale == 0.0) {
      sing = 1;
    } else {
      sum = 0;
      for (i = k; i < nr; i++) {
        ref = tm_at(R,i,k);
        *ref /= scale;
        sum += (*ref) * (*ref);
      } 
      sigma = sqrt(sum);
      sigma = (*tm_at(R,k,k) >= 0) ? sigma : (-sigma);
      ref = tm_at(R,k,k);
      *ref += sigma;
      tmp = sigma * (*ref);
      /* update R */
      for (j=k+1; j < nc; j++) {
        sum = 0.0;
        for (i = k; i < nr; i++) 
          sum += (*tm_at(R,i,k)) * (*tm_at(R,i,j));
        tau = sum/tmp;
        for (i = k; i < nr; i++)
          *tm_at(R,i,j) -= tau * (*tm_at(R,i,k));
      }
      /* update Q */
      for (j=0; j < nr; j++) {
        sum = 0.0;
        for (i = k; i < nr; i++) 
          sum += (*tm_at(R,i,k)) * (*tm_at(Qt,i,j));
        tau = sum/tmp;
        for (i = k; i < nr; i++)
          *tm_at(Qt,i,j) -= tau * (*tm_at(R,i,k));
      }

      *tm_at(R,k,k) = -scale*sigma;
    }
  }
  return sing;
}

#include <stdio.h>
// https://stackoverflow.com/questions/53489237/how-can-you-implement-householder-based-qr-decomposition-in-python
void tf_qr(tMat* Q, tMat* R, tMat* m, int* err)
{
  int e = 0, i, j;
  int nr = m->rows, nc = m->cols;
  tmVal tmp;
  
  TM_ASSERT_ARGS(m && Q && R && IS_UNIQUE3(m,Q,R), e, end_qr);

  if (tm_relevant(Q,nr,nr,&e) && tm_relevant(R,nr,nc,&e)) {
    for (i = 0; i < nr; i++) {
      for (j = 0; j < nc; j++)
        *tm_at(R,i,j) = *tm_at(m,i,j);
    }
    tm_eye(Q);

    qrdcmp(Q, R);
    /* transpose Q */
    for (i = 0; i < nr; i++) {
      for (j = i+1; j < nr; j++) {
        tmp = *tm_at(Q,i,j);
        *tm_at(Q,i,j) = *tm_at(Q,j,i);
        *tm_at(Q,j,i) = tmp;
      }
    }
    /* clear R */
    for (j = 0; j < nc; j++) {
      for (i = j+1; i < nr; i++) 
        *tm_at(R,i,j) = 0;
    }
  }

end_qr:
  if (err) *err = e;
}
