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
      if (!IS_PRIM(dst)) {
        e = TM_ERR_NOT_MAIN;
        goto end_chol;
      } 
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
  tmVal tmp = 0, *row;
  int arr[MEM_TMP_VEC] = {0}, *idx = NULL, *acc = NULL, swap;

  TM_ASSERT_ARGS(m && L && U && P && IS_UNIQUE4(m,L,U,P), e, end_lup);

  if (m->rows == m->cols) {
    n = m->rows;
    if (tm_relevant(L,n,n,&e) && tm_relevant(U,n,n,&e) && tm_relevant(P,n,n,&e)) {
      if (!IS_PRIM(L) || !IS_PRIM(U) || !IS_PRIM(P)) {
        e = TM_ERR_NOT_MAIN;
        goto end_lup;
      } 
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
        row = U->data + i*n;
        for (j = 0; j < n; j++)
          row[j] = *tm_at(m, i, j);
      }
      /* find decomposition */
      ludcmp(U,idx,&tmp,&e); 
      if (e) goto end_lup;
      for (i = 0; i < n; i++) {
        /* fill L and U */
        row = L->data + i*n;
        for (j = 0; j < n; j++) {
          if (j > i) {
            row[j] = 0;
          } else if (j == i) {
            row[j] = 1;
          } else {
            row[j] = *tm_at(U, i, j);
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
      if (!IS_PRIM(L) || !IS_PRIM(U)) {
        e = TM_ERR_NOT_MAIN;
        goto end_lu;
      }
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

/**
 * @brief Prepare QR decomposition.
 * @param Qt matrix to store transposed Q value.
 * @param R source matrix, save upper triangular result here.
 * @return 1 when matrix is singular.
 */
int qrdcmp(tMat* Qt, tMat* R)
{
  int i, j, k, singular=0;
  int nr = R->rows, nc = R->cols, ntot;
  tmVal scale, sigma, sum, tau, tmp, *ref;
  ntot = ((nr-1) < nc) ? (nr-1) : nc;

  for (k = 0; k < ntot; k++) {
    scale = 0;
    ref = tm_at(R,k,k);
    for (i=k; i < nc; i++) {
      tmp = fabs(ref[(i-k)*nc]);  // tm_at(R,i,k)
      scale = (tmp > scale) ? tmp : scale;
    }
    if (scale == 0.0) {
      singular = 1;
    } else {
      sum = 0;
      /* ref equal to R[k][k] */
      for (i = k; i < nr; i++) {
        j = (i-k)*nc;      // tm_at(R,i,k)
        ref[j] /= scale;
        sum += ref[j] * ref[j];
      } 
      sigma = sqrt(sum);
      sigma = (*ref >= 0) ? sigma : (-sigma);
      *ref += sigma;
      tmp = sigma * (*ref);
      /* update R */
      for (j=k+1; j < nc; j++) {
        sum = 0.0;
        for (i = k; i < nr; i++) {
          ref = R->data + i*nc;  // row i
          sum += ref[k] * ref[j];
        }
        tau = sum/tmp;
        for (i = k; i < nr; i++) {
          ref = R->data + i*nc;  // row i
          ref[j] -= tau * ref[k];
        }
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
  return singular;
}

/* QR decomposition */
void tf_qr(tMat* Q, tMat* R, tMat* m, int* err)
{
  int e = 0, i, j;
  int nr = m->rows, nc = m->cols;
  tmVal tmp, *row_i, *ref;
  
  TM_ASSERT_ARGS(m && Q && R && IS_UNIQUE3(m,Q,R), e, end_qr);

  if (tm_relevant(Q,nr,nr,&e) && tm_relevant(R,nr,nc,&e)) {
    if (!IS_PRIM(Q) || !IS_PRIM(R)) {
      e = TM_ERR_NOT_MAIN;
      goto end_qr;
    } 
    /* copy source matrix */
    for (i = 0; i < nr; i++) {
      row_i = R->data + i*nc;
      for (j = 0; j < nc; j++)
        row_i[j] = *tm_at(m,i,j);
    }
    tm_eye(Q);

    qrdcmp(Q, R);
    /* transpose Q */
    for (i = 0; i < nr; i++) {
      row_i = Q->data + i*nr;
      for (j = i+1; j < nr; j++) {
        tmp = row_i[j];
        ref = tm_at(Q,j,i);
        row_i[j] = *ref;
        *ref = tmp;
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

tmVal pythag(tmVal a, tmVal b)
{
  a = fabs(a);
  b = fabs(b);
  if (a > b) {
    b /= a;
    return a * sqrt(1.0 + b*b);
  } else if (b == 0.0) {
    return 0.0;
  } else {
    /* a <= b */
    a /= b;
    return b * sqrt(1.0 + a*a);
  }
}

#include <stdio.h>
void svdcmp(tMat* U, tMat* V, tMat* W, int* e)
{
  int i, j, k, l, its, flag, nm, jj;
  tmVal g = 0, scale = 0, anorm = 0, s, c, f, h;
  tmVal *rv1 = NULL, *ref = NULL, tmp, x, y, z, zprev;
  int m = U->rows, n = U->cols;

  rv1 = (tmVal*) calloc(n, sizeof(tmVal));
  if (!rv1) {
    *e = TM_ERR_NO_MEMORY;
    return;
  }

  /* householder reduction to bidiagonal form */
  for (i = 0; i < n; i++) {
    l = i+1;
    rv1[i] = scale*g;
    g = s = scale = 0.0;
    if (i < m) {
      for (k = i; k < m; k++)
        scale += fabs(*tm_at(U,k,i));
      if (scale > 0) {
        for (k = i; k < m; k++) {
          ref = tm_at(U,k,i);
          *ref /= scale;
          s += (*ref) * (*ref);
        }
        f = *tm_at(U,i,i);
        g = sqrt(s);
        if (f > 0) g = -g;
        h = f*g - s;
        *tm_at(U,i,i) = f-g;
        for (j = l; j < n; j++) {
          s = 0.0;
          for (k = i; k < m; k++) {
            s += (*tm_at(U,k,i)) * (*tm_at(U,k,j));
          }
          f = s/h;
          for (k = i; k < m; k++) {
            *tm_at(U,k,j) += f * (*tm_at(U,k,i));
          }
        }
        for (k = i; k < m; k++)
          *tm_at(U,k,i) *= scale;
      }
    }
    *tm_at(W,i,i) = scale*g;
    g = s = scale = 0.0;
    if (i < m && i != n-1) {
      for (k = l; k < n; k++) 
        scale += fabs(*tm_at(U,i,k));
      if (scale > 0) {
        for (k = l; k < n; k++) {
          ref = tm_at(U,i,k);
          *ref /= scale;
          s += (*ref) * (*ref);
        }
        ref = tm_at(U,i,l);
        f = *ref;
        g = sqrt(s);
        if (f > 0) g = -g;
        h = f*g - s;
        *ref = f-g;
        for (k = l; k < n; k++) 
          rv1[k] = *tm_at(U,i,k)/h;
        for (j = l; j < m; j++) {
          s = 0.0;
          for (k = l; k < n; k++) 
            s += (*tm_at(U,j,k)) * (*tm_at(U,i,k));
          for (k = l; k < n; k++) 
            *tm_at(U,j,k) += s*rv1[k];
        }
        for (k = l; k < n; k++) 
          *tm_at(U,i,k) *= scale;
      }
    }
    tmp = fabs(*tm_at(W,i,i)) + fabs(rv1[i]);
    anorm = (tmp > anorm) ? tmp : anorm;
  }

  /* accumulation of right-hand transformations */
  for (i = n-1; i >= 0; i--) {
    if (i < n-1) {
      if (g) {
        for (j = l; j < n; j++) {
          *tm_at(V,j,i) = (*tm_at(U,i,j)) / (*tm_at(U,i,l)) / g;
        }
        for (j = l; j < n; j++) {
          s = 0.0;
          for (k = l; k < n; k++) 
            s += (*tm_at(U,i,k)) * (*tm_at(V,k,j));
          for (k = l; k < n; k++) 
            *tm_at(V,k,j) += s * (*tm_at(V,k,i));
        }
      }
      for (j = l; j < n; j++) {
        *tm_at(V,i,j) = 0.0;
        *tm_at(V,j,i) = 0.0;
      }
    }
    *tm_at(V,i,i) = 1.0;
    g = rv1[i];
    l = i;
  }

  /* accumulation of left-hand transformations */
  i = (m < n) ? m : n;
  i -= 1;
  for (; i >= 0; i--) {
    l = i+1;
    g = *tm_at(W,i,i);
    for (j = l; j < n; j++) 
      *tm_at(U,i,j) = 0.0;
    if (g) {
      g = 1.0/g;
      for (j = l; j < n; j++) {
        s = 0.0;
        for (k = l; k < m; k++) 
          s += (*tm_at(U,k,i)) * (*tm_at(U,k,j));
        f = (s / *tm_at(U,i,i)) * g;
        for (k = i; k < m; k++) 
          *tm_at(U,k,j) += f * (*tm_at(U,k,i));
      }
      for (j = i; j < m; j++) 
        *tm_at(U,j,i) *= g;
    } else {
      for (j = i; j < m; j++) 
        *tm_at(U,j,i) = 0.0;
    }
    *tm_at(U,i,i) += 1;
  }

  /* diagonalization of the bidiagonal form */
  for (k = n-1; k >= 1; k--) {
    for (its = 1; its <= 30; its++) {
      zprev = *tm_at(W,k,k);
      flag = 1;
      for (l = k; l >= 1; l--) {
        nm = l-1;
        if ((fabs(rv1[l]) + anorm) == anorm) {
          flag = 0;
          break;
        }
        if ((fabs(*tm_at(W,nm,nm)) + anorm) == anorm) 
          break;
      }
      if (flag) {
        c = 0.0; s = 1.0;
        for (i = l; i <= k; i++) {
          f = s*rv1[i];
          rv1[i] *= c;
          if ((fabs(f) + anorm) == anorm) 
            break;
          g = *tm_at(W,i,i);
          h = pythag(f, g);
          *tm_at(W,i,i) = h;
          h = 1.0/h;
          c = g*h; s = -f*h;
          for (j = 0; j < m; j++) {
            y = *tm_at(U,j,nm);
            z = *tm_at(U,j,i);
            *tm_at(U,j,nm) = y*c + z*s;
            *tm_at(U,j,i)  =-y*s + z*c;
          }
        }
      }
      z = *tm_at(W,k,k);
      //if (its > 10 && fabs(z-zprev) < 1E-16) {
      if (l == k) {
        if (z < 0.0) {
          *tm_at(W,k,k) = -z;
          for (j = 0; j < n; j++) 
            *tm_at(V,j,k) *= -1;
        }
        break;
      }
      if (its == 30) {
        *e = TM_ERR_NO_SOLUTN;
        /* error */
        return;
      }  
      x = *tm_at(W,l,l);
      nm = k-1;
      y = *tm_at(W,nm,nm);
      g = rv1[nm];
      h = rv1[k];
      f = ((y-z)*(y+z) + (g-h)*(g+h)) / (2*h*y);
      g = pythag(f, 1.0);
      tmp = (f >= 0) ? g : (-g);
      f = ((x-z)*(x+z) + h*(y/(f+tmp)-h)) / x;
      /* next QR transform */
      c = s = 1.0;
      for (j = l; j <= nm; j++) {
        i = j+1;
        g = rv1[i];
        y = *tm_at(W,i,i);
        h = s*g;
        g = c*g;
        z = pythag(f, h);
        rv1[j] = z;
        c = f/z; s = h/z;
        f = x*c + g*s;
        g =-x*s + g*c;
        h = y*s;
        y *= c;
        for (jj = 0; jj < n; jj++) {
          x = *tm_at(V,jj,j);
          z = *tm_at(V,jj,i);
          *tm_at(V,jj,j) = x*c + z*s;
          *tm_at(V,jj,i) =-x*s + z*c;
        }
        z = pythag(f, h);
        *tm_at(W,j,j) = z;
        if (z > 0) {
          z = 1.0/z;
          c = f*z; s = h*z;
        }
        f = c*g + s*y;
        x =-s*g + c*y;
        for (jj = 0; jj < m; jj++) {
          y = *tm_at(U,jj,j);
          z = *tm_at(U,jj,i);
          *tm_at(U,jj,j) = y*c + z*s;
          *tm_at(U,jj,i) =-y*s + z*c;
        }
      }

      rv1[l] = 0.0;
      rv1[k] = f;
      *tm_at(W,k,k) = x;
    }
  }

  free(rv1);
}

void tf_svd(tMat* U, tMat* S, tMat* V, tMat* m, int* err)
{
  int e = 0, i, j;
  int nr = m->rows, nc = m->cols;
  tmVal *row_i;

  TM_ASSERT_ARGS(m && U && S && V && IS_UNIQUE4(m,U,S,V), e, end_svd);

  if (tm_relevant(U,nr,nc,&e) && tm_relevant(S,nc,nc,&e) && tm_relevant(V,nc,nc,&e)) {
    if (!IS_PRIM(U) || !IS_PRIM(S) || !IS_PRIM(V)) {
      e = TM_ERR_NOT_MAIN;
      goto end_svd;
    }
    /* copy source matrix */
    for (i = 0; i < nr; i++) {
      row_i = U->data + i*nc;
      for (j = 0; j < nc; j++)
        row_i[j] = *tm_at(m,i,j);
    }
    tm_zeros(S);
    //tm_zeros(V);
    /* find SVD */
    svdcmp(U, V, S, &e);
      
  }


end_svd:
  if (err) *err = e;
}
