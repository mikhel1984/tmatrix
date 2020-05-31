/**
 * @file tmatrix_calc.c
 * @author Stanislav Mikhel
 * @date 2019
 * @brief Arithmetic operations and transformations.
 */ 
#include <stdlib.h>
#include <math.h>
#include "tmatrix.h"
#include "tmatrix_priv.h"

#define TINY         1E-20
#define PINV_TOL_MAX 1E100
#define PINV_TOL     1E-9
#define MEM_TMP_VEC 8


int tm_add(tMat *dst, tMat* m, int* err) 
{
  int i,j,R,C, e = 0;
  tmVal *p1, *p2;

  TM_ASSERT_ARGS(dst && m, e, end_add);
   
  R = m->rows; C = m->cols;         
  if(dst->rows == R && dst->cols == C) {
    /* equal */
    if(IS_PRIM(dst) && IS_PRIM(m)) {
      /* just add element-wise */
      R *= C;    /* reuse variable */
      p1 = dst->data; p2 = m->data;
      for(i = 0; i < R; i++)  
        *p1++ += *p2++;                       
    } else {
      /* use method 'at' */        
      for(i = 0; i < R; i++) {
        for(j = 0; j < C; j++) 
          *tm_at(dst,i,j) += *tm_at(m,i,j);	          
      }
    }           
  } else 
    e = TM_ERR_WRONG_SIZE;                  
           
end_add:
  if(err) *err = e;
  return !e;
}

int tm_sub(tMat *dst, tMat* m, int* err) 
{
  int i,j,R,C, e = 0;
  tmVal *p1,*p2;

  TM_ASSERT_ARGS(dst && m, e, end_sub);
   
  R = m->rows; C = m->cols;
  if(dst->rows == R && dst->cols == C) {
    /* equal */
    p1 = dst->data; p2 = m->data;
    if(IS_PRIM(dst) && IS_PRIM(m)) {
      /* element-wise */
      R *= C;
      p1 = dst->data; p2 = m->data;
      for(i = 0; i < R; i++) 
        *p1++ -= *p2++;              
    } else {
      /* use method 'at' */
      for(i = 0; i < R; i++) {
        for(j = 0; j < C; j++) 
          *tm_at(dst,i,j) -= *tm_at(m,i,j);	          
      }
    }           
  } else 
    e = TM_ERR_WRONG_SIZE;

end_sub:
  if(err) *err = e;
  return !e;
}

int tm_scale(tMat *dst, tmVal k, int* err)
{
  int R,C,i,j, e = 0;
  tmVal *p;

  TM_ASSERT_ARGS(dst, e, end_scale);
   
  R = dst->rows; C = dst->cols;
  if(IS_PRIM(dst)) {
    R *= C;
    p = dst->data;
    for(i = 0; i < R; i++)
      *p++ *= k;    
  } else {
    for(i = 0; i < R; i++) {
      for(j = 0; j < C; j++) {
        *tm_at(dst,i,j) *= k;
      }
    }
  }
   
end_scale:
  if(err) *err = e;
  return !e;
}


int tm_mul(tMat* dst, tMat *a, tMat *b, int *err)
{
  tmSize R1,C1,C2;
  int i,j,k, e = 0;
  tmVal acc;

  TM_ASSERT_ARGS(dst && a && b, e, end_mul);
   
  if(dst != a && dst != b) {
    /* check matrices for product */
    R1 = a->rows; C1 = a->cols;
    C2 = b->cols;
    if(C1 == b->rows) {
      /* check destination */
      if(tm_relevant(dst,R1,C2,err)) {  
        /* evaluate result */
        for(i = 0; i < R1; i++) {
          for(j = 0; j < C2; j++) {
            acc = 0;
            for(k = 0; k < C1; k++) 
              acc += (*tm_at(a,i,k)) * (*tm_at(b,k,j));            
            *tm_at(dst,i,j) = acc;
          }
        } 
      }     
    } else 
      e = TM_ERR_NOT_COMPAT;         
  } else 
    e = TM_ERR_NOT_DEF;
   
end_mul:
  if(err) *err = e;  
  return !e;
}

void ludcmp(tMat *dst, int indx[], tmVal* d, int *err)
{
  int i, imax, j, k, n;
  tmVal big, dum, sum, *irow, *vv = NULL;
  tmVal arr[MEM_TMP_VEC] = {0};
   
  n = dst->rows;
  /* allocate memory */
  vv = (n <= MEM_TMP_VEC) ? arr : (tmVal*) malloc(n * sizeof(tmVal));
  if(!vv) {
    if(err) *err = TM_ERR_NO_MEMORY;
    goto end_ludcmp;
  }  
   
  *d = 1.0;
  /* get normalization coefficient */
  for(i = 0; i < n; i++) {
    big = 0.0;
    irow = dst->data + i * n;
    for(j = 0; j < n; j++) {
      dum = fabs(irow[j]);
      if(dum > big) big = dum;
    }
    if(big == 0.0) {
      if(err) *err = TM_ERR_NO_SOLUTN;
      goto end_ludcmp;
    }
    vv[i] = 1.0/big;     
  }
  /* evaluate solution */
  for(j = 0; j < n; j++) {
    for(i = 0; i < j; i++) {
      irow = dst->data + i * n;
      sum = irow[j];
      for(k = 0; k < i; k++) 
        sum -= irow[k] * dst->data[k*n + j];
      irow[j] = sum;
    }
    big = 0.0;
    for(i = j; i < n; i++) {
      irow = dst->data + i * n;
      sum = irow[j];
      for(k = 0; k < j; k++) 
        sum -= irow[k] * dst->data[k*n + j];
      irow[j] = sum;
      dum = vv[i]*fabs(sum);
      if(dum >= big) {
        big = dum;
        imax = i;
      }
    }
    if(j != imax) {
      irow = dst->data + imax*n;
      for(k = 0; k < n; k++) {
        dum = irow[k];
        irow[k] = dst->data[j*n + k];
        dst->data[j*n + k] = dum;
      }
      *d = -(*d);
      vv[imax] = vv[j];
    }
    if(indx) indx[j] = imax;
    if(dst->data[j*n + j] == 0.0) dst->data[j*n + j] = TINY;  /* avoid singularity */
    if(j != n-1) {
      dum = 1.0 / dst->data[j*n + j];
      for(i = j+1; i < n; i++) 
        dst->data[i*n+j] *= dum;         
    }
  }  

end_ludcmp:  
  if(n > MEM_TMP_VEC) free(vv);  
}

void lubksb(tMat *src, int indx[], tMat *b)
{
  int i, ii = -1, ip, j, n;
  tmVal sum, *irow;
   
  n = src->rows;
  for(i = 0; i < n; i++) {
    ip = indx[i];
    sum = *tm_at(b,ip,0);
    *tm_at(b,ip,0) = *tm_at(b,i,0);      
    if(ii > -1) {
      irow = src->data + i*n;
      for(j = ii; j < i; j++) 
        sum -= irow[j] * (*tm_at(b,j,0)); 
    } else if(sum) {
      ii = i;         
    }
    *tm_at(b,i,0) = sum;      
  }
  for(i = n-1; i >= 0; i--) {
    sum = *tm_at(b,i,0);
    irow = src->data + i*n;
    for(j = i+1; j < n; j++) 
      sum -= irow[j] * (*tm_at(b,j,0));
    *tm_at(b,i,0) = sum / irow[i];
  }   
}

tmVal tm_det(tMat *m, int *err)
{
  int e = 0, n, j, n1, n2;
  tmVal d;   
  tMat tmp = NULL_TMATRIX; 

  TM_ASSERT_ARGS(m, e, end_det);
   
  if(m->rows == m->cols) {
    n = m->rows;
    tmp = tm_copy(m,&e);                 
    if(!e) ludcmp(&tmp,NULL,&d,&e);
    if(!e) {
      n1 = n+1; n2 = n*n;
      /* product of diagonal elements */
      for(j = 0; j < n2; j += n1)
        d *= tmp.data[j];
    }                
  } else 
    e = TM_ERR_NOT_DEF;
      
end_det:
  if(err) *err = e;  
  tm_clear(&tmp);
   
  return d;
}

tMat tm_inv(tMat *m, int *err)
{
  int e = 0, n = 0, *idx = NULL, j;
  int arr[MEM_TMP_VEC] = {0};
  tmVal d;
  tMat tmp = NULL_TMATRIX, res = NULL_TMATRIX, 
       col;
   
  TM_ASSERT_ARGS(m, e, end_inv);

    if(m->rows == m->cols) {
      n = m->rows;
      tmp = tm_copy(m,&e);
      if(!e) res = tm_eye(n,n,&e);         
      if(!e) {        
        idx = (n <= MEM_TMP_VEC) ? arr : (int*) malloc(n * sizeof(int));        
        if(idx) {
          ludcmp(&tmp,idx,&d,&e);
          if(!e) {
            for(j = 0; j < n; j++) {
              col = tm_block(&res,0,j,n,1,&e);
              if(e) break;                     
              lubksb(&tmp,idx,&col);
            }
          }
        } else 
          e = TM_ERR_NO_MEMORY;            
      }
    } else 
      e = TM_ERR_NOT_DEF;
  
end_inv:
  if(n > MEM_TMP_VEC) free(idx);
  tm_clear(&tmp);
   
  if(err) *err = e;
  return res;
}

tMat pinva(tMat *src, int* transp, tmVal* tolerance, int* err)
{
  int e = 0, m, n, i;
  tmVal tol = PINV_TOL_MAX, v;
  tMat A = NULL_TMATRIX, st;
   
  TM_ASSERT_ARGS(src, e, end_pinva);

  *transp = 0;
  m = src->rows; n = src->cols;
  st = tm_T(src,&e); 		if(e) goto end_pinva;
  if(m < n) {
    /* transpose */
    *transp = 1;
    n = m;
    A = tm_new(n,n,&e); 	if(e) goto end_pinva;
    tm_mul(&A,src,&st,&e); 	if(e) goto end_pinva;        
  } else {
    /* original */
    A = tm_new(n,n,&e); 	if(e) goto end_pinva;
    tm_mul(&A,&st,src,&e); 	if(e) goto end_pinva;
  }
  /* tolerance */
  m = n+1; n *= n; 
  for(i = 0; i < n; i += m) {
    v = A.data[i];
    if(v <= 0) 
      v = PINV_TOL_MAX;
    tol = (tol < v) ? tol : v;
  }
  tol *= PINV_TOL;    

end_pinva:    
  *err = e;
  *tolerance = tol;
  return A;
}

tMat pinvl(tMat *A, tmVal tol, int *rr, int *err)
{
  int e = 0, n, r = 0, k, xx, n2;
  tmVal lkr;
  tMat L = NULL_TMATRIX, B = NULL_TMATRIX,
       tmp1 = NULL_TMATRIX, tmp2 = NULL_TMATRIX,
       blk;
  
  n = A->rows;
  n2 = n * n;
  /* prepare matrices */
  L = tm_new(n,n,&e); 				if(e) goto end_pinvl;
  B = tm_new(n,1,&e); 				if(e) goto end_pinvl;
  tmp1 = tm_new(n,1,&e); 			if(e) goto end_pinvl;
  tmp2 = tm_new(n,1,&e); 			if(e) goto end_pinvl;
  
  for(k = 0; k < n; k++) {
    r++;
    /* get column */
    blk = tm_block(A, k,k,n-k,1, &e); 		if(e) goto end_pinvl;
    B.rows = n-k; B.cols = 1;
    tm_insert(&B, &blk, &e); 			if(e) goto end_pinvl;
    
    if(r > 1) {
      /* correct B matrix */
      tmp1.rows = 1; tmp1.cols = r-1;
      blk = tm_block(&L, k,0,1,r-1, &e); 	if(e) goto end_pinvl;
      tm_insert(&tmp1, &blk, &e); 		if(e) goto end_pinvl;
      tmp1.rows = r-1; tmp1.cols = 1;    /* "transpose */
      blk = tm_block(&L, k,0,n-k,r-1, &e); 	if(e) goto end_pinvl;
      if(!(tm_mul(&tmp2,&blk,&tmp1,&e) && tm_sub(&B,&tmp2,&e))) goto end_pinvl;       
    }
    
    blk = tm_block(&L, k,r-1,n-k,1, &e); 	if(e) goto end_pinvl;
    tm_insert(&blk,&B,&e); 			if(e) goto end_pinvl;
    
    lkr = L.data[k*n+(r-1)];
    if(lkr > tol) {
      /* normalize column */
      lkr = sqrt(lkr);
      L.data[k*n+(r-1)] = lkr;
      if(k < n) {
        for(xx = (k+1)*n+(r-1); xx < n2; xx += n) {
          L.data[xx] /= lkr;
        } 
      }
    } else 
      r--; 
  } 
  
end_pinvl:
  /* free memory */
  tm_clear(&B);
  tm_clear(&tmp1);
  tm_clear(&tmp2);
  
  *err = e;
  *rr = r;
  
  return L;  
}

tMat tm_pinv(tMat *src, int *err)
{
  int e = 0, transp = 0, r, n;
  tmVal tol = 0;
  tMat A = NULL_TMATRIX, L = NULL_TMATRIX, prod1 = NULL_TMATRIX,
       LL = NULL_TMATRIX, M = NULL_TMATRIX, prod2 = NULL_TMATRIX,
       blk, Lt;
  /* main matrices */
  A = pinva(src,&transp,&tol,&e); 	if(e) goto end_pinv;
  n = A.rows;
  L = pinvl(&A,tol,&r,&e); 		if(e) goto end_pinv;
  
  /* auxiliary matrices */
  blk = tm_block(&L, 0,0,n,r, &e); 	if(e) goto end_pinv;
  LL = tm_copy(&blk,&e); 		if(e) goto end_pinv;
  Lt = tm_T(&LL,&e); 			if(e) goto end_pinv;
    
  prod1 = tm_new(r,r,&e); 		if(e) goto end_pinv; 
  tm_mul(&prod1, &Lt, &LL, &e); 	if(e) goto end_pinv;
  M = tm_inv(&prod1, &e); 		if(e) goto end_pinv;
  blk = tm_T(src,&e); 			if(e) goto end_pinv;
    
  /* find pseudo inverse*/
  prod2 = tm_simp();
  if(transp) {
    /* blk*LL*M*M*Lt  */
    if(!(tm_mul(&prod2,&blk,&LL,&e)    && tm_mul(&prod1,&prod2,&M,&e) &&
         tm_mul(&prod2,&prod1,&M,&e)  && tm_mul(&prod1,&prod2,&Lt,&e))) { goto end_pinv; }
  } else {
    /* LL*M*M*Lt*blk */
    if(!(tm_mul(&prod2,&LL,&M,&e)     && tm_mul(&prod1,&prod2,&M,&e) &&
         tm_mul(&prod2,&prod1,&Lt,&e) && tm_mul(&prod1,&prod2,&blk,&e))) { goto end_pinv; }
  } 
  
end_pinv:
  tm_clear(&A);
  tm_clear(&L);
  tm_clear(&LL);
  tm_clear(&M);
  tm_clear(&prod2);
  
  if(err) *err = e; 
  
  return prod1;
}

int tm_rank(tMat* m, int* err)
{
  int e = 0, res = 0, i,j, R = 0,C, q;
  tMat cp = NULL_TMATRIX, tmp; 
  tmVal k, *swp, **ptr = NULL, *arr[MEM_TMP_VEC];

  TM_ASSERT_ARGS(m, e, end_rank);

    /* make copy for modification */
  if(m->rows > m->cols) {
    tmp = tm_T(m,&e); 
    if(!e) cp = tm_copy(&tmp,&e);
  } else 
    cp = tm_copy(m,&e);
  if(!e) {
    /* use pointers to rows instead of copying elements */
    R = cp.rows; C = cp.cols;
    ptr = (R <= MEM_TMP_VEC) ? arr : (tmVal**) malloc(sizeof(tmVal*) * R);
    if(!ptr) {
      if(err) *err = TM_ERR_NO_MEMORY;
      return 0;
    }
    /* initialize rows pointers */
    ptr[0] = cp.data;
    for(i = 1; i < R; i++) ptr[i] = ptr[i-1] + C;
    /* triangulate */
    for(i = 0; i < R; i++) {
      /* looking for nonzero elements */
      for(j = i+1; j < R && ptr[i][i] == 0; j++) {
        if(ptr[j][i] != 0) {
          swp = ptr[j]; ptr[j] = ptr[i]; ptr[i] = swp;
        }
      }
      k = ptr[i][i];
      if(k != 0) {
        /* normalize to k */
        for(j = i; j < C; j++) ptr[i][j] /= k;
          /* subtract */
          for(q = i+1; q < R; q++) {
            k = ptr[q][i];
            for(j = i; j < C; j++) ptr[q][j] -= k * ptr[i][j];
          }
          res++;  /* increase rank counter */
      } else 
        break;
    }
  }

end_rank:
  if(err) *err = e;
  tm_clear(&cp);
  if(R > MEM_TMP_VEC) free(ptr);

  return res;
}
