
#include <stdlib.h>
#include <math.h>
#include "tmatrix.h"
#include "tmatrix_priv.h"

#define TINY         1E-20
#define PINV_TOL_MAX 1E100
#define PINV_TOL     1E-9
#define PINV_TMP_VEC 8

int tm_add(tMat *dst, tMat* m, int* err) 
{
  int i,j,R,C;
  int e = 0;
  tmVal *p1, *p2;
   
  if(dst && m) {         
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
  } else    
    e = TM_ERR_EMPTY_ARGS;
           
  if(err) *err = e;
   
  return !e;
}

int tm_sub(tMat *dst, tMat* m, int* err) 
{
  int i,j,R,C;
  int e = 0;
  tmVal *p1,*p2;
   
  if(dst && m) {
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
  } else    
    e = TM_ERR_EMPTY_ARGS;
     
  if(err) *err = e;
   
  return !e;
}

int tm_scale(tMat *dst, tmVal k, int* err)
{
  int R,C,i,j;
  int e = 0;
  tmVal *p;
   
  if(dst) {
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
  } else 
    e = TM_ERR_EMPTY_ARGS;
   
  if(err) *err = e;
   
  return !e;
}


int tm_mul(tMat* dst, tMat *a, tMat *b, int *err)
{
  tmSize R1,C1,R2,C2,R3,C3;
  int N,i,j,k;
  tmVal* data, acc;
  int e = 0;
   
  if(dst && a && b) {
    if(dst != a && dst != b) {
      /* check matrices for product */
      R1 = a->rows; C1 = a->cols;
      R2 = b->rows; C2 = b->cols;
      if(C1 == R2) {
        /* check destination */
        R3 = dst->rows; C3 = dst->cols; 
        if(R3 != R1 || C3 != C2) { 
          /* if the size is equal, any type is acceptible, otherwise the size will be changed */
          if(IS_PRIM(dst)) {
            N = R1 * C2;
            if(R3 * C3 < N) {
              /* if the memory is static, it cannot be expanded... */
              if(dst->type == TM_MAIN) {
                /* get new memory */
                data = (tmVal*) malloc(N*sizeof(tmVal));
                if(data) {
                  free(dst->data);
                  dst->data = data;  
                } else 
                  e = TM_ERR_NO_MEMORY;                
              } else 
                e = TM_ERR_NOT_MAIN;                     
            }               
          } else 
            e = TM_ERR_NOT_MAIN;
          if(!e) {
            /* change size */
            dst->rows = R1;
            dst->cols = C2;
            dst->width = C2;
          } else {
            if(err) *err = e;
            return !e;
          }
        }   
        /* evaluate result */
        for(i = 0; i < R1; i++) {
          for(j = 0; j < C2; j++) {
            acc = 0;
            for(k = 0; k < C1; k++) 
              acc += (*tm_at(a,i,k)) * (*tm_at(b,k,j));            
            *tm_at(dst,i,j) = acc;
          }
        }      
      } else 
        e = TM_ERR_NOT_COMPAT;         
    } else 
      e = TM_ERR_NOT_DEF;
  } else 
    e = TM_ERR_EMPTY_ARGS;
   
  if(err) *err = e;  
   
  return !e;
}

void ludcmp(tMat *dst, int indx[], tmVal* d, int *err)
{
  int i, imax, j, k, n;
  tmVal big, dum, sum, tmp, *irow, *vv = NULL;
  tmVal arr[PINV_TMP_VEC] = {0};
   
  n = dst->rows;
  /* allocate memory */
  vv = (n <= PINV_TMP_VEC) ? arr : (tmVal*) malloc(n * sizeof(tmVal));
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
      tmp = fabs(irow[j]);
      if(tmp > big) big = tmp;
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
  if(n > PINV_TMP_VEC) free(vv);  
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
   
  if(m) {
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
  } else 
    e = TM_ERR_EMPTY_ARGS;
      
  if(err) *err = e;  
  tm_clear(&tmp);
   
  return d;
}

tMat tm_inv(tMat *m, int *err)
{
  int e = 0, n = 0, *idx = NULL, j;
  int arr[PINV_TMP_VEC] = {0};
  tmVal d;
  tMat tmp = NULL_TMATRIX, res = NULL_TMATRIX, 
       col;
   
  if(m) {
    if(m->rows == m->cols) {
      n = m->rows;
      tmp = tm_copy(m,&e);
      if(!e) res = tm_eye(n,n,&e);         
      if(!e) {        
        idx = (n <= PINV_TMP_VEC) ? arr : (int*) malloc(n * sizeof(int));        
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
  } else 
    e = TM_ERR_EMPTY_ARGS;
  
  if(n > PINV_TMP_VEC) free(idx);
  tm_clear(&tmp);
   
  if(err) *err = e;
     
  return res;
}

tMat pinva(tMat *src, int* transp, tmVal* tolerance, int* err)
{
  int e = 0, m, n, n1, n2, i;
  tmVal tol = PINV_TOL_MAX, v;
  tMat A = NULL_TMATRIX, st;
   
  if(src) {
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
    n2 = n*n; n1 = n+1;
    for(i = 0; i < n2; i += n1) {
      v = A.data[i];
      if(v <= 0) 
        v = PINV_TOL_MAX;
      tol = (tol < v) ? tol : v;
    }
    tol *= PINV_TOL;    
  } else
    e = TM_ERR_EMPTY_ARGS;

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
    /* get culumn */
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
       blk, Lt, st;
  /* main matrices */
  st = tm_T(src,&e); 			if(e) goto end_pinv;
  A = pinva(src,&transp,&tol,&e); 	if(e) goto end_pinv;
  n = A.rows;
  L = pinvl(&A,tol,&r,&e); 		if(e) goto end_pinv;
  
  /* auxilary matrices */
  blk = tm_block(&L, 0,0,n,r, &e); 	if(e) goto end_pinv;
  LL = tm_copy(&blk,&e); 		if(e) goto end_pinv;
  Lt = tm_T(&LL,&e); 			if(e) goto end_pinv;
    
  prod1 = tm_new(r,r,&e); 		if(e) goto end_pinv; 
  tm_mul(&prod1, &Lt, &LL, &e); 	if(e) goto end_pinv;
  M = tm_inv(&prod1, &e); 		if(e) goto end_pinv;
    
  /* find pseudoinverse */
  prod2 = tm_new(1,1,&e); 		if(e) goto end_pinv;
  if(transp) {
    /* st*LL*M*M*Lt  */
    if(!(tm_mul(&prod2,&st,&LL,&e)    && tm_mul(&prod1,&prod2,&M,&e) &&
         tm_mul(&prod2,&prod1,&M,&e)  && tm_mul(&prod1,&prod2,&Lt,&e))) { goto end_pinv; }
  } else {
    /* LL*M*M*Lt*st */
    if(!(tm_mul(&prod2,&LL,&M,&e)     && tm_mul(&prod1,&prod2,&M,&e) &&
         tm_mul(&prod2,&prod1,&Lt,&e) && tm_mul(&prod1,&prod2,&st,&e))) { goto end_pinv; }
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
