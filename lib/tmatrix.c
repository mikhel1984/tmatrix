/**
 * @file tmatrix.c
 * @author Stanislav Mikhel
 * @date 2020
 * @brief Main structures and functions for manipulation with matrices.
 */ 
#include <stdlib.h>
#include "tmatrix.h"
#include "tmatrix_priv.h"

/* "Private" method for access to matrix elements without additional checking. */
tmVal* tm_at(tMat* m, tmSize r, tmSize c) {
  return (m->type == TM_TRANSPOSE) ? (m->data + (c * m->width + r)) 
                                   : (m->data + (r * m->width + c));
}

/* "Private", check size and allocate memory if need */
int tm_relevant(tMat* m, tmSize R, tmSize C, int* err)
{
  int e = 0, N;
  tmVal *data;
  
  if(m->rows != R || m->cols != C) {
    if(m->type == TM_ALLOC) {
      N = R * C;
      if(m->rows * m->cols < N) {
        data = (tmVal*) malloc(N * sizeof(tmVal));
        if(data) {
          free(m->data);
          m->data = data;
        } else 
          e = TM_ERR_NO_MEMORY;
      }
      if(!e) {
        m->rows = R;
        m->width = m->cols = C;
      }
    } else 
      e = TM_ERR_NOT_MAIN;
  }
  
  if(err) *err = e;
  return !e;
}

/* New zero matrix */ 
tMat tm_new(tmSize r, tmSize c, int* err) 
{
  int e = 0;
  tMat res = NULL_TMATRIX;   
  res.type = TM_ALLOC;
   
  if(r && c) {
    res.data = (tmVal*) calloc(r*c, sizeof(tmVal));
    if(res.data) {
      res.rows = r;
      res.width = res.cols = c;    
    } else 
      e = TM_ERR_NO_MEMORY;      
  } else if(!r && !c) {    
    /* allow empty matrix for further usage in tm_mul */
  } else 
    e = TM_ERR_WRONG_SIZE;

  if(err) *err = e;
  return res;
}

/* New identity matrix */
tMat tm_eye(tmSize r, tmSize c, int *err)
{
  int e = 0, n1 = c+1, n2 = r*c, i;
  tMat res = tm_new(r,c,&e);
   
  if(!e) {
    for(i = 0; i < n2; i += n1)
      res.data[i] = 1;
  }
   
  if(err) *err = e;
  return res;
}

/* Matrix with static preallocated memory */
tMat tm_static(tmSize r, tmSize c, tmVal dat[], int* err)
{
  tMat res = NULL_TMATRIX;
  int e = 0;
  
  TM_ASSERT_ARGS(dat, e, end_static);
  TM_ASSERT_INDEX(r && c, e, end_static);
   
  res.data = dat;
  res.rows = r;
  res.width = res.cols = c;
  res.type = TM_STATIC;      
   
end_static:
  if(err) *err = e;
  return res;
}

/* Free dynamically allocated memory */
void tm_clear(tMat* matrix) 
{
  if(matrix && matrix->type == TM_ALLOC) {
    free(matrix->data);
    matrix->data = NULL;
    matrix->rows = 0;
    matrix->cols = 0;
  }
}

/* "Protected" getter */
tmVal tm_get(tMat* m, tmSize r, tmSize c, int* err)
{
  int e = 0;
  tmVal res = 0;
  
  TM_ASSERT_ARGS(m, e, end_get);
  TM_ASSERT_INDEX(r < m->rows && c < m->cols, e, end_get);

  res = *tm_at(m,r,c);
      
end_get:
  if(err) *err = e;
  return res;
}

/* "Protected" setter */
void tm_set(tMat* m, tmSize r, tmSize c, tmVal v, int* err)
{
  int e = 0;

  TM_ASSERT_ARGS(m, e, end_set);
  TM_ASSERT_INDEX(r < m->rows && c < m->cols, e, end_set);
   
  *tm_at(m,r,c) = v;
   
end_set:
  if(err) *err = e;  
}

/* Get row number */
tmSize tm_rows(tMat* m) 
{ 
  return m ? m->rows : 0;
}

/* Get column number */
tmSize tm_cols(tMat* m) 
{
  return m ? m->cols : 0;
}

/* Create matrix copy */
tMat tm_copy(tMat* src, int* err)
{
  tmVal* data = NULL, *p;
  int e = 0, i, j;
  tMat res = NULL_TMATRIX;

  TM_ASSERT_ARGS(src, e, end_copy);
   
  res.rows = src->rows;
  res.cols = src->cols;
  res.width = res.cols;
  j = res.rows*res.cols;
  TM_ASSERT_INDEX(j, e, end_copy);

  data = (tmVal*) malloc(j * sizeof(tmVal));
  if(data) {
    res.data = data;
    res.type = TM_ALLOC;        
    if(IS_PRIM(src)) {
      p = src->data;
      for(i = 0; i < j; i++) 
        *data++ = *p++;               
    } else {
      for(i = 0; i < res.rows; i++) {
        for(j = 0; j < res.cols; j++) 
          *data++ = *tm_at(src,i,j);	          
      }
    }
  } else 
    e = TM_ERR_NO_MEMORY;         
      
end_copy:
  if(err) *err = e;
  return res;
}


/* Get transposed version of the matrix */
tMat tm_T(tMat* src, int *err) 
{
  tMat res = NULL_TMATRIX;
  int e = 0;    

  TM_ASSERT_ARGS(src, e, end_T);
  
  if(IS_PRIM(src)) {
    res = *src;
    res.cols = src->rows;
    res.rows = src->cols;
    res.type = TM_TRANSPOSE;
  } else 
    e = TM_ERR_NOT_MAIN;
   
end_T:
  if(err) *err = e;
  return res;
}

/* Get sub-matrix */
tMat tm_block(tMat* src, tmSize r0, tmSize c0, tmSize Nr, tmSize Nc, int *err)
{
  tMat res = NULL_TMATRIX;
  int e = 0;
   
  TM_ASSERT_ARGS(src, e, end_block);
  TM_ASSERT_INDEX(Nr > 0 && (r0+Nr) <= src->rows && Nc > 0 && (c0+Nc) <= src->cols, e, end_block);

  if(IS_PRIM(src)) {
    res = *src;
    res.type = TM_SUB;
    res.rows = Nr;
    res.cols = Nc;
    res.data = src->data + (r0*src->width + c0);
  } else 
    e = TM_ERR_NOT_MAIN;

end_block:
  if(err) *err = e;  
  return res;
}

/* Copy values from one matrix to another */
int tm_insert(tMat *dst, tMat *src, int* err) 
{
  int e = 0, i,j,R,C;

  TM_ASSERT_ARGS(dst && src, e, end_insert);
   
  R = dst->rows;
  C = dst->cols;
  if(R == src->rows && C == src->cols) {
    for(i = 0; i < R; i++) {
      for(j = 0; j < C; j++) 
        *tm_at(dst,i,j) = *tm_at(src,i,j);            
    }
  } else 
    e = TM_ERR_DIFF_SIZE;
   
end_insert:
  if(err) *err = e;  
  return !e; 
}

/* Set all to zero */
void tm_zeros(tMat* dst)
{
  int i, N;
  tmVal* ptr;
  if(dst && dst->data) {
    N = dst->rows * dst->cols;
    ptr = dst->data;
    for(i = 0; i < N; i++) *ptr++ = 0;
  }
}

/* Create matrix from list of matrices and some rule */
tMat tm_make(tMat src[], tmSize N, tmSize R, tmSize C, 
             tmVal (*rule)(tMat*,tmSize,tmSize,tmSize,int*), int* err)
{
  int e = 0,i,j;
  tmVal *data = NULL;
  tMat res = NULL_TMATRIX;  
  res.type = TM_ALLOC;

  TM_ASSERT_ARGS(src && rule, e, end_make);
  TM_ASSERT_INDEX(R && C, e, end_make);
  
  data = (tmVal*) calloc(R*C, sizeof(tmVal));
  if(data) {
    /* initialize matrix */
    res.data = data;
    res.rows = R;
    res.width = res.cols = C; 
    /* find values */ 
    for(i = 0; !e && i < R; i++) {
      for(j = 0; !e && j < C; j++) 
        *data++ = rule(src,N,i,j,&e);          
    }   
  } else 
    e = TM_ERR_NO_MEMORY; 
    
end_make:
  if(err) *err = e;
  return res;
}

/*		Concatenation function		*/
/* Vectical concatenation */
tmVal concatv(tMat src[], tmSize N, tmSize r, tmSize c, int* err)
{
  int i = 0;
  
  while(i < N && r >= src[i].rows) r -= src[i++].rows;
  
  if(i >= N) {
    if(err) *err = TM_ERR_WRONG_SIZE;
    return 0;
  }
  
  return *tm_at(src+i,r,c);  
}

/* Horizontal concatenation */
tmVal concath(tMat src[], tmSize N, tmSize r, tmSize c, int* err)
{
  int i = 0;
  
  while(i < N && c >= src[i].cols) c -= src[i++].cols;
  
  if(i >= N) {
    if(err) *err = TM_ERR_WRONG_SIZE;
    return 0;
  }
  
  return *tm_at(src+i,r,c);  
}

/* General case */
tMat tm_concat(tMat src[], int N, int dir, int* err)
{
  int e = 0, i, r,c;
  tMat res = NULL_TMATRIX;  

  TM_ASSERT_ARGS(src, e, end_concat);
  TM_ASSERT_INDEX(N > 0, e, end_concat);
  
  switch(dir) {
  case 0: /* horizontal concatenation */
    c = src[0].cols;
    r = src[0].rows;
    for(i = 1; i < N; i++) {
      if(src[i].rows != r) {
        e = TM_ERR_NOT_COMPAT;
        goto end_concat;
      }
      c += src[i].cols;
    }
    res = tm_make(src,N,r,c,concath,&e);
    break;
  case 1: /* vertical concatenation */
    r = src[0].rows;
    c = src[0].cols;
    for(i = 1; i < N; i++) {
      if(src[i].cols != c) {
        e = TM_ERR_NOT_COMPAT;
        goto end_concat;
      }
      r += src[i].rows;
    }
    res = tm_make(src,N,r,c,concatv,&e);
    break;
  default:
    e = TM_ERR_WRONG_SIZE;
    break;
  }

end_concat:
  if(err) *err = e;
  return res;
} 
/*		End concatenation 		*/
