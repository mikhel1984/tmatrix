/**
 * @file tmatrix.c
 * @author Stanislav Mikhel
 * @date 2019
 * @brief Main structures and functions for manipulation with matrices.
 */ 
#include <stdio.h>
#include <stdlib.h>
#include "tmatrix.h"
#include "tmatrix_priv.h"

/* "Private" method for access to matrix elements without additional checking. */
tmVal* tm_at(tMat* m, tmSize r, tmSize c) {
  return (m->type == TM_TRANSPOSE) ? (m->data + (c * m->width + r)) 
                                   : (m->data + (r * m->width + c));
}

tMat tm_new(tmSize r, tmSize c, int* err) 
{
  int e = 0;
  tMat res = NULL_TMATRIX;   
  res.type = TM_MAIN;
   
  if(r && c) {
    res.data = (tmVal*) calloc(r*c, sizeof(tmVal));
    if(res.data) {
      res.rows = r;
      res.cols = c;
      res.width = c;      
    } else 
      e = TM_ERR_NO_MEMORY;      
  } else if(!r && !c) {}    /* allow empty matrix for further usage in tm_mul */
    else 
    e = TM_ERR_WRONG_SIZE;

  if(err) *err = e;

  return res;
}

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

tMat tm_static(tmSize r, tmSize c, tmVal dat[], int* err)
{
  tMat res = NULL_TMATRIX;
  int e = 0;
   
  if(r && c) {
    if(dat) {
      res.data = dat;
      res.rows = r;
      res.cols = c;
      res.width = c;
      res.type = TM_STATIC;      
    } else 
      e = TM_ERR_EMPTY_ARGS;
  } else 
    e = TM_ERR_WRONG_SIZE;
   
  if(err) *err = e;
   
  return res;
}

void tm_clear(tMat* matrix) 
{
  if(matrix && matrix->type == TM_MAIN) {
    free(matrix->data);
    matrix->data = NULL;
    matrix->rows = 0;
    matrix->cols = 0;
  }
}

int tm_init(tMat* dst, tmVal src[], int* err) {
  int i,j,R,C, e=0;
  tmVal *data;
   
  if(dst && src) {
    R = dst->rows;
    C = dst->cols;    

    if(IS_PRIM(dst)) {
      /* just copy array */
      R *= C;
      data = dst->data;
      for(i = 0; i < R; i++) 
        *data++ = *src++;         
    } else {
      /* copy for each row and column */
      for(i = 0; i < R; i++) {
        for(j = 0; j < C; j++) 
          *tm_at(dst,i,j) = *src++;          
      }
    }
  } else 
    e = TM_ERR_EMPTY_ARGS;
   
  if(err) *err = e; 
   
  return !e;    
}

tmVal tm_get(tMat* m, tmSize r, tmSize c, int* err)
{
  int e = 0;
  tmVal res = 0;
  
  if(m) {
    if(r < m->rows && c < m->cols)
      res = *tm_at(m,r,c);
    else 
      e = TM_ERR_WRONG_SIZE;     
  } else 
    e = TM_ERR_EMPTY_ARGS;
      
  if(err) *err = e;
   
  return res;
}

void tm_set(tMat* m, tmSize r, tmSize c, tmVal v, int* err)
{
  int e = 0;
   
  if(m) {
    if(r < m->rows && c < m->cols) 
      *tm_at(m,r,c) = v;
    else 
      e = TM_ERR_WRONG_SIZE;    
  } else 
    e = TM_ERR_EMPTY_ARGS;
   
  if(err) *err = e;  
}

tmSize tm_rows(tMat* m) 
{ 
  return m ? m->rows : 0;
}

tmSize tm_cols(tMat* m) 
{
  return m ? m->cols : 0;
}

tMat tm_copy(tMat* src, int* err)
{
  tmVal* data = NULL, *p;
  int e = 0, i, j;
  tMat res = NULL_TMATRIX;
   
  if(src) {
    res.rows = src->rows;
    res.cols = src->cols;
    res.width = res.cols;
    j = res.rows*res.cols;
    if(j) {
      data = (tmVal*) malloc(j * sizeof(tmVal));
      if(data) {
        res.data = data;
        res.type = TM_MAIN;        
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
    } else
      e = TM_ERR_WRONG_SIZE;
  } else 
    e = TM_ERR_EMPTY_ARGS;
      
  if(err) *err = e;
   
  return res;
}

void tm_print(tMat *m)
{
  int i, j, R, C;
  if(!m) {
    printf("No matrix!");  
    return;    
  }
  R = m->rows;
  C = m->cols;

  for(i = 0; i < R; i++) {
    printf("[ ");
    for(j = 0; j < C; j++) {
      printf("%.3f ", *tm_at(m, i, j));
    }
    printf("]\n");
  }

  switch(m->type) {
  case TM_STATIC:
    puts("(static)");
    break;
  case TM_TRANSPOSE:
    puts("(transposed)");
    break;
  case TM_SUB:
    puts("(submatrix)");
    break;  
  default:
    break;
  }
}

tMat tm_T(tMat* src, int *err) 
{
  tMat res = NULL_TMATRIX;
  int e = 0;    
  
  if(src) {
    if(IS_PRIM(src)) {
      res = *src;
      res.cols = src->rows;
      res.rows = src->cols;
      res.type = TM_TRANSPOSE;
    } else 
      e = TM_ERR_NOT_MAIN;
  } else 
    e = TM_ERR_EMPTY_ARGS;
   
  if(err) *err = e;
   
  return res;
}

tMat tm_block(tMat* src, tmSize r0, tmSize c0, tmSize Nr, tmSize Nc, int *err)
{
  tMat res = NULL_TMATRIX;
  int e = 0;
   
  if(src) {
    if(IS_PRIM(src)) {
      if(Nr > 0 && (r0+Nr) <= src->rows && Nc > 0 && (c0+Nc) <= src->cols) {
        res = *src;
        res.type = TM_SUB;
        res.rows = Nr;
        res.cols = Nc;
        res.data = src->data + (r0*src->cols + c0);
      } else 
        e = TM_ERR_WRONG_SIZE;
   } else 
     e = TM_ERR_NOT_MAIN;
  } else 
    e = TM_ERR_EMPTY_ARGS;

  if(err) *err = e;  

  return res;
}

int tm_insert(tMat *dst, tMat *src, int* err) 
{
  int e = 0, i,j,R,C;
   
  if(dst && src) {
    R = dst->rows;
    C = dst->cols;
    if(R == src->rows && C == src->cols) {
      for(i = 0; i < R; i++) {
        for(j = 0; j < C; j++) 
          *tm_at(dst,i,j) = *tm_at(src,i,j);            
      }
    } else 
      e = TM_ERR_DIFF_SIZE;
  } else 
    e = TM_ERR_EMPTY_ARGS;
   
  if(err) *err = e;  
   
  return !e; 
}

const char* tm_error(int code)
{
  static char *list[TM_ERR_TOTAL] = {
    "wrong size definition",
    "cannot allocate memory",
    "obligatory arguments are missed",
    "matrix with statically or dynamically allocated memory is expected",
    "source and destination have different size",
    "matrices are not compatible for current operation",
    "operation is not defined",
    "solution cannot be found",
    "homogeneous matrix is expected",
    "vector is expected"
  };
  
  if(code == 0)
    return "No errors";
  else if(code >= TM_ERR_TOTAL || code < 0)
    return "Wrong error code";
    
  return list[--code];
}
