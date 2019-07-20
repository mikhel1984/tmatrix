
#include <stdlib.h>
#include "tmatrix_vec.h"
#include "tmatrix_priv.h"

tmSize vec_len(tMat *m)
{
  if(m) {
    if(m->cols == 1) {
      return m->rows;
    } else if(m->rows == 1) {
      return m->cols;
    }
  }
  return 0;
}

tmVal vec_get(tMat *m, tmSize k, int *err)
{
  int e = 0;
  tmVal res = 0;
  if(m) {
    if(m->cols == 1) {
      res = tm_get(m,k,0,&e);
    } else if(m->rows == 1) {
      res = tm_get(m,0,k,&e);
    } else 
      e = TM_ERR_NOT_VEC;
  } else 
    e = TM_ERR_EMPTY_ARGS;
    
  if(err) *err = e;
  
  return res;
}

void vec_set(tMat *m, tmSize k, tmVal v, int *err)
{
  int e = 0;
  
  if(m) {
    if(m->cols == 1) {
      tm_set(m,k,0,v,&e);
    } else if(m->rows == 1) {
      tm_set(m,0,k,v,&e);      
    } else 
      e = TM_ERR_NOT_VEC;
  } else 
    e = TM_ERR_EMPTY_ARGS;
    
  if(err) *err = e;
}

tmVal vec_dot(tMat *a, tMat *b, int *err)
{
  int e = 0,i,N1,N2;
  tmVal sum = 0;
  
  N1 = vec_len(a); N2 = vec_len(b);
  if(N1 && N2) {
    if(N1 == N2) {
      if(IS_PRIM(a) && IS_PRIM(b)) {
        for(i = 0; i < N1; i++) 
          sum += a->data[i] * b->data[i];        
      } else {
        for(i = 0; i < N1 && !e; i++)
          sum += vec_get(a,i,&e) * vec_get(b,i,&e);
      }
    } else 
      e = TM_ERR_NOT_COMPAT;
  } else 
    e = TM_ERR_NOT_VEC;
  
  if(err) *err = e;
  
  return sum;
}

int vec_cross(tMat *res, tMat *a, tMat *b, int *err)
{
  int e = 0;
  tmVal *data = 0, a0,a1,a2,b0,b1,b2;
  
  
  if(res && a && b) {
    if(vec_len(a) == 3 && vec_len(b) == 3) {
      if(vec_len(res) != 3) {
        if(IS_PRIM(res)) {
          if(res->rows * res->cols < 3) {
            if(res->type == TM_MAIN) {
              data = (tmVal*) malloc(3*sizeof(tmVal));
              if(data) {
                free(res->data);
                res->data = data;
              } else 
                e = TM_ERR_NO_MEMORY;
            } else 
              e = TM_ERR_NOT_MAIN;
            if(!e) {
              res->rows = 3;
              res->cols = 1;
              res->width = 1;
            } else {
              if(err) *err = e;
              return !e;
            }
          }          
        } else 
          e = TM_ERR_NOT_MAIN; 
      }
      a0 = vec_get(a,0,0); a1 = vec_get(a,1,0); a2 = vec_get(a,2,0);
      b0 = vec_get(b,0,0); b1 = vec_get(b,1,0); b2 = vec_get(b,2,0);
      if(!e) vec_set(res,0,a1*b2-b1*a2,&e);
      if(!e) vec_set(res,1,a2*b0-a0*b2,&e);
      if(!e) vec_set(res,2,a0*b1-b0*a1,&e);      
    } else
      e = TM_ERR_NOT_DEF;
  } else 
    e = TM_ERR_EMPTY_ARGS;
    
  if(err) *err = e;
  
  return !e;  
}

