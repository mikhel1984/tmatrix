/**
 * @file tmatrix_print.c
 * @author Stanislav Mikhel
 * @date 2020
 * @brief Visualization for matrices and errors
 */ 
#include <stdio.h>
#include "tmatrix_io.h"
#include "tmatrix_priv.h"

/* Simple matrix visualization */
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

/* Error desctiption */
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
    return "Unknown error code";
    
  return list[--code];
}
