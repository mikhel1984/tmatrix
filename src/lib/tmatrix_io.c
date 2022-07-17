/**
 * @file tmatrix_print.c
 * @author Stanislav Mikhel
 * @date 2020
 * @brief Visualization for matrices and errors
 */ 
#include <stdio.h>
#include <stdlib.h>
#include "tmatrix_io.h"
#include "tmatrix_priv.h"

#define BUFF_SIZE 80

/* Simple matrix visualization */
void tm_print(tMat *m)
{
  int i, j, R, C;
  if(!m) {
    puts("No matrix!");  
    return;    
  }
  R = m->rows;
  C = m->cols;

  for(i = 0; i < R; i++) {
    printf("[ ");
    for(j = 0; j < C; j++) {
      printf("%G ", *tm_at(m, i, j));
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

/* Matrix to file in csv-like form.*/
int tm_to_file(tMat* m, char* fname, char sep)
{
  int r, c, r1, c1; 
  FILE* pFile = 0;

  if(!m || !fname) return 0; 
  if((pFile = fopen(fname, "w")) == 0) return 0;

  r1 = m->rows-1; 
  c1 = m->cols-1;
  for(r = 0; r <= r1 ; r++) {
    for(c = 0; c <= c1 ; c++) {
      fprintf(pFile, "%G", *tm_at(m,r,c));
      if(c != c1) fputc(sep, pFile);
    }
    if(r != r1) fputc('\n', pFile);
  }

  fclose(pFile);

  return 1;
}

/* Initialize matrix from file */
int tm_from_file(tMat* dst, char* fname, char sep)
{
  int i = 0, j = 0, c;
  char buff[BUFF_SIZE], *pb;
  FILE* pFile = 0;

  if(!dst || !fname) return 0;
  if((pFile = fopen(fname,"r")) == 0) return 0;

  // estimate matrix size 
  while((c = fgetc(pFile)) != EOF) {
    if(c == sep) {
      j++; 
    } else if(c == '\n') {
      i++; j = 0; 
    }
  }

  // prepare matrix
  if(!tm_relevant(dst,++i,++j,0)) {
    fclose(pFile);
    return 0;
  }

  rewind(pFile); 
  i = 0; j = 0; 
  pb = buff;
  while((c = fgetc(pFile)) != EOF) {
    if(c == sep || c == '\n') {
      *pb = '\0';
      *tm_at(dst,i,j++) = (tmVal) strtod(buff, &pb);
      pb = buff;
      if(c == '\n') {
        i ++; j = 0;
      }
    } else {
      *pb++ = (char) c;
    }
  }
  *pb = '\0';
  *tm_at(dst,i,j++) = (tmVal) strtod(buff, &pb);

  fclose(pFile);
  return 1;
}

/* Show quaternion */
void qn_print(tQn* q)
{
  if(q) {
    printf("{w:%f x:%f y:%f z:%f}\n", q->w, q->x, q->y, q->z);
  } else {
    puts("No quaternion");
  }
}

