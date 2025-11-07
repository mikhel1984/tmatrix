

#include <stdlib.h>
#include <math.h>
#include "tmatrix.h"
#include "tmatrix_priv.h"


int tf_chol(tMat* dst, tMat* m, int* err)
{
  int e = 0, n = 0, i, j, k, s, pos_def = 1;

  TM_ASSERT_ARGS(m && dst && m != dst, e, end_chol);

  if (m->rows == m->cols) {
    n = m->rows;
    if (tm_relevant(dst,n,n,&e)) {
      /* clear */
      for (i = 0; i < n; i++) {
        for (j = 0; j <= n; j++) 
          *tm_at(dst,i,j) = 0;
      }
      /* find lower left part */
      for (i = 0; i < n; i++) {
        for (j = 0; j < i; j++) {
          s = 0;
          for (k = 0; k < j; k++) {
            s += (*tm_at(dst,j,k)) * (*tm_at(dst,i,k));
          }
          s = *tm_at(m,i,j) - s;
          if (j < i) {
            *tm_at(dst,i,j) = s / (*tm_at(dst,j,i));
          } else {
            if (s < 0) {
              pos_def = 0;
              goto end_chol;
            }
            *tm_at(dst,i,j) = sqrt(s);
          }
        }
      }
    }
  } else 
    e = TM_ERR_NOT_DEF;
  
end_chol:
  if(err) *err = e;

  return pos_def;
}
