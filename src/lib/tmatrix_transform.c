

#include <stdlib.h>
#include <math.h>
#include "tmatrix.h"
#include "tmatrix_priv.h"

int tf_chol(tMat* dst, tMat* m, int* err)
{
  int e = 0, pos_def = 1, i, j, k, n;
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
              pos_def = 0; goto end_chol;
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

  return pos_def;
}
