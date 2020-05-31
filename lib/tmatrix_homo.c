/**
 * @file tmatrix_homo.c
 * @author Stanislav Mikhel
 * @date 2019
 * @brief Homogeneous transformations and specific simplifications.
 */ 
#include <stdlib.h>
#include <math.h>
#include "tmatrix_homo.h"
#include "tmatrix_priv.h"

#define HOMO_SIDE  4
#define HOMO_MASK {1,0,0,0, 0,1,0,0, 0,0,1,0, 0,0,0,1}

int is_homogenous(tMat *m, int *err)
{
  int e = 0;

  TM_ASSERT_ARGS(m, e, end_homogenous);

  if(IS_PRIM(m)) {
    if(m->rows != HOMO_SIDE || m->cols != HOMO_SIDE) 
      e = TM_ERR_NOT_HOMO;
  } else 
    e = TM_ERR_NOT_MAIN;
  
end_homogenous:
  if(err) *err = e;
  return !e;
}

int copy_data(tMat* m, tmVal v[])
{
  int i;
  tmVal *dst = m->data;
  for(i = 0; i < TM_HOMO_SIZE; i++)
    *dst++ = *v++;

  return 1;
}

int h_Txyz(tMat* dst, tmVal x, tmVal y, tmVal z, int *err) 
{
  tmVal arr[TM_HOMO_SIZE] = HOMO_MASK;
  arr[3]  = x;
  arr[7]  = y;
  arr[11] = z;

  return is_homogenous(dst, err) && copy_data(dst, arr);
}

int h_Rz(tMat* dst, tmVal a, int *err) 
{
  tmVal c, s;
  tmVal arr[TM_HOMO_SIZE] = HOMO_MASK;
  c = cos(a); s = sin(a);
  arr[0] = c; arr[1] = -s;
  arr[4] = s; arr[5] = c;

  return is_homogenous(dst, err) && copy_data(dst, arr);
}

int h_Ry(tMat* dst, tmVal b, int *err) 
{
  tmVal c, s;
  tmVal arr[TM_HOMO_SIZE] = HOMO_MASK;
  c = cos(b); s = sin(b);
  arr[0] = c; arr[2] = s;
  arr[8] = -s; arr[10] = c;

  return is_homogenous(dst, err) && copy_data(dst, arr);
}

int h_Rx(tMat* dst, tmVal v, int *err) 
{
  tmVal c, s;
  tmVal arr[TM_HOMO_SIZE] = HOMO_MASK;
  c = cos(v); s = sin(v);
  arr[5] = c; arr[6] = -s;
  arr[9] = s; arr[10] = c;

  return is_homogenous(dst, err) && copy_data(dst, arr);
}

int h_DH(tMat* dst, tmVal a, tmVal alpha, tmVal d, tmVal theta, int *err) 
{
  tmVal ca, sa, cth, sth;
  tmVal arr[TM_HOMO_SIZE] = HOMO_MASK;
  ca = cos(alpha); sa = sin(alpha);
  cth = cos(theta); sth = sin(theta);
  arr[0] = cth; arr[1] = -sth*ca; arr[2]  =  sth*sa; arr[3] = a*cth;
  arr[4] = sth; arr[5] =  cth*ca; arr[6]  = -cth*sa; arr[7] = a*sth;
                arr[9] = sa;      arr[10] = ca;     arr[11] = d;

  return is_homogenous(dst, err) && copy_data(dst, arr);
}

int h_mul(tMat *dst, tMat *m, int* err) 
{
  tmVal arr[TM_HOMO_SIZE] = HOMO_MASK;
  tmVal *a, *b;
  
  if(is_homogenous(dst,err) && is_homogenous(m,err)) {
    a = dst->data; b = m->data;
    /* row 1 */
    arr[0] = a[0]*b[0]+a[1]*b[4]+a[2]*b[8];
    arr[1] = a[0]*b[1]+a[1]*b[5]+a[2]*b[9];
    arr[2] = a[0]*b[2]+a[1]*b[6]+a[2]*b[10];
    arr[3] = a[0]*b[3]+a[1]*b[7]+a[2]*b[11]+a[3];
    /* row 2 */
    arr[4] = a[4]*b[0]+a[5]*b[4]+a[6]*b[8];
    arr[5] = a[4]*b[1]+a[5]*b[5]+a[6]*b[9];
    arr[6] = a[4]*b[2]+a[5]*b[6]+a[6]*b[10];
    arr[7] = a[4]*b[3]+a[5]*b[7]+a[6]*b[11]+a[7];
    /* row 3 */
    arr[8] = a[8]*b[0]+a[9]*b[4]+a[10]*b[8];
    arr[9] = a[8]*b[1]+a[9]*b[5]+a[10]*b[9];
    arr[10] = a[8]*b[2]+a[9]*b[6]+a[10]*b[10];
    arr[11] = a[8]*b[3]+a[9]*b[7]+a[10]*b[11]+a[11];
  } else
    return 0;
  
  return copy_data(dst, arr);   
}

int h_inv(tMat *dst, int *err)
{
  tmVal *a;
  tmVal arr[TM_HOMO_SIZE] = HOMO_MASK; 
  if(is_homogenous(dst,err)) {
    a = dst->data;
    /* row 1 */
    arr[0] = a[0]; arr[1] = a[4]; arr[2] = a[8];
    arr[3] = -(a[0]*a[3]+a[4]*a[7]+a[8]*a[11]);
    /* row 2 */
    arr[4] = a[1]; arr[5] = a[5]; arr[6] = a[9];
    arr[7] = -(a[1]*a[3]+a[5]*a[7]+a[9]*a[11]);
    /* row 3 */
    arr[8] = a[2]; arr[9] = a[6]; arr[10] = a[10];
    arr[11] = -(a[2]*a[3]+a[6]*a[7]+a[10]*a[11]);
  } else 
    return 0;
    
  return copy_data(dst,arr);
}

int h_T(tMat *dst, int *err)
{
  tmVal *a;
  tmVal arr[TM_HOMO_SIZE]; 
  if(is_homogenous(dst,err)) {
    a = dst->data;
    arr[0]  = a[0];  arr[1] = a[4];  arr[2] = a[8];   arr[3] = a[12];
    arr[4]  = a[1];  arr[5] = a[5];  arr[6] = a[9];   arr[7] = a[13];
    arr[8]  = a[2];  arr[9] = a[6]; arr[10] = a[10]; arr[11] = a[14];
    arr[12] = a[3]; arr[13] = a[7]; arr[14] = a[11]; arr[15] = a[15];
  } else
    return 0;
    
  return copy_data(dst,arr);
}

int h_eye(tMat *dst, int *err)
{
  static tmVal arr[TM_HOMO_SIZE] = HOMO_MASK;
  
  return is_homogenous(dst,err) && copy_data(dst,arr);
}
