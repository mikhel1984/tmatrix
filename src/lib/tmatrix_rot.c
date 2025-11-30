/**
 * @file tmatrix_rot.c
 * @author Stanislav Mikhel.
 * @date 2022
 * @brief Rotation matrices and quaternions.
 */
#include "tmatrix_rot.h"
#include "tmatrix_priv.h"
#include "tmatrix_vec.h"
//#define _USE_MATH_DEFINES
#include <math.h>

#define MAX_ERR 1E-6

#define qn_dot_(a,b) (a->w*b->w + a->x*b->x + a->y*b->y + a->z*b->z)

/* Find sign (-1/0/1) */
int sign_(tmVal v)
{
  return (v > 0) ?  1:
         (v < 0) ? -1: 0;
}

/* Find rotation matrix inversion */
int rot_inv(tMat* dst, tMat* m, int* err)
{
  int e = 0, i, j, k;
  TM_ASSERT_ARGS(dst && m, e, end_rinv);
  TM_ASSERT_INDEX(m->rows == 3 && m->cols == 3, e, end_rinv);

  if(tm_relevant(dst, 3, 3, &e)) {
    /* copy transposed values */
    if(IS_PRIM(dst)) {
      k = 0;
      for(i = 0; i < 3; ++i) {
        for(j = 0; j < 3; ++j)
          dst->data[k++] = *tm_at(m, j, i);
      }
    } else {
      for(i = 0; i < 3; ++i) {
        for(j = 0; j < 3; ++j)
          *tm_at(dst,i,j) = *tm_at(m,j,i);
      }
    }
  }

end_rinv:
  if(err) *err = e;

  return !e;
}

/* Get matrix from roll-pitch-yaw angles */
int rot_rpy(tMat* m, tmVal roll, tmVal pitch, tmVal yaw, int *err)
{
  int e = 0;
  tmVal cr, sr, cp, sp, cy, sy;

  TM_ASSERT_ARGS(m, e, end_rpy);

  if(tm_relevant(m, 3, 3, &e)) {
    sr = sin(roll);  cr = cos(roll);
    sp = sin(pitch); cp = cos(pitch);
    sy = sin(yaw);   cy = cos(yaw);

    if(IS_PRIM(m)) {
      m->data[0] = cp*cy; m->data[1] = sp*sr*cy-cr*sy; m->data[2] = sr*sy+sp*cr*cy;
      m->data[3] = cp*sy; m->data[4] = sp*sr*sy+cr*cy; m->data[5] = sp*cr*sy-sr*cy;
      m->data[6] = -sp;   m->data[7] = cp*sr;          m->data[8] = cp*cr;
    } else {
      *tm_at(m,0,0) = cp*cy; *tm_at(m,0,1) = sp*sr*cy-cr*sy; *tm_at(m,0,2) = sr*sy+sp*cr*cy;
      *tm_at(m,1,0) = cp*sy; *tm_at(m,1,1) = sp*sr*sy+cr*cy; *tm_at(m,1,2) = sp*cr*sy-sr*cy;
      *tm_at(m,2,0) = -sp;   *tm_at(m,2,1) = cp*sr;          *tm_at(m,2,2) = cp*cr;
    }
  }

end_rpy:
  if(err) *err = e;

  return !e;
}

/* Get matrix from axis and angle of rotation */
int rot_aa(tMat* m, tMat* k, tmVal a, int* err)
{
  int e = 0;
  tmVal c, s, v, kx, ky, kz;
  tmVal kxx, kxyv, kxzv, kyy, kyzv, kzz;

  TM_ASSERT_ARGS(m && k, e, end_aa);
  TM_ASSERT_INDEX(vec_len(k) == 3, e, end_aa);

  if(tm_relevant(m, 3, 3, &e)) {
    kx = vec_get(k, 0, &e);
    ky = vec_get(k, 1, &e);
    kz = vec_get(k, 2, &e);
    kxx = kx*kx; kyy = ky*ky; kzz = kz*kz;

    /* unit vector is expected */
    if(fabs(kxx + kyy + kzz - 1) < MAX_ERR) {
      c = cos(a); s = sin(a);
      v = 1 - c;
      kxyv = kx*ky*v; kxzv = kx*kz*v; kyzv = ky*kz*v;
      kx *= s; ky *= s; kz *= s;

      if(IS_PRIM(m)) {
        m->data[0] = kxx*v+c; m->data[1] = kxyv-kz; m->data[2] = kxzv+ky;
        m->data[3] = kxyv+kz; m->data[4] = kyy*v+c; m->data[5] = kyzv-kx;
        m->data[6] = kxzv-ky; m->data[7] = kyzv+kx; m->data[8] = kzz*v+c;
      } else {
        *tm_at(m,0,0) = kxx*v+c; *tm_at(m,0,1) = kxyv-kz; *tm_at(m,0,2) = kxzv+ky;
        *tm_at(m,1,0) = kxyv+kz; *tm_at(m,1,1) = kyy*v+c; *tm_at(m,1,2) = kyzv-kx;
        *tm_at(m,2,0) = kxzv-ky; *tm_at(m,2,1) = kyzv+kx; *tm_at(m,2,2) = kzz*v+c;
      }
    } else {
      e = TM_ERR_NOT_DEF;
    }
  }

end_aa:
  if(err) *err = e;

  return !e;
}

/* Get matrix from unit quaternion */
int rot_qn(tMat* m, tQn* q, int* err)
{
  int e = 0;
  tmVal xx, xy, xz, xw, yy, yz, yw, zz, zw, ww;

  TM_ASSERT_ARGS(m && q, e, end_qn);

  if(tm_relevant(m, 3, 3, &e)) {
    ww = q->w * q->w;
    xx = q->x * q->x;
    yy = q->y * q->y;
    zz = q->z * q->z;

    /* unit quaternion is expected */
    if(fabs(ww + xx + yy + zz - 1) < MAX_ERR) {
      xy = q->x * q->y;
      xz = q->x * q->z;
      xw = q->x * q->w;
      yz = q->y * q->z;
      yw = q->y * q->w;
      zw = q->z * q->w;

      if(IS_PRIM(m)) {
        m->data[0] = 1-2*(yy+zz); m->data[1] = 2*(xy-zw);   m->data[2] = 2*(xz+yw);
        m->data[3] = 2*(xy+zw);   m->data[4] = 1-2*(xx+zz); m->data[5] = 2*(yz-xw);
        m->data[6] = 2*(xz-yw);   m->data[7] = 2*(yz+xw);   m->data[8] = 1-2*(xx+yy);
      } else {
        *tm_at(m,0,0) = 1-2*(yy+zz); *tm_at(m,0,1) = 2*(xy-zw);   *tm_at(m,0,2) = 2*(xz+yw);
        *tm_at(m,1,0) = 2*(xy+zw);   *tm_at(m,1,1) = 1-2*(xx+zz); *tm_at(m,1,2) = 2*(yz-xw);
        *tm_at(m,2,0) = 2*(xz-yw);   *tm_at(m,2,1) = 2*(yz+xw);   *tm_at(m,2,2) = 1-2*(xx+yy);
      }
    } else {
      e = TM_ERR_NOT_DEF;
    }
  }

end_qn:
  if(err) *err = e;

  return !e;
}

/* Get unit quaternion from rotation matrix */
int rot_toqn(tQn* q, tMat* r, int* err)
{
  int e = 0;
  tmVal v, xx, yy, zz;

  TM_ASSERT_ARGS(q && r, e, end_toqn);
  TM_ASSERT_INDEX(r->rows == 3 && r->cols == 3, e, end_toqn);

  xx = *tm_at(r,0,0);
  yy = *tm_at(r,1,1);
  zz = *tm_at(r,2,2);
  if(xx + yy + zz > 0) {
    v = 2 * sqrt(1 + xx + yy + zz);
    q->w = 0.25 * v;
    q->x = ((*tm_at(r,2,1)) - (*tm_at(r,1,2))) / v;
    q->y = ((*tm_at(r,0,2)) - (*tm_at(r,2,0))) / v;
    q->z = ((*tm_at(r,1,0)) - (*tm_at(r,0,1))) / v;
  } else if(xx > yy && xx > zz) {
    v = 2 * sqrt(1 + xx - yy - zz);
    q->w = ((*tm_at(r,2,1)) - (*tm_at(r,1,2))) / v;
    q->x = 0.25 * v;
    q->y = ((*tm_at(r,0,1)) + (*tm_at(r,1,0))) / v;
    q->z = ((*tm_at(r,0,2)) + (*tm_at(r,2,0))) / v;
  } else if(yy > zz) {
    v = 2 * sqrt(1 + yy - xx - zz);
    q->w = ((*tm_at(r,0,2)) - (*tm_at(r,2,0))) / v;
    q->x = ((*tm_at(r,0,1)) + (*tm_at(r,1,0))) / v;
    q->y = 0.25 * v;
    q->z = ((*tm_at(r,2,1)) + (*tm_at(r,1,2))) / v;
  } else {
    v = 2 * sqrt(1 + zz - xx - yy);
    q->w = ((*tm_at(r,2,1)) - (*tm_at(r,1,2))) / v;
    q->x = ((*tm_at(r,0,2)) + (*tm_at(r,2,0))) / v;
    q->y = ((*tm_at(r,2,1)) + (*tm_at(r,1,2))) / v;
    q->z = 0.25 * v;
  }

end_toqn:
  if(err) *err = e;

  return !e;
}

/* Get roll-pitch-yaw from rotation matrix */
int rot_torpy(tmVal* roll, tmVal* pitch, tmVal* yaw, tMat* r, int *err)
{
  int e = 0;
  tmVal s;

  TM_ASSERT_ARGS(roll && pitch && yaw && r, e, end_torpy);
  TM_ASSERT_INDEX(r->rows == 3 && r->cols == 3, e, end_torpy);

  s = -(*tm_at(r,2,0));
  if(fabs(s-1) > MAX_ERR && fabs(s+1) > MAX_ERR) {
    /* s != 1 && s != -1 */
    *roll  = atan2(*tm_at(r,2,1), *tm_at(r,2,2));
    *pitch = asin(s);
    *yaw   = atan2(*tm_at(r,1,0), *tm_at(r,0,0));
  } else {
    /* s == -1 or s == 1, can find yaw +/- roll, set yaw */
    *roll = 0;
    *pitch = (s > 0) ? M_PI_2 : (-M_PI_2);
    *yaw = atan2(*tm_at(r,1,2), *tm_at(r,1,1));
    e = TM_ERR_NO_SOLUTN;
  }

end_torpy:
  if(err) *err = e;

  return !e;
}

/* Get axis and angle from rotation matrix */
int rot_toaa(tMat* k, tmVal* a, tMat* r, int* err)
{
  int e = 0;
  tmVal s, c, p0, p1, p2;

  TM_ASSERT_ARGS(k && a && r, e, end_toaa);
  TM_ASSERT_INDEX(r->rows == 3 && r->cols == 3, e, end_toaa);

  if(tm_relevant(k,3,1,&e)) {
    c = ((*tm_at(r,0,0)) + (*tm_at(r,1,1)) + (*tm_at(r,2,2)) - 1) * 0.5;

    if(fabs(c-1) > MAX_ERR && fabs(c+1) > MAX_ERR) {
      p0 = (*tm_at(r,2,1)) - (*tm_at(r,1,2));
      p1 = (*tm_at(r,0,2)) - (*tm_at(r,2,0));
      p2 = (*tm_at(r,1,0)) - (*tm_at(r,0,1));
      s = sqrt(p0*p0 + p1*p1 + p2*p2) * 0.5;
      *a = atan2(s, c);

      if(c >= 0) {
        s *= 2;
        *tm_at(k,0,0) = p0 / s;
        *tm_at(k,1,0) = p1 / s;
        *tm_at(k,2,0) = p2 / s;
      } else {
        s = 1 - c;
        *tm_at(k,0,0) = sign_(p0) * sqrt((*tm_at(r,0,0)-c)/s);
        *tm_at(k,1,0) = sign_(p1) * sqrt((*tm_at(r,1,1)-c)/s);
        *tm_at(k,2,0) = sign_(p2) * sqrt((*tm_at(r,2,2)-c)/s);
      }
    } else {
      e = TM_ERR_NO_SOLUTN;
      *a = 0;
      /* arbitrary unit axis */
      *tm_at(k,0,0) = 0;
      *tm_at(k,1,0) = 0;
      *tm_at(k,2,0) = 1;
    }
  }

end_toaa:
  if(err) *err = e;

  return !e;
}

/* q1 + q2 */
tQn qn_add(tQn* q1, tQn* q2, int* err)
{
  int e = 0;
  tQn sum = UNIT_QTN;

  TM_ASSERT_ARGS(q1 && q2, e, end_qadd);

  sum.w = q1->w + q2->w;
  sum.x = q1->x + q2->x;
  sum.y = q1->y + q2->y;
  sum.z = q1->z + q2->z;

end_qadd:
  if(err) *err = e;

  return sum;
}

/* q1 - q2 */
tQn qn_sub(tQn* q1, tQn* q2, int* err)
{
  int e = 0;
  tQn sum = UNIT_QTN;

  TM_ASSERT_ARGS(q1 && q2, e, end_qsub);

  sum.w = q1->w - q2->w;
  sum.x = q1->x - q2->x;
  sum.y = q1->y - q2->y;
  sum.z = q1->z - q2->z;

end_qsub:
  if(err) *err = e;

  return sum;
}

/* q1 * q2 */
tQn qn_mul(tQn* q1, tQn* q2, int* err)
{
  int e = 0;
  tQn res = UNIT_QTN;

  TM_ASSERT_ARGS(q1 && q2, e, end_qmul);

  res.w = q1->w*q2->w - q1->x*q2->x - q1->y*q2->y - q1->z*q2->z;
  res.x = q1->w*q2->x + q1->x*q2->w + q1->y*q2->z - q1->z*q2->y;
  res.y = q1->w*q2->y - q1->x*q2->z + q1->y*q2->w + q1->z*q2->x;
  res.z = q1->w*q2->z + q1->x*q2->y - q1->y*q2->x + q1->z*q2->w;

end_qmul:
  if(err) *err = e;

  return res;
}

/* Quaternion conjugation */
tQn qn_conj(tQn* q, int* err)
{
  int e = 0;
  tQn res = UNIT_QTN;

  TM_ASSERT_ARGS(q, e, end_conj);

  res.w = q->w;
  res.x = -q->x;
  res.y = -q->y;
  res.z = -q->z;

end_conj:
  if(err) *err = e;

  return res;
}

/* 1 / q */
tQn qn_inv(tQn* q, int* err)
{
  int e = 0;
  tmVal n;
  tQn res = UNIT_QTN;

  TM_ASSERT_ARGS(q, e, end_qinv);

  n = qn_dot_(q, q);
  if(n > MAX_ERR * MAX_ERR) {
    res.w = q->w / n;
    res.x = -q->x / n;
    res.y = -q->y / n;
    res.z = -q->z / n;
  } else {
    e = TM_ERR_NO_SOLUTN;
  }

end_qinv:
  if(err) *err = e;

  return res;
}

/* Obtain unit quaternion */
int qn_normalize(tQn* q)
{
  if(q) {
    tmVal n = sqrt(qn_dot_(q, q));
    if(n > MAX_ERR * MAX_ERR) {
      q->w /= n;
      q->x /= n;
      q->y /= n;
      q->z /= n;
      return 1;
    }
  }
  return 0;
}

/* Absolute value */
tmVal qn_abs(tQn* q)
{
  return q ? sqrt(qn_dot_(q,q)) : 0;
}

/* Product with scalar */
tQn qn_scale(tQn* q, tmVal k, int* err)
{
  int e = 0;
  tQn res = UNIT_QTN;

  TM_ASSERT_ARGS(q, e, end_qscale);

  res.w = q->w * k;
  res.x = q->x * k;
  res.y = q->y * k;
  res.z = q->z * k;

end_qscale:
  if(err) *err = e;

  return res;
}

/* Interpolation between two unit quaternions */
tQn qn_slerp(tQn* q1, tQn* q2, tmVal t, int *err)
{
  int e = 0;
  tmVal fTheta, theta, t1 = 0.5, t2 = 0.5;
  tQn res = UNIT_QTN;

  TM_ASSERT_ARGS(q1 && q2, e, end_slerp);
  TM_ASSERT_INDEX(t >= 0 && t <= 1, e, end_slerp);

  /* cos(theta/2) */
  fTheta = qn_dot_(q1, q2);
  if(fabs(fTheta-1) > MAX_ERR && fabs(fTheta+1) > MAX_ERR) {
    theta = acos(fTheta);
    fTheta = sin(theta);
    if(fabs(fTheta) > MAX_ERR) {
      t1 = sin((1-t)*theta) / fTheta;
      t2 = sin(t*theta) / fTheta;
    }

    res.w = t1*q1->w + t2*q2->w;
    res.x = t1*q1->x + t2*q2->x;
    res.y = t1*q1->y + t2*q2->y;
    res.z = t1*q1->z + t2*q2->z;
  } else {
    /* return q1 */
    res = *q1;
  }

end_slerp:
  if(err) *err = e;

  return res;
}


