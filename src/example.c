/*	example.h

  Consider an example of library application for creating forward
  kinematics and Jacobean functions for "typical" manipulator structure 
  Rz - Tz - Ry - Tz - Ry - Tz - Rx - Tz - Ry - Tz - Rx - Tz
*/
#include <stdio.h>
#include "lib/tmatrix_homo.h"
#include "lib/tmatrix_vec.h"
#include "lib/tmatrix_io.h"

#define JOINT_NO 6   /* number of joints */
#define CART_NO  6   /* number of coordinates in Cartesian space */

/* Forward kinematics for current joint state */
tMat FK(tmVal q[], int *err);
/* Jacobean for current joint state */
tMat Jac(tmVal q[], int *err);
/* Simplify Jacobean evaluation */
int copy_cols(tMat *dst, tMat *src, int c, int axe, int *err);
/* Links */
tmVal d[JOINT_NO] = {0.3,0.4,0.5,0.5,0.2,0.1};

int main()
{
  int err = 0;
  /* joint state */
  tmVal q[JOINT_NO] = {0.1,-0.2,0.3,-0.4,0.5,-0.6};
  /* matrices */
  tMat fk = NULL_TMATRIX, jac = NULL_TMATRIX, 
       cvel = NULL_TMATRIX, qvel = NULL_TMATRIX, 
       ct;
  
  puts("Forward kinematics");
  fk = FK(q, &err);
  if(err) printf("ERROR: %s\n", tm_error(err));
  tm_print(&fk);

  puts("Jacobean");
  jac = Jac(q, &err);
  if(err) printf("ERROR: %s\n", tm_error(err));
  tm_print(&jac);
  
  /* assume that joint velocity is also q 
     and find Cartesian velocities */
  puts("Cartesian velocity");
  qvel = vec_static(JOINT_NO,q,&err);
  if(err) printf("ERROR: %s\n", tm_error(err));
  
  if(!tm_mul(&cvel,&jac,&qvel,&err)) printf("ERROR: %s\n", tm_error(err));
  /* transpose for visualization */
  ct = tm_T(&cvel,&err); 
  if(err) printf("ERROR: %s\n", tm_error(err));
  tm_print(&ct);  
  
  /* free memory */
  /* obligatory */
  tm_clear(&fk);
  tm_clear(&jac);
  tm_clear(&cvel);
  /* optional (do nothing) */
  tm_clear(&qvel);  /* static memory */
  tm_clear(&ct);    /* reference */
  
  return 0;
}

tMat FK(tmVal q[], int *err)
{
  tmVal a1[TM_HOMO_SIZE];  /* memory for intermediate operations */
  tMat tmp = NULL_TMATRIX, res = NULL_TMATRIX;

  *err = 0;
  tmp = h_static(a1,err); /* matrix for temporary states */
  if(*err) return res; 
  res = h_new(err);       /* allocate memory */
  if(*err) return res;

  h_Rz(&res,q[0],err);    /* initialize res as Z rotation matrix */
  if(*err) return res;

  /* transformations (in place) */
  /*   initialize     &&    multiply       */
  if(
  h_Tz(&tmp,d[0],err) && h_mul(&res,&tmp,err) &&
  h_Ry(&tmp,q[1],err) && h_mul(&res,&tmp,err) &&
  h_Tz(&tmp,d[1],err) && h_mul(&res,&tmp,err) &&
  h_Ry(&tmp,q[2],err) && h_mul(&res,&tmp,err) && 
  h_Tz(&tmp,d[2],err) && h_mul(&res,&tmp,err) &&
  h_Rx(&tmp,q[3],err) && h_mul(&res,&tmp,err) &&
  h_Tz(&tmp,d[3],err) && h_mul(&res,&tmp,err) &&
  h_Ry(&tmp,q[4],err) && h_mul(&res,&tmp,err) &&
  h_Tz(&tmp,d[4],err) && h_mul(&res,&tmp,err) &&
  h_Rx(&tmp,d[5],err) && h_mul(&res,&tmp,err) &&
  h_Tz(&tmp,d[5],err) && h_mul(&res,&tmp,err)) {}  /* just call sequential execution */

  return res;
}

int copy_cols(tMat *dst, tMat *src, int c, int axe, int *err)
{
  tMat sub1, sub2;
  
  /* position */
  sub1 = tm_block(dst, 0, c, 3, 1, err); if(*err) return 0;
  sub2 = tm_block(src, 0, 3, 3, 1, err); if(*err) return 0;
  if(!tm_insert(&sub1, &sub2, err)) return 0;
  /* axis */
  sub1 = tm_block(dst, 3, c, 3, 1, err); if(*err) return 0;
  sub2 = tm_block(src, 0, axe, 3, 1, err); if(*err) return 0;
  if(!tm_insert(&sub1, &sub2, err)) return 0;
  
  return 1;  
}

tMat Jac(tmVal q[], int *err)
{
  tmVal a1[TM_HOMO_SIZE], a2[TM_HOMO_SIZE]; /* local memory */
  tMat tmp = NULL_TMATRIX, acc = NULL_TMATRIX, res = NULL_TMATRIX,
       ee, col, axis;
  int success, i;

  *err = 0;
  /* intermediate states */
  tmp = h_static(a1,err); if(*err) return res; 
  acc = h_static(a2,err); if(*err) return res;
  /* allocate memory for Jacobean (6x6) */
  res = tm_new(CART_NO,JOINT_NO,err);
  if(*err) return res;
  
  /* prepare Jacobean columns */
  h_eye(&acc,err); if(*err) return res;
  /* joint 1 position and Z axis */
  if(!copy_cols(&res,&acc,0,2,err)) return res;
  success = 
  /* initialize            transform             save to Jacobean */
  /* joint 2 position and Y axis */  
  h_Rz(&tmp,q[0],err) && h_mul(&acc,&tmp,err) &&
  h_Tz(&tmp,d[0],err) && h_mul(&acc,&tmp,err) && copy_cols(&res,&acc,1,1,err) && 
  /* joint 3 position and Y axis */
  h_Ry(&tmp,q[1],err) && h_mul(&acc,&tmp,err) &&
  h_Tz(&tmp,d[1],err) && h_mul(&acc,&tmp,err) && copy_cols(&res,&acc,2,1,err) &&
  /* joint 4 position and X axis */
  h_Ry(&tmp,q[2],err) && h_mul(&acc,&tmp,err) &&
  h_Tz(&tmp,d[2],err) && h_mul(&acc,&tmp,err) && copy_cols(&res,&acc,3,0,err) &&
  /* joint 5 position and Y axis */
  h_Rx(&tmp,q[3],err) && h_mul(&acc,&tmp,err) && 
  h_Tz(&tmp,d[3],err) && h_mul(&acc,&tmp,err) && copy_cols(&res,&acc,4,1,err) &&
  /* joint 6 position and X axis */
  h_Ry(&tmp,q[4],err) && h_mul(&acc,&tmp,err) &&
  h_Tz(&tmp,d[4],err) && h_mul(&acc,&tmp,err) && copy_cols(&res,&acc,5,0,err) &&
  /* end-effector position */
  h_Rx(&tmp,q[5],err) && h_mul(&acc,&tmp,err) &&
  h_Tz(&tmp,d[5],err) && h_mul(&acc,&tmp,err); 
  if(!success) return res;
  
  /* recalculate linear velocity components */
  ee = tm_block(&acc,0,3,3,1,err); if(*err) return res; /* read from final state */
  tmp = vec_static(3,a1,err);      if(*err) return res; /* reuse variable */ 
  
  for(i = 0; i < JOINT_NO; i++) {    
    col  = tm_block(&res,0,i,3,1,err); if(*err) return res;
    axis = tm_block(&res,3,i,3,1,err); if(*err) return res;
    /* axis x (ee - col) */
    success = tm_insert(&tmp,&ee,err) && tm_sub(&tmp,&col,err) && vec_cross(&col,&axis,&tmp,err);
    if(!success) break;
  }
  
  return res;
}

