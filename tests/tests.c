
#include <stdio.h>
#include <math.h>
#include <time.h>
#include "../lib/tmatrix_homo.h"
#include "../lib/tmatrix_vec.h"
#include "../lib/tmatrix_print.h"
#include "minunit.h"

#define EQL(X,Y) fabs((X)-(Y)) < 1E-6

int tests_run = 0;
char buf[BUFSIZE];

/* ~~~~~~~~~~~~~~~~ Test list ~~~~~~~~~~~~~~~~~~~~~~ */

static char* all_tests();

static char* test_init();
static char* test_transpose();
static char* test_submatrix();
static char* test_cols();
static char* test_sum();
static char* test_prod();
static char* test_inv();
static char* test_pinv();
static char* test_homo();
static char* test_vec();
static char* test_rank();
static char* test_make();

/* ~~~~~~~~~~~~~~~~~~~ Main ~~~~~~~~~~~~~~~~~~~~~~~~ */

int main()
{  
  float t;
  clock_t beg = clock();
  char *result = all_tests();
  t = (clock()-beg) / ((float) CLOCKS_PER_SEC);
  if (result) 
    printf("ERROR %s\n", result);
  else 
    puts("ALL TESTS PASSED");
  
  printf("Tests No: %d\n", tests_run);
  printf("Duration: %f s\n", t);

  return result != 0;
}

/* ~~~~~~~~~~~~~~~~~~~ Tests ~~~~~~~~~~~~~~~~~~~~~~ */ 

static char* all_tests() 
{
  mu_run(test_init);
  mu_run(test_transpose);
  mu_run(test_submatrix);
  mu_run(test_cols);
  mu_run(test_sum);
  mu_run(test_prod);
  mu_run(test_inv);
  mu_run(test_pinv);
  mu_run(test_homo);
  mu_run(test_vec);
  mu_run(test_rank);
  mu_run(test_make);
  
  return 0;
}

static char* test_init()
{
  int err = 0, i;
  tmVal v;
  tmVal arr[9] = {
    0,1,2,
    3,4,5,
    6,7,0};
  tMat cp, m = tm_new(3,3,&err), m2;
  mu_check("Init (new):", err);
   
  tm_set(&m, 1, 1, 2, &err);
  mu_check("Init (set):", err);

  v = tm_get(&m, 1, 1, &err);
  mu_check("Init (get):", err);
  mu_assert("Init (get): v != 2", EQL(v,2));
   
  m2 = tm_static(3,3, arr, &err);
  mu_check("Init (static):", err);   
  for(i = 0; i < 9; i++) {
    mu_assert("Init (static): copy error", EQL(arr[i],m2.data[i])); 
  }
   
  cp = tm_copy(&m2,&err);
  mu_check("Init (copy):", err);
  for(i = 0; i < 9; i++) {
    mu_assert("Init (copy): wrong result", EQL(cp.data[i],m2.data[i])); 
  }
  
  tm_print(&m);

  tm_clear(&m);
  tm_clear(&cp);

  return 0;
}

static char* test_transpose() 
{
  int err = 0, i,j;
  tmVal v;
  tmVal arr[9] = {
    0,1,2,
    3,4,5,
    6,7,0};

  tMat t, cp, m = tm_static(3,3,arr,&err);
  mu_check("Transpose (static):", err);
  
  t = tm_T(&m, &err);
  mu_check("Transpose (T):", err);
  for(i = 0; i < 3; i++) {
    for(j = 0; j < 3; j++) {
      mu_assert("Transpose (T): wrong result", EQL(tm_get(&t,i,j,0), tm_get(&m,j,i,0)));
    }
  }

  v = tm_get(&t, 0, 1, &err);
  mu_check("Transpose (get):", err);
  mu_assert("Transpose (get): m(x,y) != m'(y,x)", EQL(v,tm_get(&m,1,0,NULL)));   

  tm_set(&t, 2, 1, 10, &err);
  mu_check("Transpose (set):", err);
  mu_assert("Transpose (set): m(x,y) != m'(y,x)", EQL(10,tm_get(&m,1,2,NULL)));  

  cp = tm_copy(&t,&err);
  mu_check("Transpose (copy):", err);
  for(i = 0; i < 3; i++) {
    for(j = 0; j < 3; j++) {
      mu_assert("Transpose (copy): wrong result", EQL(tm_get(&t,i,j,0), tm_get(&cp,i,j,0)));
    }
  }

  tm_clear(&m);
  tm_clear(&t);
  tm_clear(&cp);
  
  return 0;
}

static char* test_submatrix()
{
  int err = 0,i,j;
  tmVal v;
  tmVal arr[9] = {
    0,1,2,
    3,4,5,
    6,7,0};

  tMat s, cp, m = tm_static(3,3,arr,&err);
  mu_check("Submatrix (static):", err);

  s = tm_block(&m, 1,1,2,2, &err);
  mu_check("Submatrix (block):", err);
  
  v = tm_get(&s, 0, 1, &err);
  mu_check("Submatrix (get):", err);  
  for(i = 0; i < 2; i++) {
    for(j = 0; j < 2; j++) {
      mu_assert("Submatrix (get): not equal", EQL(tm_get(&m,i+1,j+1,0), tm_get(&s,i,j,0)));
    }
  }
   
  tm_set(&s, 1, 0, 10, &err);
  mu_check("Submatrix (set):", err);  
  v = tm_get(&s,1,0,0); 
  mu_assert("Submatrix (s): v != 10", EQL(v,10));
  
  cp = tm_copy(&s,&err);
  mu_check("Submatrix (copy):", err);
  for(i = 0; i < 2; i++) {
    for(j = 0; j < 2; j++) {
      mu_assert("Submatrix (copy): wrong result", EQL(tm_get(&cp,i,j,0), tm_get(&s,i,j,0)));
    }
  }
  

  tm_clear(&m);
  tm_clear(&s);
  tm_clear(&cp);
  
  return 0;
}

static char* test_cols()
{
  int err = 0,i;
  tmVal arr[9] = {
    0,1,2,
    3,4,5,
    6,7,0};
  tMat c1,c2,m = tm_static(3,3,arr,&err);
  mu_check("Cols (static):", err);
  
  c1 = tm_col(&m,1,&err);
  mu_check("Cols (col):", err);
  //c2 = tm_col(&m,3,&err);
  //mu_assert("Cols (col): wrong indexation", err != 0);
  c2 = tm_col(&m,2,&err);
  mu_check("Cols (col):", err);
  mu_assert("Cols (get): not equal", EQL(tm_get(&c1,1,0,0), tm_get(&m,1,1,0)));
  
  tm_insert(&c1,&c2,&err);
  mu_check("Cols (insert):", err);
  for(i = 0; i < 3; i++) {
      mu_assert("Cols (insert): wrong result", EQL(tm_get(&m,i,1,0), tm_get(&m,i,2,0)));    
  }
    
  return 0;  
}

static char* test_sum()
{
  int err = 0, i;
  tmVal a1[9] = {1,2,3,4,5,6,7,8,9};
  tmVal a2[9] = {3,5,7,2,4,8,0,6,9};
  tMat m1, m2, m0;
   
  m0 = tm_static(3,3,a1,&err);
  mu_check("Sum (static):", err);
  m1 = tm_copy(&m0,&err);
  mu_check("Sum (copy):", err);
  m2 = tm_static(3,3,a2,&err);
  mu_check("Sum (static):", err);
   
  tm_add(&m1,&m2,&err);
  mu_check("Sum (add):", err);
  for(i = 0; i < 9; i++) {
    mu_assert("Sum (add): wrong sum", EQL(m1.data[i], a1[i]+a2[i])); 
  }
   
  tm_sub(&m1,&m2,&err);
  mu_check("Sum (sub):", err);
  for(i = 0; i < 9; i++) {
    mu_assert("Sum (add): wrong sum", EQL(m1.data[i], a1[i])); 
  }
  
  tm_scale(&m1,2,&err);
  mu_check("Sum (scale):", err);
  for(i = 0; i < 9; i++) {
    mu_assert("Sum (add): wrong sum", EQL(m1.data[i], a1[i]*2)); 
  }
   
  tm_clear(&m1);
  tm_clear(&m2);
   
  return 0;   
}

static char* test_prod()
{
  int err = 0, i;
  tmVal a1[9] = {1,2,3,4,5,6,7,8,9};
  tmVal a2[9];
  tmVal a3[9] = {14,32,50,32,77,122,50,122,194};
  tMat m1, m2, m3;
   
  m1 = tm_static(3,3,a1,&err);
  mu_check("Prod (static):", err);
  m2 = tm_static(3,3,a2,&err);
  mu_check("Prod (static):", err);
  m3 = tm_T(&m1,&err);
  mu_check("Prod (transpose):", err);
   
  tm_mul(&m2, &m1,&m3, &err);
  mu_check("Prod (mul):", err);
  for(i = 0; i < 9; i++) {
    mu_assert("Prod (mul): wrong product", EQL(m2.data[i], a3[i]));
  }
    
  return 0;     
}

static char* test_inv()
{
  int err = 0,i,j;
  tmVal a1[9] = {1,3,4,5,7,2,8,6,1}, d;
  tMat m, im, pr; 
   
  m = tm_static(3,3,a1,&err);
  mu_check("Inv (static):", err);
   
  d = tm_det(&m, &err);
  mu_check("Inv (det):", err);
  mu_assert("Inv (det): wrong value", EQL(d,-76));
   
  im = tm_new(0,0,&err);
  tm_inv(&im, &m, &err);
  mu_check("Inv (det):", err);
  pr = tm_simp();
  
  tm_mul(&pr,&m,&im,&err); 
  mu_check("Inv (mul):", err);
  for(i = 0; i < 3; i++) {
    for(j = 0; j < 3; j++) {
      d = (i == j) ? 1 : 0;
      mu_assert("Inv (inv): wrong matrix", EQL(tm_get(&pr,i,j,0),d));
    }
  }
  
  tm_mul(&pr,&im,&m,&err); 
  mu_check("Inv (mul):", err);
  for(i = 0; i < 3; i++) {
    for(j = 0; j < 3; j++) {
      d = (i == j) ? 1 : 0;
      mu_assert("Inv (inv): wrong matrix", EQL(tm_get(&pr,i,j,0),d));
    }
  }
  
  tm_clear(&pr);   
  tm_clear(&m);
  tm_clear(&im);
   
  return 0;
}

static char* test_pinv()
{
  int err = 0,i,j;
  tmVal a1[8] = {1,3,4,5,7,2,8,6},v;
  tMat m, im, pr; 
   
  m = tm_static(2,4,a1,&err);
  mu_check("Pinv (static):", err);
   
  im = tm_new(0,0, &err);
  tm_pinv(&im, &m, &err);
  mu_check("Pinv (pinv):", err);
   
  pr = tm_simp();     
  tm_mul(&pr,&m,&im,&err);
  mu_check("Pinv (mul):", err);   
  for(i = 0; i < 2; i++) {
    for(j = 0; j < 2; j++) {
      v = (i == j) ? 1 : 0;
      mu_assert("Pinv (pinv): wrong matrix", EQL(tm_get(&pr,i,j,0),v));
    }
  }
  
  m = tm_static(4,2,a1,NULL);
  tm_pinv(&im, &m,&err);
  mu_check("Pinv (pinv):", err);
  tm_mul(&pr,&im,&m,&err);
  mu_check("Pinv (mul):", err);   
  for(i = 0; i < 2; i++) {
    for(j = 0; j < 2; j++) {
      v = (i == j) ? 1 : 0;
      mu_assert("Pinv (pinv): wrong matrix", EQL(tm_get(&pr,i,j,0),v));
    }
  }   
   
  tm_clear(&pr);
  tm_clear(&m);
  tm_clear(&im);
   
  return 0;
}

static char* test_homo()
{
  int err = 0, i, j;
  tmVal a1[TM_HOMO_SIZE], v;
  tMat m1, m2, cp, pr;
     
  m1 = h_static(a1,&err);
  mu_check("Homo (static):", err);   
  h_Tx(&m1,10,&err);
  mu_check("Homo (Tx):", err);  
      
  m2 = h_new(&err);
  mu_check("Homo (new):", err);     
  h_Ry(&m2,0.3,&err);
  mu_check("Homo (Ry):", err);  
   
  cp = tm_copy(&m2,&err);
  mu_check("Homo (copy):", err);  
   
  h_inv(&m2,&err);
  mu_check("Homo (inv):", err);
  pr = h_new(&err);
  mu_check("Homo (new):", err); 
  tm_mul(&pr,&m2,&cp,&err);
  for(i = 0; i < 4; i++) {
    for(j = 0; j < 4; j++) {
      v = (i == j) ? 1 : 0;
      mu_assert("Homo (inv): wrong matrix", EQL(tm_get(&pr,i,j,0),v));       
    }
  }   
  tm_clear(&cp);
   
  cp = tm_copy(&m2,&err);
  mu_check("Homo (copy):", err);  
   
  h_T(&m2,&err);
  mu_check("Homo (T):", err);  
  for(i = 0; i < 4; i++) {
    for(j = 0; j < 4; j++) {
      mu_assert("Homo (T): wrong matrix", EQL(tm_get(&m2,i,j,0),tm_get(&cp,j,i,0))); 
    }
  }
   
  tm_mul(&pr,&m1,&m2,&err);
  mu_check("Homo (mul):", err);     
  h_mul(&m1,&m2,&err);
  mu_check("Homo (mul):", err); 
  for(i = 0; i < 4; i++) {
    for(j = 0; j < 4; j++) {
      mu_assert("Homo (mul): wrong matrix", EQL(tm_get(&m1,i,j,0),tm_get(&pr,i,j,0))); 
    }
  }
     
  tm_clear(&pr);
  tm_clear(&cp);
  tm_clear(&m1);
  tm_clear(&m2);
   
  return 0;
}

static char* test_vec()
{
  int err = 0,i;
  tmVal a1[9] = {1,3,4,5,7,2,8,6,1}, d, a2[] = {-10,6,-2};
  tMat v1, v2, pr; 
  
  v1 = vec_static(3,a1,&err);
  mu_check("Vec (static):", err);
  
  mu_assert("Vec (get): wrong value", EQL(vec_get(&v1,0,0),1));
  
  d = vec_norm2(&v1,&err);
  mu_check("Vec (norm2):", err);
  mu_assert("Vec (norm2):", EQL(d, a1[0]*a1[0]+a1[1]*a1[1]+a1[2]*a1[2]));
  
  v2 = tm_static(3,3,a1,&err);
  mu_check("Vec (static):", err);
  v2 = tm_col(&v2,1,&err);
  mu_check("Vec (col):", err);
  
  d = vec_dot(&v1,&v2,&err);
  mu_check("Vec (dot):", err);
  mu_assert("Vec (dot): wrong result", EQL(d,48));
  
  pr = tm_simp();
  vec_cross(&pr,&v1,&v2,&err);
  mu_check("Vec (cross):", err);
  
  for(i = 0; i < 3; i++) {
    mu_assert("Vec (cross): wrong product", EQL(vec_get(&pr,i,0),a2[i]));
  }
  
  vec_normalize(&v1,&err);
  mu_check("Vec (normalize):", err);
  mu_assert("Vec (normalize):", EQL(1.0, a1[0]*a1[0]+a1[1]*a1[1]+a1[2]*a1[2]));
    
  tm_clear(&pr);
  
  return 0;
}

static char* test_rank()
{
  int err = 0, rnk;
  tmVal a1[6] = {1,2,3,1,2,3}, a2[9] = {1,2,3,0,0,6,7,8,9};
  tMat m1;

  m1 = tm_static(2,3,a1,&err);
  mu_check("Rank (static):",err);

  rnk = tm_rank(&m1,&err);
  mu_check("Rank (rank):", err);
  mu_assert("Rank (rank): wrong result", rnk == 1);

  m1 = tm_static(3,3,a2,&err);
  mu_check("Rank (static):",err);

  rnk = tm_rank(&m1,&err);
  mu_check("Rank (rank):", err);
  mu_assert("Rank (rank): wrong result", rnk == 3);

  return 0;
}

/* one of minor values */
tmVal minor00(tMat *m, tmSize N, tmSize r, tmSize c, int* err)
{
  /* N is used as a matrix index in this example */
  return tm_get(m+N, r+1, c+1, err);
}

static char* test_make()
{
  int err = 0, i, j;
  tmVal e1,e2;
  tmVal a1[4] = {1,2,3,4}, a2[4] = {5,6,7,8};
  tMat m[2], c1,c2,c3;
  
  m[0] = tm_static(2,2,a1,0);
  m[1] = tm_static(2,2,a2,0);
  
  c1 = tm_concat(m,2,TM_HORIZONTAL,&err);
  mu_check("Make (concat h):", err);
    
  for(i = 0; i < 2; i++) {
    for(j = 0; j < 2; j++) {
      e1 = tm_get(m,i,j,&err);
      mu_check("Make (get):", err);
      e2 = tm_get(&c1,i,j,&err);
      mu_check("Make (get):", err);
      mu_assert("Make (concat): wrong horizontal concatenation", EQL(e1,e2));
      e1 = tm_get(m+1,i,j,&err);
      mu_check("Make (get):", err);
      e2 = tm_get(&c1,i,j+2,&err);
      mu_check("Make (get):", err);
      mu_assert("Make (concat): wrong horizontal concatenation", EQL(e1,e2));
    }
  }
  
  c1 = tm_concat(m,2,TM_VERTICAL,&err);
  mu_check("Make (concat v):", err);
  
  for(i = 0; i < 2; i++) {
    for(j = 0; j < 2; j++) {
      e1 = tm_get(m,i,j,&err);
      mu_check("Make (get):", err);
      e2 = tm_get(&c1,i,j,&err);
      mu_check("Make (get):", err);
      mu_assert("Make (concat): wrong horizontal concatenation", EQL(e1,e2));
      e1 = tm_get(m+1,i,j,&err);
      mu_check("Make (get):", err);
      e2 = tm_get(&c1,i+2,j,&err);
      mu_check("Make (get):", err);
      mu_assert("Make (concat): wrong horizontal concatenation", EQL(e1,e2));
    }
  }
  
  c3 = tm_make(m,0,m[0].rows-1,m[0].cols-1,minor00,&err);
  mu_check("Make (make):", err);
  mu_assert("Make (make): wrong value", EQL(tm_get(m,1,1,0),tm_get(&c3,0,0,0)));   
  
  tm_clear(&c1);
  tm_clear(&c2);
  tm_clear(&c3);
  
  return 0;
} 
