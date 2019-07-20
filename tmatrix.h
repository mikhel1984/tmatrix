#ifndef T_MATRIX_H
#define T_MATRIX_H

#define NULL_TMATRIX {0,0,0,0,0}

#define tm_simp()          tm_new(0,0,0)
#define tm_col(mat,c,err)  tm_block(mat, 0, c, (mat)->rows, 1, err)
#define tm_row(mat,r,err)  tm_block(mat, r, 0, 1, (mat)->cols, err)

/* Data types */
typedef unsigned char tmSize;
typedef double tmVal;

/* Matrix object */
typedef struct tMat_ {
   tmVal* data;   
   tmSize rows;
   tmSize cols;         
   tmSize width;
   tmSize type;

}  tMat;

/** Create zero matrix of given size 
    @fn tm_new
    @param r - number of rows
    @param c - number of columns
    @param err - error code
    @return Matrix structure */
tMat tm_new(tmSize r, tmSize c, int* err);

/** Create matrix using static array of values
    @fn tm_static
    @param r - number of rows
    @param c - number of columns
    @param dat - array of values
    @param err - error code
    @return Matrix structure */
tMat tm_static(tmSize r, tmSize c, tmVal dat[], int* err);

/** Create identity matrix of given size 
    @fn tm_new
    @param r - number of rows
    @param c - number of columns
    @param err - error code
    @return Matrix structure */
tMat tm_eye(tmSize r, tmSize c, int *err);

/** Clear allocated memory (only for main type) 
    @fn tm_clear
    @param m - pointer to matrix */
void tm_clear(tMat* m);

/** Copy matrix content from array 
    @fn tm_init
    @param dst - destination matrix
    @param src - array of numbers
    @param err - error code
    @return 1 in case of success */
int tm_init(tMat* dst, tmVal src[], int* err);

/** Get matrix element 
    @fn tm_get
    @param m - matrix object
    @param r - row number
    @param c - column number
    @param err - error code
    @return element value */
tmVal tm_get(tMat* m, tmSize r, tmSize c, int* err);

/** Set matrix element 
    @fn tm_set
    @param m - matrix object
    @param r - row number
    @param c - column number
    @param v - new value
    @param err - error code */
void tm_set(tMat* m, tmSize r, tmSize c, tmVal v, int* err);

/** Get number of rows 
    @fn tm_rows
    @param m - matrix object
    @return number of rows */
tmSize tm_rows(tMat* m);

/** Get number of columns 
    @fn tm_cols
    @param m - matrix object
    @return number of columns */
tmSize tm_cols(tMat* m);

/** Create deep copy of the matrix 
    @fn tm_copy
    @param src - source matrix
    @param err - error code
    @return matrix copy */
tMat tm_copy(tMat* src, int *err);

/** Create 'interface' equal to matrix transposition
    @fn tm_T
    @param src - source matrix
    @param err - error code
    @return transposed matrix */
tMat tm_T(tMat* src, int *err);

/** Create 'interface' equal to matrix part
    @fn tm_block
    @param src - source matrix
    @param r0 - begin row number
    @param c0 - begin column number
    @param Nr - number of rows
    @param Nc - number of columns
    @param err - error code
    @return submatrix */
tMat tm_block(tMat* src, tmSize r0, tmSize c0, tmSize Nr, tmSize Nc, int *err);

/** Copy values from source to destination 
    if the size is equal
    @fn tm_insert
    @param dst - destination matrix
    @param src - source matrix
    @param err - error code
    @return 1 in case of success */
int tm_insert(tMat *dst, tMat *src, int* err);

/** Simple matrix visualization
    @fn tm_print
    @param m - matrix object */
void tm_print(tMat *m);

/******* Arithmetic module ********/

/** Get sum of two matrices, save result to dst:
    dst += m
    @fn tm_add
    @param dst - matrix to change (main or static)
    @param m - second matrix
    @param err - error code
    @return 1 in case of success */
int tm_add(tMat *dst, tMat* m, int* err);

/** Get difference of two matrices, save result to dst:
    dst -= m
    @fn tm_sub
    @param dst - matrix to change (main or static)
    @param m - second matrix
    @param err - error code
    @return 1 in case of success */
int tm_sub(tMat *dst, tMat* m, int* err);

/** Multiply matrix to constant value
    @fn tm_scale
    @param dst - matrix to change (main or static)
    @param k - multiplier
    @param err - error code
    @return 1 in case of success */
int tm_scale(tMat *dst, tmVal k, int* err);

/** Get product of two matrices, save result to dst
    @fn tm_mul
    @param dst - matrix for result
    @param a - first matrix
    @param b - second matrix
    @param err - error code
    @return 1 in case of error */
int tm_mul(tMat* dst, tMat *a, tMat *b, int *err);

/** Find determinant of a square matrix
    @fn tm_det
    @param m - matrix object
    @param err - error code
    @return determinant value */
tmVal tm_det(tMat *m, int *err);

/** Get inversion of a square matrix
    @fn tm_inv
    @param m - matrix object
    @param err - error code
    @return inverted matrix */
tMat tm_inv(tMat *m, int *err);

/** Get pseudo-inverse of an arbitrary matrix
    @fn tm_pinv
    @param m - matrix object
    @param err - error code 
    @return pseudo-invertion of the matrix */
tMat tm_pinv(tMat *m, int *err);

const char* tm_error(int code);


#endif /* T_MATRIX_H */
