#ifndef T_MATRIX_PRIVATE_H
#define T_MATRIX_PRIVATE_H

/* Matrix type */
#define TM_STATIC    0
#define TM_MAIN      1
#define TM_TRANSPOSE 2
#define TM_SUB       3

/* Error types */
enum TM_ERR {
     TM_ERR_WRONG_SIZE = 1,   /* Wrong matrix size definition  */
     TM_ERR_NO_MEMORY,        /* Can't allocate memory */
     TM_ERR_EMPTY_ARGS,       /* Some of pointers in argument list are empty */
     TM_ERR_NOT_MAIN,         /* Main matrix (or sometimes static) is expected */
     TM_ERR_DIFF_SIZE,        /* Src and dst have different size */
     TM_ERR_NOT_COMPAT,       /* Matrices are not compatible for current operation */
     TM_ERR_NOT_DEF,          /* Operation is not defined */
     TM_ERR_NO_SOLUTN,        /* Solution cannot be found */
     TM_ERR_NOT_HOMO,         /* Not a homogenous transform matrix */
     TM_ERR_NOT_VEC,          /* Not a vector */
     TM_ERR_TOTAL             /* Number of errors */   
};

#define IS_PRIM(X) ((X)->type == TM_MAIN || (X)->type == TM_STATIC)

tmVal* tm_at(tMat* m, tmSize r, tmSize c);

#endif /* T_MATRIX_PRIVATE_H */
