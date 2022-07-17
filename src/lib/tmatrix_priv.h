/**
 * @file tmatrix_priv.h
 * @author Stanislav Mikhel
 * @date 2020
 * @brief "Private" matrix parameters and methods.
 */
#ifndef T_MATRIX_PRIVATE_H
#define T_MATRIX_PRIVATE_H

/** 
 * @brief Matrix types. 
 */
#define TM_STATIC    0   /**< Matrix contains pointer to static memory. */
#define TM_ALLOC     1   /**< Matrix contains pointer to dynamic memory. */
#define TM_TRANSPOSE 2   /**< Matrix has link to another matrix data and define transposition access. */
#define TM_SUB       3   /**< Matrix has link to another matrix data and define sub matrix access. */
/**
 * @brief Check if the matrix has direct access to the data.
 * @param X pointer to matrix.
 * @return 1 if the matrix has static or dynamic types.
 */
#define IS_PRIM(X) ((X)->type == TM_ALLOC || (X)->type == TM_STATIC)

/** 
 * @brief Error types. 
 */
enum TM_ERR {
     TM_ERR_WRONG_SIZE = 1,   /**< Wrong matrix size definition. */
     TM_ERR_NO_MEMORY,        /**< Can't allocate memory. */
     TM_ERR_EMPTY_ARGS,       /**< Some of mandatory pointers in argument list are empty. */
     TM_ERR_NOT_MAIN,         /**< Main matrix (or sometimes static) is expected for current operation. */
     TM_ERR_DIFF_SIZE,        /**< Src and dst have different size. */
     TM_ERR_NOT_COMPAT,       /**< Matrices are not compatible for current operation. */
     TM_ERR_NOT_DEF,          /**< Operation is not defined. */
     TM_ERR_NO_SOLUTN,        /**< Solution cannot be found. */
     TM_ERR_NOT_HOMO,         /**< Not a homogeneous transform matrix. */
     TM_ERR_NOT_VEC,          /**< Not a vector. */
     TM_ERR_TOTAL             /**< Number of errors in list. */   
};

/**
 * @brief Macro for safety parameters.
 */
#ifdef TM_CHECK_ARGS
#define TM_ASSERT_ARGS(cond,var,lbl)  if(!(cond)) {var = TM_ERR_EMPTY_ARGS; goto lbl;}
#else
#define TM_ASSERT_ARGS(cond,var,lbl)  ((void)0)
#endif

#ifdef TM_CHECK_INDEX
#define TM_ASSERT_INDEX(cond,var,lbl) if(!(cond)) {var = TM_ERR_WRONG_SIZE; goto lbl;}
#else
#define TM_ASSERT_INDEX(cond,var,lbl)  ((void)0)
#endif

/**
 * @brief Fast memory access without additional checking.
 * @param m matrix object.
 * @param r row index.
 * @param c column index.
 * @return Pointer to matrix element.
 */
tmVal* tm_at(tMat* m, tmSize r, tmSize c);
/**
 * @brief Check matrix size, correct if need and possible.
 * @param m matrix object.
 * @param R desired number of rows.
 * @param C desired number of columns.
 * @param err error code.
 * @return 1 in case of success.
 */
int tm_relevant(tMat* m, tmSize R, tmSize C, int* err);

#endif /* T_MATRIX_PRIVATE_H */
