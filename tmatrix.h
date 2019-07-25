/**
 * @file tmatrix.h
 * @author Stanislav Mikhel
 * @date 2019
 * @brief Definitions of the main structures and functions for matrices.
 * 
 * This library is written in pure C and focused on coordinate transformations
 * in robotics applications. The header contains main data types and functions 
 * for the matrix operations in general.  
 */

#ifndef T_MATRIX_H
#define T_MATRIX_H

/** 
 * @brief Empty matrix initialization 
 */
#define NULL_TMATRIX {0,0,0,0,0}
/** 
 * @brief Return empty matrix with dynamically allocated memory.
 */
#define tm_simp()          tm_new(0,0,0)
/**
 * @brief Get column.
 * @param mat pointer to matrix object.
 * @param c column number.
 * @param err error pointer.
 */
#define tm_col(mat,c,err)  tm_block(mat, 0, c, (mat)->rows, 1, err)
/**
 * @brief Get row.
 * @param mat pointer to matrix object.
 * @param r row number.
 * @param err error pointer.
 */
#define tm_row(mat,r,err)  tm_block(mat, r, 0, 1, (mat)->cols, err)

/** 
 * @brief Internal data types.
 */
typedef unsigned char tmSize;  /**< Number of rows/columns (assumed < 256).  */
typedef double tmVal;          /**< Type of the matrix elements. */

/** 
 * @brief Matrix object. 
 *
 * Matrix could contain allocated memory or include pointer to another object
 * with definition of access rules (for transposed matrix of sub matrix).
 */
typedef struct tMat_ {
   tmVal* data;                /**< Pointer to data array. */
   tmSize rows;                /**< Number of rows.        */
   tmSize cols;                /**< Number of columns.     */
   tmSize width;               /**< Parameter is used for index evaluation. */
   tmSize type;                /**< Type of memory / element access. */
}  tMat;

/** 
 * @brief Create matrix of given size.
 *
 * Method allocate memory and initialize it with zeros. 
 * @param r number of rows.
 * @param c number of columns.
 * @param err error code.
 * @return New matrix. 
 * @note Free memory with @a tm_clear.
 */
tMat tm_new(tmSize r, tmSize c, int* err);
/** 
 * @brief Create matrix from static array.
 *
 * Method pointer into a matrix structure without modifications.
 * @param r number of rows.
 * @param c number of columns.
 * @param dat data array.
 * @param err error code.
 * @return New matrix. 
 * @note Be shure that the data array is more of equal to r * c.
 */
tMat tm_static(tmSize r, tmSize c, tmVal dat[], int* err);
/** 
 * @brief Create identity matrix.
 *
 * Method allocate memory, initialize it with zeros and sets diagonal elements
 * equal to 1.
 * @param r number of rows.
 * @param c number of columns.
 * @param err error code.
 * @return New matrix. 
 * @note Free memory with @a tm_clear.
 */
tMat tm_eye(tmSize r, tmSize c, int *err);
/** 
 * @brief Clear allocated memory (if need).
 * @param m pointer to matrix.
 */
void tm_clear(tMat* m);
/**  
 * @brief Copy numbers from array into the matrix.
 * 
 * It is assumed that the array represents data row by row.
 * @param dst destination matrix.
 * @param src array of numbers.
 * @param err error code.
 * @return 1 in case of success. 
 */
int tm_init(tMat* dst, tmVal src[], int* err);
/** 
 * @brief Get matrix element.
 * @param m matrix object.
 * @param r row number.
 * @param c column number.
 * @param err error code.
 * @return Element value.
 */
tmVal tm_get(tMat* m, tmSize r, tmSize c, int* err);
/** 
 * @brief Set matrix element.
 * @param m matrix object.
 * @param r row number.
 * @param c column number.
 * @param v new value.
 * @param err error code.
 */
void tm_set(tMat* m, tmSize r, tmSize c, tmVal v, int* err);
/** 
 * @brief Get number of rows.
 * @param m matrix object.
 * @return Number of rows.
 */
tmSize tm_rows(tMat* m);
/** 
 * @brief Get number of columns.
 * @param m matrix object.
 * @return Number of columns.
 */
tmSize tm_cols(tMat* m);
/** 
 * @brief Create deep copy of the matrix.
 * @param src source matrix.
 * @param err - error code.
 * @return Matrix copy.
 * @note Free memory with @a tm_clear.
 */
tMat tm_copy(tMat* src, int *err);
/** 
 * @brief Create transposed copy of the matrix.
 *
 * Method doesn't allocate new memory, it is just a reference to the original matrix.
 * @param src source matrix.
 * @param err error code.
 * @return Transposed matrix. 
 */
tMat tm_T(tMat* src, int *err);
/** 
 * @brief Create sub matrix.
 * 
 * Method doesn't allocate new memory, is is just a reference to the original matrix.
 * @param src source matrix.
 * @param r0 begin row number.
 * @param c0 begin column number.
 * @param Nr number of rows.
 * @param Nc number of columns.
 * @param err error code.
 * @return Sub matrix.
 */
tMat tm_block(tMat* src, tmSize r0, tmSize c0, tmSize Nr, tmSize Nc, int *err);
/** 
 * @brief Copy values from source to destination.
 *
 * Matrix size must be equal.
 * @param dst destination matrix.
 * @param src source matrix.
 * @param err error code.
 * @return 1 in case of success. 
 */
int tm_insert(tMat *dst, tMat *src, int* err);
/** 
 * @brief Simple matrix visualization.
 * @param m matrix object. 
 */
void tm_print(tMat *m);

/*============== Arithmetic methods ==============*/

/** 
 * @brief Add two matrices.
 * 
 * Save result to dst: <i> dst += m </i>.
 * @param dst matrix to change.
 * @param m second matrix.
 * @param err error code.
 * @return 1 in case of success. 
 */
int tm_add(tMat *dst, tMat* m, int* err);
/** 
 * @brief Subtract two matrices.
 * 
 * Save result to dst: <i> dst -= m </i>.
 * @param dst matrix to change.
 * @param m second matrix.
 * @param err error code.
 * @return 1 in case of success. 
 */
int tm_sub(tMat *dst, tMat* m, int* err);
/** 
 * @brief Multiply matrix to scalar value.
 * 
 * Save result to dst: <i> dst *= k </i>.
 * @param dst matrix to change.
 * @param k multiplier.
 * @param err error code.
 * @return 1 in case of success. 
 */
int tm_scale(tMat *dst, tmVal k, int* err);
/** 
 * @brief Multiply matrices.
 * 
 * Save result to dst: <i> dst = a * b </i>.
 * Destination matrix must has equal size or be dynamic.
 * @param dst matrix for result.
 * @param a first matrix.
 * @param b second matrix.
 * @param err error code.
 * @return 1 in case of success. 
 */
int tm_mul(tMat* dst, tMat *a, tMat *b, int *err);
/** 
 * @brief Find determinant.
 * 
 * Matrix must be square.
 * @param m matrix object.
 * @param err error code.
 * @return Determinant value. 
 */
tmVal tm_det(tMat *m, int *err);
/** 
 * @brief Find inverted matrix.
 * 
 * Source matrix must be equal.
 * @param m matrix object.
 * @param err error code.
 * @return Inverted matrix.
 * @note Free memory with @a tm_clear.
 */
tMat tm_inv(tMat *m, int *err);
/** 
 * @brief Find pseudo-inverted matrix.
 * 
 * @param m matrix object.
 * @param err error code.
 * @return Pseudo-inverted matrix.
 * @note Free memory with @a tm_clear.
 */
tMat tm_pinv(tMat *m, int *err);
/**
 * @brief Error description.
 * @param code error value.
 * @return Description string.
 */
const char* tm_error(int code);
/**
 * @brief Find rank of the matrix.
 * @param m matrix object.
 * @param err error code.
 * @return Rank value.
 */
int tm_rank(tMat *m, int *err);

#endif /* T_MATRIX_H */
