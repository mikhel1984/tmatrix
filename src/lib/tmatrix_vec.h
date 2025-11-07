/**
 * @file tmatrix_vec.h
 * @author Stanislav Mikhel
 * @date 2020
 * @brief Definitions of vector specific methods.
 *
 * Vector here is a matrix with number of rows or columns equal to 1.
 * Provides all basic matrix operations.
 */
#ifndef T_MATRIX_VECTOR_H
#define T_MATRIX_VECTOR_H

#include "tmatrix.h"

/**
 * @brief Create vector with dynamic memory.
 * @param n vector length.
 * @param err error code.
 * @return New vector.
 * @note Free memory with @a tm_clear.
 */
#define vec_new(n,err)        tm_new(n,1,err)

/**
 * @brief Create vector with static memory.
 * @param n vector length.
 * @param dat data array.
 * @param err error code.
 * @note Data length must be equal to @a n.
 */
#define vec_static(n,dat,err) tm_static(n,1,dat,err)

/**
 * @brief Get vector length.
 * @param m matrix object.
 * @return Vector length or 0 if the object is not a vector.
 */
tmSize vec_len(tMat *m);

/**
 * @brief Get vector element.
 * @param m matrix object.
 * @param k index.
 * @param err error code.
 * @return Vector element.
 */
tmVal vec_get(tMat *m, tmSize k, int *err);

/**
 * @brief Update vector element.
 * @param m matrix object.
 * @param k index.
 * @param v new value.
 * @param err error code.
 */
void vec_set(tMat *m, tmSize k, tmVal v, int *err);

/**
 * @brief Dot product of two vectors.
 * @param a first vector.
 * @param b second vector.
 * @param err error code.
 * @return Product value.
 */
tmVal vec_dot(tMat *a, tMat *b, int *err);

/**
 * @brief Cross product of two vectors.
 * 
 * Save result to res. Vectors must be of the length 3.
 * @param res matrix for result (with length equal to 3 or dynamic memory).
 * @param a first vector.
 * @param b second vector.
 * @param err error code.
 * @return 1 in case of success.
 */
int vec_cross(tMat *res, tMat *a, tMat *b, int *err);

/**
 * @brief Find square of Eucledian norm.
 *
 * The function calculates sum of squares for vector elements.
 * @param m matrix object.
 * @param err error code.
 * @return Square norm.
 */
tmVal vec_norm2(tMat *m, int *err);

/**
 * @brief Normalize vector.
 * 
 * Obtaine unit vector if its norm is not 0.
 * @param m matrix object.
 * @param err error code.
 * @return 1 in case of success.
 */
int vec_normalize(tMat *m, int *err);

#endif /* T_MATRIX_VECTOR_H */
