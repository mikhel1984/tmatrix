/**
 * @file tmatrix_transform.h
 * @author Stanislav Mikhel
 * @date 2025
 * @brief Matrix decompositions.
 */
#ifndef T_MATRIX_TRANSFORM_H
#define T_MATRIX_TRANSFORM_H

#include "tmatrix.h"

/**
 * @brief Cholesky decomposition.
 *
 * Find such matrix C for M that M = C*C^T.
 * Error code is TM_ERR_NOT_POS_DEF when matrix is not positive definite.
 * @param dst found lower left matrix.
 * @param m source matrix.
 * @param err error code.
 */
void tf_chol(tMat* dst, tMat* m, int* err);

/**
 * @brief LU decomposition with permutations.
 *
 * Find such matrices L, U, P for the given matrix M that 
 * L*U = P*M.
 * @param L lower left triangle matrix.
 * @param U upper right triangle matrix.
 * @param P permutation matrix.
 * @param m source matrix.
 * @param err error code.
 */
void tf_lup(tMat* L, tMat* U, tMat* P, tMat* m, int* err);

/**
 * @brief LU decomposition.
 *
 * Find such matrices L, U for the given matrix M that L*U = M.
 * @param L lower left triangle matrix.
 * @param U upper right triangle matrix.
 * @param m source matrix.
 * @param err error code.
 */
void tf_lu(tMat* L, tMat* U, tMat* m, int* err);

void tf_qr(tMat* Q, tMat* R, tMat* m, int* err);

#endif 
