/**
 * @file tmatrix_homo.h
 * @author Stanislav Mikhel
 * @date 2019
 * @brief Definition of methods for homogeneous transformations.
 *
 * Homogeneous transformations assumes manipulations with matrices 4x4. 
 * Provides all basic matrix operations.
 */
#ifndef T_MATRIX_HOMO_H
#define T_MATRIX_HOMO_H

#include "tmatrix.h"

/**
 * @brief Homogeneous matrix size. 
 */
#define TM_HOMO_SIZE 16
/**
 * @brief Get matrix for X translation.
 * @param dst pointer to matrix.
 * @param x displacement.
 * @param err error code.
 * @return 1 in case of success.
 */
#define h_Tx(dst,x,err)  h_Txyz(dst,x,0,0,err)
/**
 * @brief Get matrix for Y translation.
 * @param dst pointer to matrix.
 * @param y displacement.
 * @param err error code.
 * @return 1 in case of success.
 */
#define h_Ty(dst,y,err)  h_Txyz(dst,0,y,0,err)
/**
 * @brief Get matrix for Z translation.
 * @param dst pointer to matrix.
 * @param z displacement.
 * @param err error code.
 * @return 1 in case of success.
 */
#define h_Tz(dst,z,err)  h_Txyz(dst,0,0,z,err)
/**
 * @brief Create homogeneous matrix with dynamic memory.
 * @param err error code.
 * @return New matrix.
 * @note Free memory with @a tm_clear.
 */
#define h_new(err)       tm_new(4,4,err)
/**
 * @brief Create homogeneous matrix with static memory.
 * @param err error code.
 * @return New matrix.
 */
#define h_static(d,err)  tm_static(4,4,d,err)
/**
 * @brief Get matrix for combined X-Y-Z translations.
 * @param dst destination matrix.
 * @param x displacement in X.
 * @param y displacement in Y.
 * @param z displacement in Z.
 * @param err error code.
 * @return 1 in case of success.
 */
int h_Txyz(tMat* dst, tmVal x, tmVal y, tmVal z, int *err);
/**
 * @brief Get matrix for Z rotation.
 * @param dst destination matrix.
 * @param a angle (rad).
 * @param err error code.
 * @return 1 in case of success.
 */
int h_Rz(tMat* dst, tmVal a, int *err);
/**
 * @brief Get matrix for Y rotation.
 * @param dst destination matrix.
 * @param b angle (rad).
 * @param err error code.
 * @return 1 in case of success.
 */
int h_Ry(tMat* dst, tmVal b, int *err);
/**
 * @brief Get matrix for X rotation.
 * @param dst destination matrix.
 * @param v angle (rad).
 * @param err error code.
 * @return 1 in case of success.
 */
int h_Rx(tMat* dst, tmVal v, int *err);
/** 
 * @brief Find matrix for Denavit-Hartenberg parameters.
 * @param dst destination matrix.
 * @param a length of common normal.
 * @param alpha angle about common normal.
 * @param d offset along previous Z to common normal.
 * @param theta angle about previous Z.
 * @param err error code.
 * @return 1 in case of success.
 */
int h_DH(tMat* dst, tmVal a, tmVal alpha, tmVal d, tmVal theta, int *err);
/** 
 * @brief Multiply matrices.
 * 
 * Save result to dst: <i> dst *= m </i>.
 * @param dst matrix for result.
 * @param m second matrix.
 * @param err error code.
 * @return 1 in case of success. 
 */
int h_mul(tMat *dst, tMat *m, int* err);
/** 
 * @brief Find inverted matrix.
 * 
 * Save result to dst.
 * @param dst source (and destination) matrix.
 * @param err error code.
 * @return 1 in case of success.
 */
int h_inv(tMat *dst, int *err);
/** 
 * @brief Find transposed matrix.
 * 
 * Save result to dst.
 * @param dst source (and destination) matrix.
 * @param err error code.
 * @return 1 in case of success.
 */
int h_T(tMat *dst, int *err);
/** 
 * @brief Initialize identity matrix.
 * 
 * @param dst source (and destination) matrix.
 * @param err error code.
 * @return 1 in case of success.
 */
int h_eye(tMat *dst, int *err);

#endif 
