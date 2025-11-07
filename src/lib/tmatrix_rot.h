/**
 * @file tmatrix_rot.h
 * @author Stanislav Mikhel
 * @date 2022
 * @brief Operations with quaternions and rotation matrices.
 * 
 * Rotation matrices must be 3x3. Quaternions can have arbitrary 
 * elements, but in the case of rotation unit quaternion are 
 * required.
 */
#ifndef T_MATRIX_ROT_H
#define T_MATRIX_ROT_H

#include "tmatrix.h"

/*=============== Rotation representation =================*/

/**
 * @brief Get matrix from Euler angles.
 *
 * Rotate X(roll)-Y(pitch)-Z(yaw) in fixed frame or 
 * Z(yaw)-Y(pitch)-X(roll) in current frames.
 * @param m matrix to change.
 * @param roll x rotation angle.
 * @param pitch y rotation angle.
 * @param yaw z rotation angle.
 * @param err error code.
 * @return 1 in case of success.
 */
int rot_rpy(tMat* m, tmVal roll, tmVal pitch, tmVal yaw, int* err);

/**
 * @brief Get Euler angles from the matrix.
 *
 * @param roll x rotation angle.
 * @param pitch y rotation angle.
 * @param yaw z rotation angle.
 * @param r rotation matrix.
 * @param err error code.
 * @return 1 in case of success.
 */
int rot_torpy(tmVal* roll, tmVal* pitch, tmVal* yaw, tMat* r, int *err);

/**
 * @brief Get matrix from axis-angle representation.
 *
 * @param m matrix to change.
 * @param k unit vector 3x1.
 * @param a angle of rotation.
 * @param err error code.
 * @return 1 in case of success.
 */
int rot_aa(tMat* m, tMat* k, tmVal a, int* err);

/**
 * @brief Get axin and angle from the matrix.
 * 
 * @param k axis vector to change.
 * @param a angle of rotation.
 * @param r rotation matrix.
 * @param err error code.
 * @return 1 in case of success.
 */
int rot_toaa(tMat* k, tmVal* a, tMat* r, int* err);

/**
 * @brief Get matrix from unit quaternion.
 *
 * @param m matrix to change.
 * @param q unit quaternion.
 * @param err error code.
 * @return 1 in case of success.
 */
int rot_qn(tMat* m, tQn* q, int* err);

/**
 * @brief Get quaternion from the matrix.
 *
 * @param q quaternion to change.
 * @param r rotation matrix.
 * @param err error code.
 * @return 1 in case of success.
 */
int rot_toqn(tQn* q, tMat* r, int* err);

/**
 * @brief Inversion of rotation matrix.
 *
 * @param dst matrix to store the result.
 * @param r rotation matrix.
 * @param err error code.
 * @return 1 in case of success.
 */
int rot_inv(tMat* dst, tMat* r, int* err);

/*============= Quaternions ==============*/

/**
 * @brief Get sum of two quaternions.
 *
 * @param q1 first quaternion.
 * @param q2 second quaternion.
 * @param err error code.
 * @return operation result.
 */
tQn qn_add(tQn* q1, tQn* q2, int* err);

/**
 * @brief Get difference of two quaterions.
 *
 * @param q1 first quaternion.
 * @param q2 second quaternion.
 * @param err error code.
 * @return operation result.
 */
tQn qn_sub(tQn* q1, tQn* q2, int* err);

/**
 * @brief Get product of two quaterions.
 *
 * @param q1 first quaternion.
 * @param q2 second quaternion.
 * @param err error code.
 * @return operation result.
 */
tQn qn_mul(tQn* q1, tQn* q2, int* err);

/** 
 * @brief Multiply quaternion with scalar value.
 *
 * @param q source quaternion.
 * @param k factor.
 * @param err error code.
 * @return scaled quaternion.
 */
tQn qn_scale(tQn* q, tmVal k, int* err);

/**
 * @brief Get quaternion conjugate.
 * 
 * @param q source quaternion.
 * @param err error code.
 * @return conjugation.
 */
tQn qn_conj(tQn* q, int* err);

/**
 * @brief Get quaternion inversion.
 *
 * Find such quaternion that q * inv(q) = 1
 * @param q source quaternion.
 * @param err error code.
 * @return inverted quaternion.
 */
tQn qn_inv(tQn* q, int* err);

/**
 * @brief In-place quaternion normalization.
 *
 * Obtain unit quaternion.
 * @param q quaternion.
 * @return 1 in case of success.
 */
int qn_normalize(tQn* q);

/**
 * @brief Find absolute value.
 * 
 * @param q quaternion.
 * @return absolute value.
 */
tmVal qn_abs(tQn* q);

/**
 * @briev Use SLERP to rotate unit quaternion.
 *
 * @param q1 initial unit quaternion.
 * @param q2 final unit quaternion.
 * @param t interpolation ratio (0 - 1).
 * @param err error code.
 * @return intermediate quaternion.
 */
tQn qn_slerp(tQn* q1, tQn* q2, tmVal t, int *err);

#endif /* T_MATRIX_ROT_H */
