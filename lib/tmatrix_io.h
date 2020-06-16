/**
 * @file tmatrix_io.h
 * @author Stanislav Mikhel
 * @date 2020
 * @brief Functions for data input and output.
 */
#ifndef T_MATRIX_IO_H
#define T_MATRIX_IO_H 

#include "tmatrix.h"

/** 
 * @brief Simple matrix visualization.
 * @param m matrix object. 
 */
void tm_print(tMat *m);
/**
 * @brief Error description.
 * @param code error value.
 * @return Description string.
 */
const char* tm_error(int code);

#endif /* T_MATRIX_IO_H */
