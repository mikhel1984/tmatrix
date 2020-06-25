/**
 * @file tmatrix_io.h
 * @author Stanislav Mikhel
 * @date 2020
 * @brief Functions for data input and output.
 */
#ifndef T_MATRIX_IO_H
#define T_MATRIX_IO_H 

#include "tmatrix.h"

#define tm_to_csv(mat,fname)    tm_to_file(mat,fname,',')

#define tm_from_csv(mat,fname)  tm_from_file(mat,fname,',')

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
/**
 * @brief Save table into file.
 * @param m matrix to save.
 * @param fname file name.
 * @parem sep separator character.
 * @return 1 in case of success.
 */
int tm_to_file(tMat* m, char* fname, char sep);
/**
 * @brief Initialize matrix from file 
 * @param dst matrix object.
 * @param fname file name.
 * @param sep separator character.
 * @return 1 in case of success.
 */
int tm_from_file(tMat* dst, char* fname, char sep);

#endif /* T_MATRIX_IO_H */
