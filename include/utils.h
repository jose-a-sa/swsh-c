#ifndef UTILS_H
#define UTILS_H

#include "const.h"
#include <stdio.h>

void print_smatrix_block(matrix_t, size_t);
void print_smatrix_block_full(matrix_t, size_t);
void output_vector_h(FILE *, vector_t);
void output_matrix(FILE *, matrix_t);
void output_block_matrix(FILE *, matrix_t *, size_t, size_t);
void output_block_matrix_mathematica(FILE *, matrix_t *, size_t, size_t);

#endif // UTILS_H
