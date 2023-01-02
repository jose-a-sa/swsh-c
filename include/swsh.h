#ifndef SWSH_H
#define SWSH_H

#include <stdbool.h>
#include <stdint.h>
#include "const.h"

bool possibly_zero(ldouble_t);
bool possibly_same(ldouble_t, ldouble_t);

uint64_t factorial(int);
int64_t binomial(int, int);

ldouble_t powl_t(ldouble_t, ldouble_t);
ldouble_t sgn_t(int);

ldouble_t Yslm(int, int, int, ldouble_t);
ldouble_t YslmPrime(int, int, int, ldouble_t);

ldouble_t lim_Yslm(int, int, int, int);
ldouble_t lim_YslmPrime(int, int, int, int);
ldouble_t limratio_prime_Yslm(int, int, int, ldouble_t, ldouble_t, int);

ldouble_t g(int, vector_t, ldouble_t, int, int, int, ldouble_t);
ldouble_t Dg(int, int, vector_t, ldouble_t, int, int, int, ldouble_t);

ldouble_t Ek(int, int, vector_t, matrix_t, int, int, int, ldouble_t);
ldouble_t Sk(int, int, int, vector_t, matrix_t, int, int, int, ldouble_t);

void generate_smatrix_block(matrix_t, int, vector_t, matrix_t, int, int, int, ldouble_t);
void block_gauss_pivot(matrix_t);
void block_reduce(matrix_t, matrix_t, size_t);
void block_back_substituition(matrix_t *, size_t, size_t);
void generate_initial_Y(matrix_t, vector_t, int, int, int, ldouble_t *);
void desolve_relaxation(vector_t, matrix_t, int, int, int, ldouble_t, ldouble_t *, u_int32_t);

#endif // SWSH_H
